##############################################################################
# Functions to model and forecast COVID infections in the US using an SIR 
# model and positive test data at the county level.
##############################################################################
library(tidyverse)  # most functions require dplyr, stringr, tidyr, etc

# fit_SIR() 
#-----------------------------------------------------------------------------#
# Estimate S,I,R values from positive test data.
#
# Inputs:
#   data: data frame with columns date, cases, population, FIPS
#   p.test: probability an infected individual has positive test
#   edf: effective degrees of freedom for case count smoothing
#   deg: smoothing with B-splines of degree=deg
# Output:
#   adds columns S, I, R to the data. These are the estimated counts of the 
#   population in each compartment. 
#-----------------------------------------------------------------------------#
fit_SIR <- function(data, p.test=1/10, edf=8, deg=3) {
  data %>% arrange(date) %>% 
    transmute(FIPS, date, cases, N=population,
              C = smooth_counts(cases, edf=edf, deg=deg),
              I = C/p.test) %>% 
    mutate(
      r = impute_r(I, pars=7), 
      R = pmin(N, cumsum(r)),  # force R in [0,N]
      S = N-I-R
    ) %>% 
    select(-r)
}


# impute_r() 
#-----------------------------------------------------------------------------#
# Estimate the number newly recovered time series
# basic aoristic type method; number recovered at time t a known function.
#  of the number infected at time s (s<t). It is a type of kernel smoothing/
#  or convolution. Uses poisson pmf for days to recovery. 
# Inputs:
#   I: time series of infected counts
#   pars: expected time from testing to recovery (in days)
# Notes:
#   - could use more flexible time to recovery (e.g., nbinom)
#-----------------------------------------------------------------------------#
impute_r <- function(I, pars=7) {
  wt = dpois(0:50, lambda=pars-1) 
  wt = c(0, wt[wt>1E-9])  # don't allow same day recovery (shifted poisson)
  wt = wt/sum(wt)         # re-sum to 1 
  n = length(wt)
  pad = rep(0, n)
  stats::filter(c(pad, I), filter=wt, method="conv", sides=1)[-(1:n)] 
}


# forecast_SIR() 
#-----------------------------------------------------------------------------#
# Estimate parameters and forecast SIR model using positive case count data
#-----------------------------------------------------------------------------#
forecast_SIR <- function(data, k=30, p.test=1/10, k.beta=21, gamma=1/14, edf=8, deg=3) {
  
  #- initial estimate of S, I, R from positive case counts
  data_SIR = data %>% arrange(date) %>% 
    fit_SIR(p.test=p.test, edf=edf, deg=deg)
  
  #- estimate SIR transmission parameter (beta) using recent data
  beta = data_SIR %>% 
    tail(k.beta) %>% 
    fit_beta(gamma=gamma)
  
  #- use fitted SIR to make forecasts
  data_SIR %>%
    tail(1) %>% {SIR(beta=beta, gamma=gamma, N=.$N, I=.$I, S=.$S, R=.$R, nT = k)} %>% 
    mutate(N = round(S + I + R), 
           date = max(data_SIR$date) + time - 1, 
           C = I*p.test) %>% 
    filter(time > 1) %>% # remove duplicate observation
    select(-time) %>% 
    bind_rows(data_SIR) %>% arrange(date) %>% 
    mutate(FIPS = FIPS[1]) %>% 
    mutate(p.inf = I/N, p.rec=R/N, p.sus=S/N) 
}  




# fit_beta() 
#-----------------------------------------------------------------------------#
# Estimate beta (transmission parameter) in SIR model
#
# Find that beta that minimizes the RMSE between "observed" I_t and estimated I_t(beta).
#   - The "observed" I_t is actually estimated from the observed positive test 
#     counts.
#   - Uses a penalized RMSE to encourage small value of beta in case of sparse data.
#     The `lambda` parameter controls the penalty
# Inputs:
#   data_SIR: data.frame with columns: date, S, I, R
#   gamma: the gamma parameter of SIR model. gamma = Pr(recovery)
#   lambda: strength of penalty. Large values force hat(beta) toward 0
# Outputs:
#   estimated beta value
# Notes:
#   - there are lots of ways to improve this, but its a start. 
#-----------------------------------------------------------------------------#
fit_beta <- function(data_SIR, gamma=1/14, lambda=.1) {
  #- tidy up data
  nT = nrow(data_SIR)                   # number of days of data
  data_SIR = data_SIR %>% arrange(date) # ensure arranged by date
  initial = data_SIR %>% slice(1)       # starting values for SIR model
  
  #- optimization function (RMSE with small penalty)
  opt_fun <- function(beta) {
    hat_I = SIR(beta=beta, gamma=gamma, N=initial$N, I=initial$I, S=initial$S, 
                R=initial$R, nT = nT)$I
    sqrt( mean( (data_SIR$I - hat_I)^2 )) + lambda*beta # add penalty to discourage large betas
  }
  
  #- optimize to estimate beta
  optimize(f=opt_fun, interval=c(1E-5, .15))$minimum
}



# SIR()
#-----------------------------------------------------------------------------#
# Simulate infection dynamics with S-I-R model (discrete time)
# 
# Inputs:
#   beta: transmission parameter; the expected number of people an infectious 
#         person comes into contact with each day. Can be a vector of length nT
#         to allow time varying values.
#   gamma: recovery parameter; an infected person recovers on a given day with
#           probability gamma
#   N: population size
#   S,I,R: initial number of susceptible (S), infectious (I), and 
#           recovered/removed (R)
#   nT: number of time periods to simulate
#   format: ('wide' or 'long'). Specifies if a long or wide table be returned. 
#   force: ('linear' or 'exponential'). The Force of Infection specification. 
#          If "linear", then probability of moving from S to I is beta*I/N. 
#          If "exponential", then its 1-exp(-beta*I/N).
#   method: ('deterministic' or 'stochastic'). Specifies if transitions occur
#           according to expected value or probability. 
#
# Outputs:
#   tibble (data.frame) of counts in each S-I-R compartment by day
#-----------------------------------------------------------------------------#
SIR <- function(beta, gamma, N=1E5, I=1, S=N-I, R=0, nT = 30, 
                format=c("wide", "long"), 
                force=c("linear", "exponential"),   # force of infection form 
                method=c("deterministic", "stochastic")) {
  format = match.arg(format)
  force = match.arg(force)
  method = match.arg(method)
  if(gamma<=0 | gamma>=1) stop("gamma must be between 0 and 1")
  if(min(beta)<0) stop("beta must be > 0")
  #if(min(beta, gamma)<0 | max(beta, gamma)>1) stop("beta and gamma must be between 0 and 1")
  
  if(length(beta) == 1) beta = rep(beta, nT)
  else if(length(beta) != nT) stop("beta must be of length 1 or nT")
  
  #- Initialize vectors
  #initial = c(S=S, I=I, R=R)
  # if(sum(initial != N)) stop('S + I + R must = N')
  initial = c(S = N-I-R, I=I, R=R)
  S = I = R = rep(0, nT)
  S[1] = initial[1]
  I[1] = initial[2]
  R[1] = initial[3]
  
  for(t in 2:nT) {
    if(force == "linear"){
      p.s2i =  beta[t-1] * I[t-1] / N  
      p.i2r = gamma 
    } else{    #force == "exponential"
      p.s2i = (1-exp(-beta[t-1]*I[t-1]/N))
      p.i2r = (1-exp(-gamma))
    }
    if(method == "deterministic") {
      s2i = p.s2i * S[t-1]
      i2r = p.i2r * I[t-1]
    } else {   # method == "stochastic"
      s2i = rbinom(n=1, size=S[t-1], prob=p.s2i)
      i2r = rbinom(n=1, size=I[t-1], prob=p.i2r)
    }
    
    I[t] = I[t-1] + s2i - i2r
    R[t] = R[t-1] + i2r
    S[t] = S[t-1] - s2i
  }
  out = tibble(time=1:nT, S, I, R)
  if(format == "long") out = out %>% gather(state, count, -time) %>% 
    mutate(state = factor(state, levels=c("S", "I", "R")))
  return(out)
}



# smooth_counts()
#-----------------------------------------------------------------------------#
# Smooth count data using mgcv::gam() function
# Notes: ignoring ... 
# edf is now estimated in a data driven manner using mcgv
#-----------------------------------------------------------------------------#
smooth_counts <- function(x, ...) {
  #- if mostly 0's, return a constant fit
  if(n_distinct(x) < 7) {
    yhat = rep(mean(x), length(x))
    attr(yhat, "edf") = 1
    return(yhat)
  }
  
  t = seq_along(x)
  fit = mgcv::gam(x ~ s(t), family="poisson")
  yhat = fit$fitted.values
  attr(yhat, "edf") = sum(fit$edf)
  return(yhat)
}
  


# 
# # smooth_counts()
# #-----------------------------------------------------------------------------#
# # Smooth count data
# # Notes:
# #   - if number of distinct values is less than edf+2, then fit a constant
# #   - better, perhaps, to reduce edf or deg to accommodate
# #-----------------------------------------------------------------------------#
# smooth_counts <- function(x, edf=5, phi=NULL, deg=3) {
#   n = length(x)
#   
#   #- if mostly 0's, return a constant fit
#   if(n_distinct(x) < edf + 2) {
#     yhat = rep(mean(x), n)
#     attr(yhat, "edf") = 1
#     return(yhat)
#   }
#   penMat = crossprod(diff(diag(n), diff=deg))
#   if(is.null(phi)) { # edf specified
#     if(is.null(edf)) stop("phi or edf must be specified")
#     if(edf < deg) stop("edf must be > deg")
#     tol.edf = .1
#     if(edf < deg + tol.edf/2) edf = deg + tol.edf/2
#     int = c(1E1, 1E5)
#     opt_fun <- function(phi){
#       tryCatch({
#         est_df = penfit(x, penMat = phi*n*penMat, show.warnings = FALSE)$edf
#         edf - est_df
#       }, error = function(e) -1E9
#       ) 
#     }
#     opt = uniroot(opt_fun, interval=int, tol=tol.edf, extendInt = 'upX')  
#     phi = opt$root 
#     fit = penfit(x, penMat = phi*n*penMat)
#     if(abs(fit$edf-edf)>tol.edf){ 
#       warning(paste("edf of model is more than", tol.edf, "from target edf."))
#     } 
#   } else { # phi specified
#     if(is.null(phi)) stop("phi or edf must be specified")
#     fit = penfit(x, penMat = phi*n*penMat)
#   }
#   yhat = fit$fit
#   attr(yhat, "edf") = fit$edf
#   attr(yhat, "phi") = phi
#   return(yhat)
# }
# 
# #-- penfit
# #------------------------------------------------------------------------------#
# # Estimate intensity with penalized poisson-like log-likelihood:
# # sum_i {n_i beta_i - m_i exp(beta_i)} - t(beta) %*% (penMat/2) %*% beta
# #
# # the intensity lambda_i (expected number infected) at time i is exp(beta_i)
# #
# # Inputs:
# #   n, m: vector of length t
# #   penMat: t x t penalty matrix ( phi*t(D)%*%D )
# #   beta.start: starting values for beta
# #   maxit: maximum number of iterations
# #   tol: tolerance. Stop when sum(abs(b-b.old)) < tol
# #
# # Output:
# #   list with elements: 
# #     - fit: estimated lambda
# #     - beta: estimated beta = log(lambda)
# #     - w: weight vector at convergence
# #     - edf: estimated effective degrees of freedom
# #     - niters: number of iterations to converge
# # 
# # TODO:
# #   - make faster
# #   - compute edf optional
# #   - add Std Errors
# #   - make convergence warning optional
# #------------------------------------------------------------------------------#
# penfit <- function(n, m=1, penMat, beta.start=NULL, maxit=200, tol=1e-5, out.edf=TRUE, show.warnings=TRUE) {
#   if(is.null(beta.start)) beta.start = log(n/m) 
#   if(length(beta.start) != length(n)) stop('beta.start must be same length as n')
#   if(length(m)>1 && length(m) != length(n)) stop('m must be same length as n')
#   b = pmax(beta.start, -50)  # prevent initialization too close to 0
#   P = penMat                   
#   for(i in 1:maxit) {
#     b.old = b
#     mu = exp(b)
#     r = n - m*mu
#     w = m*mu
#     #w = pmax(w, 1e-5)
#     diag(P) = diag(penMat) + w  
#     b = solve(P, r + b*w)
#     #- Check convergence
#     if(sum(abs(b-b.old)) < tol) break
#   }
#   if(show.warnings && i == maxit) warning("did not converge in maxit iterations")
#   b = as.vector(b)
#   
#   out = list(fit = exp(b), beta=b, w=as.vector(w), niters=i)
#   #-- Get effective df and SE
#   if(out.edf) {
#     V = diag(solve(P))
#     edf = sum(V*w)
#     out = c(out, list(edf=edf, se=sqrt(V)))
#   }
#   return(out)
# }


