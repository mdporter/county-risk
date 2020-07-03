###############################################################################
# Estimate COVID-19 Infection Risk in US by county
# 
# Created by:
#   Brad Howlett (SDS), Taylor Derby (SDS), Brittany Durkin (SDS), 
#     Aaron Oliver (SDS), Michael D. Porter (SDS/ESE)
#   School of Data Science (SDS) & Engineering Systems and Environment (ESE), 
#   University of Virginia
#
# Initial Version: 2020/06/29
#
# Notes:
#   - (covid, population) data from USAFacts; demographic data from USDA
#   - Estimates S, I, R values using smoothed positive case data. 
#     Models I(t) = p.test * c(t), where c(t) is cases count and p.test is 
#     probability that an infected individual gets tested. Current treats p.test
#     as fixed over time. 
#   - The initial S, I, R estimates are used to estimate the parameters of a SIR 
#     model (beta, gamma) using recent data (e.g., last 3 weeks). These 
#     parameters are used to make forecasts. It would probably be better to 
#     estimate everything together allowing for time-varying beta, but this is
#     a good starting point. 
#   - The probability of being infectious at time t is I(t)/N 
#     (number infectious at time t / population). 
#   - Because I(t+1) depends on the number of susceptible at time t, counties 
#     where S(t) is low due to previous outbreak(s) (and most people are 
#     recovered/immune) will also have lower forecasted infections.
#   - These models will provide estimates for a population/county:
#     a. proportion currently infected 
#     b. proportion immune (recovered)
#     c. proportion still susceptible
# 
# Limitations:
#   - Highly dependent on model settings like probability of being tested. 
#   - Model parameters fixed over time and same for all counties. 
#   - This analysis does not include any spatial information. This would be 
#     helpful for better capturing cross-county contagion. 
#   - Forecasts are based on the assumption that the current situation will 
#     persist through the forecasting period. If conditions change (e.g., 
#     more/less restrictions), then forecasts can change substantially.
#   - Far from perfect, but will hopefully give decent estimates. 
#   - Need to include uncertainty!
###############################################################################

dir_save = "forecasts"      # directory for saving forecasts

#-------------------------------------------------------------------------#
#-- Load Required Packages and Functions
#-------------------------------------------------------------------------#
library(tidyverse)
source("R/functions.R")  
source("R/load-data.R")  

#-------------------------------------------------------------------------#
#-- Load Data
#-------------------------------------------------------------------------#

#-- County Data (Static)
county = load_population("covid_county_population_usafacts.csv")

#-- County Level Covid data (by date)
cases = load_usafacts("covid_confirmed_usafacts.csv", "cases_total") %>% 
  select(FIPS, date, cases_total)
deaths= load_usafacts("covid_deaths_usafacts.csv", "deaths_total") %>% 
  select(FIPS, date, deaths_total)

covid = full_join(cases, deaths, by=c("FIPS", "date")) %>% 
  group_by(FIPS) %>% arrange(date) %>% 
  mutate(cases = get_daily_counts(cases_total), 
         deaths = get_daily_counts(deaths_total)) %>% ungroup()  

#-- Add population and get ratios
covid2 = covid %>% 
  left_join(county, by="FIPS") %>% 
  group_by(FIPS) %>% arrange(date) %>% 
  mutate(deaths_perc = deaths_total/population, 
         cases_perc = cases_total/population) %>% 
  ungroup() %>% 
  filter(population > 0) # remove regions with no recorded population

#-- Demographic Data
url_demographic = "https://www.ers.usda.gov/webdocs/DataFiles/48747/Unemployment.csv?v=8130.5"
demographics = read_csv(url_demographic) %>% 
  select(FIPS = FIPStxt, everything(), -Stabr, -area_name) %>% 
  filter(FIPS %in% county$FIPS) %>% 
  mutate_at(vars(Rural_urban_continuum_code_2013,
                 Urban_influence_code_2013, Metro_2013), as.factor)

#-- Date of most recent data
current_date = max(covid$date) %>% as.character()

#-------------------------------------------------------------------------#
#-- Model and Forecast
#-------------------------------------------------------------------------#

#-- Set Parameter Values:
pars = list(
  k=60,          # length of forecast period
  p.test=1/10,   # probability infected person is tested 
  k.beta=21,     # number of days to use for estimating beta
  gamma=1/14,    # gamma parameter for SIR (prob recover)
  edf=8, deg=3   # parameters for case smoothing
)

#-- Estimate and Make Forecasts
covid_fitted = covid2 %>%  
  group_by(FIPS) %>% 
  do(forecast_SIR(., k=pars$k, p.test=pars$p.test, k.beta=pars$k.beta, 
                     gamma=pars$gamma, edf=pars$edf, deg=pars$deg)) 



#-- Save Results
write_rds(covid_fitted, 
          path=file.path(dir_save, paste0(current_date, ".rds")), 
          compress="gz")

#-- Load Past Results
# covid_fitted = read_rds(file.path(dir_save, paste0('2020-06-30', ".rds")))

#-------------------------------------------------------------------------#
#-- Plots and Tables
my_theme <- function(...) theme_bw() + theme(...)
#-------------------------------------------------------------------------#

forecast_date = "2020-08-01"


#-- Table at forecast period
covid_table = covid_fitted %>% ungroup() %>% 
  filter(date == as.Date(forecast_date)) %>% 
  select(-cases) %>% 
  left_join(county %>% select(FIPS, County, State), by='FIPS') %>% 
  mutate(County = str_replace(County, "County", "Co.")) %>% 
  arrange(-p.inf) %>% 
  select(FIPS, County, State, date, p.inf, p.rec, p.sus, est.cases=C, 
         population=N) %>% 
  left_join(demographics %>% select(FIPS, Unemployment_rate_2019, 
                        Median_Household_Income_2018,
                        Rural_urban_continuum_code_2013), by="FIPS")


#-- Epicurves for most risky places
top_FIPS = covid_table %>% arrange(-p.inf) %>% slice(1:25) %>% pull(FIPS)

covid_fitted %>% ungroup() %>% 
  filter(FIPS %in% top_FIPS) %>% 
  left_join(county %>% select(FIPS, County, State), by='FIPS') %>%
  mutate(County = str_replace(County, "County", "Co.")) %>% 
  mutate(place_name = paste(County, State, sep=', ')) %>% 
  ggplot(aes(date, color=date>current_date)) +
  geom_vline(xintercept = as.Date(current_date), lty=3, color="black") + 
  geom_line(aes(y=p.inf)) +
  scale_color_manual(values=c("black", "red"), guide=FALSE) + 
  coord_cartesian(xlim=c(as.Date('2020-05-01'), as.Date(forecast_date))) + 
  facet_wrap(~place_name) + 
  my_theme(axis.title.x=element_blank()) + 
  labs(y="Pr(Infectious)", title="Most Risky Counties") 


#-- Epi-Curve for a state
state = "VA"

covid_fitted %>% ungroup() %>% 
  filter(between(date, as.Date('2020-05-01'), as.Date(forecast_date))) %>% 
  left_join(county %>% select(FIPS, County, State), by='FIPS') %>%
  filter(State == !!state) %>% 
  mutate(County = str_replace(County, "County", "Co.")) %>% 
  mutate(place_name = paste(County, State, sep=', ')) %>% 
  ggplot(aes(date, color=date>current_date)) +
  geom_vline(xintercept = as.Date(current_date), lty=3, color="black") + 
  geom_line(aes(y=p.inf)) +
  scale_color_manual(values=c("black", "red"), guide=FALSE) + 
  facet_wrap(~place_name) + 
  my_theme(axis.title.x=element_blank()) + 
  labs(y="Pr(Infectious)", title=paste('Counties in', state))


#-- Bar type heat plot for a state
state = "VA"

covid_fitted %>% 
  filter(between(date, as.Date('2020-05-01'), as.Date(forecast_date))) %>% 
  left_join(county %>% select(FIPS, County, State), by='FIPS') %>% 
  filter(State == !!state) %>% 
  mutate(County = str_replace(County, "County", "Co.")) %>% 
  mutate(place_name = paste(County, State, sep=', ')) %>% 
  ggplot(aes(x=date, xend=date, y=0, yend=1, color=p.inf)) + 
  geom_segment() + 
  geom_vline(xintercept = as.Date(current_date), lty=3, color="black") + 
  facet_wrap(~place_name) + 
  scale_color_viridis_c(direction=-1, option = "C", name='Pr(Infectious)', 
                        limits=c(0, NA)) + 
  labs(title=paste('Counties in', state)) + 
  my_theme(axis.title = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y= element_blank(), 
        panel.grid = element_blank(),
        legend.position=c(1,0), legend.direction="horizontal",
        legend.justification=c(1,0))


#-- Distribution of Risk Scores
covid_table %>% 
  ggplot(aes(x=100*p.inf)) + geom_histogram(binwidth=.25, boundary=0) + 
  scale_x_continuous(breaks=seq(0, 100, by=1)) + 
  labs(x=paste("Percent Infected on", forecast_date), 
       title="Distribution of Risk Scores") +
  my_theme()

#-- Risk Scores by Demographics

# P.inf vs. Population
covid_table %>% 
  ggplot() + 
  geom_point(aes(population, p.inf)) + 
  scale_x_log10(labels = scales::comma) + 
  labs(y="Pr(Infectious)") + my_theme()

# P.inf vs. Median Household Income   
covid_table %>% 
  ggplot() + 
  geom_point(aes(Median_Household_Income_2018, p.inf)) +
  scale_x_continuous(labels = scales::comma) + 
  labs(y="Pr(Infectious)")  + my_theme() 


#-- Percent change (future-today)/today
covid_fitted %>% 
  mutate(p.change =  (p.inf - p.inf[date == Sys.Date()])) %>% 
  filter(date == as.Date(forecast_date)) %>% 
  arrange(p.change)

#-- Current vs. Forecasted
# TODO


#-------------------------------------------------------------------------#
#-- Risk Map
#-------------------------------------------------------------------------#

# Important note, the county polygons are loaded using "urbnmapr" R package 
#   which can be installed using: 
#   devtools::install_github("UrbanInstitute/urbnmapr")
library(urbnmapr)

#-- get county polygons as sf object
county_sf = urbnmapr::get_urbn_map(map="counties", sf=TRUE) %>% 
  st_transform(crs="+proj=merc")

#-- Make Map
state = "VA"

covid_fitted %>% ungroup() %>% 
  filter(date == as.Date(forecast_date)) %>% 
  full_join(county_sf, by=c('FIPS' = 'county_fips')) %>% 
  filter(state_abbv == !!state) %>% 
  ggplot() + 
  geom_sf(aes(geometry=geometry, fill=p.inf), color=alpha("black", .10), size=.4) + 
  scale_fill_gradient2(low = "green", high = "red", mid = "yellow", 
                       name="Pr(Infectious)") + 
  labs(title=paste('Forecasted Risk for', 
                   format(as.Date(forecast_date), '%b%e, %Y'),
                   'in', state)) + 
  my_theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(), 
        panel.grid=element_blank(), 
        panel.background=element_rect(fill="black"))

