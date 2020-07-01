##############################################################################
# Functions to load COVID and population data from USAFacts
##############################################################################
library(tidyverse)  # most functions require readr, dplyr, stringr, tidyr, etc


# load_usafacts()
#-----------------------------------------------------------------------------#
# Load county level data from USAFacts 
# https://usafacts.org/visualizations/coronavirus-covid-19-spread-map/
#-----------------------------------------------------------------------------#
load_usafacts <- function(file, value="cases_total") {
  ## retrieve raw data
  url_base = "https://usafactsstatic.blob.core.windows.net/public/data/covid-19"
  read_csv(file.path(url_base, file), 
           col_types = cols(.default="d",
                            countyFIPS="c", 
                            `County Name`="c", State="c", stateFIPS="c")) %>% 
    ## add proper 5 digit FIPS  
    mutate(FIPS = str_pad(countyFIPS, width=5, pad="0", side = 'left')) %>%     
    ## remove statewide unallocated counts (**note: ignoring these counts**)  
    filter(countyFIPS != '0', countyFIPS != '1') %>% 
    ## keep FIPS    
    select(-countyFIPS, -stateFIPS) %>%     
    ## convert to long format for date
    pivot_longer(matches("\\d{1,2}/\\d{1,2}/\\d{2}"), # date
                 names_to="date", values_to=value) %>% 
    ## make date a proper date object      
    mutate(date = as.Date(date, "%m/%d/%y")) %>% 
    ## arrange and return      
    arrange(FIPS, date) %>% rename(County = `County Name`)
}

# load_population()
#-----------------------------------------------------------------------------#
# Load population data from USAFacts 
# https://usafacts.org/visualizations/coronavirus-covid-19-spread-map/
#-----------------------------------------------------------------------------#
load_population <- function(file="covid_county_population_usafacts.csv") {
  ## retrieve raw data  
  url_base = "https://usafactsstatic.blob.core.windows.net/public/data/covid-19"
  read_csv(file.path(url_base, file),
           col_types = cols(population="d",
                            countyFIPS="c", 
                            `County Name`="c", State="c")) %>%
    ## remove statewide unallocated rows 
    filter(countyFIPS != '0', countyFIPS != '1') %>% 
    ## add proper 5 digit FIPS     
    mutate(FIPS = str_pad(countyFIPS, width=5, pad="0", side = 'left')) %>% 
    ## arrange and return    
    select(-countyFIPS) %>% arrange(FIPS) %>% rename(County = `County Name`)
}

# get_daily_counts()
#-----------------------------------------------------------------------------#
# Convert cumulative counts to daily counts 
# Notes:
#   - todo: better method for handling decreases in cumulative counts
#-----------------------------------------------------------------------------#
get_daily_counts <- function(x) {
  lag.x = dplyr::lag(x)
  y = x - ifelse(is.na(lag.x), 0, lag.x)  # treat NA's as zeros
  y = pmax(y, 0)                          # treat negative counts as zeros
}