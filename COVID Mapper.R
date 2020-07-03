# A script which allows for county-level plotting, uses functions from the load-data.R file

# Important note, requires package "urbnmapr" which can be installed using:
# devtools::install_github("UrbanInstitute/urbnmapr")
library(tidyverse)
library(ggmap)
library(urbnmapr)
library(scales)

# This bit is not required, loads in NYT data on COVID case count, used in example case
{
case_data <- read.csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/live/us-counties.csv")
names(case_data)[4] <- "county_fips"
names(county_pop)[4] <- "county_fips"
}

# Function which loads in the county population orginially from load-data.r file
{
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
county_pop <-load_population()
names(county_pop)[4] <- "county_fips"
}

#-------------------------------- IMPORTANT NOTE ------------------------------------
# This file joins databases together on a column named "county_fips", both dataframes
# must have that as the title of the FIPS column to ensure proper functionality
#-------------------------------- IMPORTANT NOTE ------------------------------------

# replace the X below with the path to the location of the county level data stored 
# in a csv file and run the below three lines:
file <- (X)
case_data <- read.csv(file)

# Rename whichever column stores your zip codes/FIPS info to "county_fips"
# Replace the X below with the index of your data's FIPS column
names(case_data)[X] <- "county_fips"

# Pads FIPS column with leading 0 if needed
case_data$county_fips <- sprintf("%05d",case_data$county_fips)

# Joins both dataframes on county_fips
pop_covid <- left_join(case_data, county_pop, by = "county_fips")

# Joins plotting data from the urbnmapr package
county_covid_data <- left_join(pop_covid, counties, by = "county_fips") 

# Creates infection rate variable
#county_covid_data$inf_rate <- county_covid_data$cases / county_covid_data$population



#------------------------------------------------------------------------------
########                            PLOTS                             #########
#------------------------------------------------------------------------------

# State plot
# Continuous variable
# Filtered on VIRGINIA as is, can put any state
# Showing DEATHS as is, can put any variable

county_covid_data %>%
  filter(state_name =="Virginia") %>% 
  ggplot(aes(long, lat, group = group, fill = deaths)) +
  geom_polygon(color = NA) +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
  labs(fill = "Deaths") +
  scale_fill_gradient(low = "darkgreen", high = "red") +
  ggtitle("Deaths by county VA") +
  geom_polygon(mapping = aes(long, lat,group = group),
               fill = NA, color = "#000000", size = 0.4) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())




# Percentage variable (Creates a plot and casts a numeric variable as a percentage)

county_covid_data %>%
  filter(state_name =="Virginia") %>% 
  ggplot(aes(long, lat, group = group, fill = inf_rate)) +
  geom_polygon(color = NA) +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
  labs(fill = "Infection rate") +
  scale_fill_gradient2(low = "green", high = "red", mid = "yellow", labels = percent) +
  ggtitle("Infection rate (cases / pop) by county Virginia") +
  geom_polygon(mapping = aes(long, lat,group = group),
               fill = NA, color = "#000000", size = 0.4) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
