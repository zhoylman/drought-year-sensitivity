#######################################################################

# Script accomponying Hoylman et al., 2022, Nat. Comm.  
# Drought assessment has been outpaced by climate change: 
# empirical arguments for a paradigm shift
# https://doi.org/10.1038/s41467-022-30316-5

# Author: Dr. Zachary Hoylman
# Contact: zachary.hoylman@umontana.edu

#######################################################################

# This script computes the preliminary filtering criteria
# for what will be considered a "Valid Station" for further ananlysis. 
# This initial filtering requires the data to meet 2 primary constraints. 
# Firstly, the GHCN station should start recording no later than 1951 
# which provides a minimum of 70 years of data (potentially). 
# Then we need to compute the number of stations with complete "Seasonal" data.
# This is defined using months of interest for summer and the number of observations
# within these months of interest for each site. 
# Months of interest (summer) Logic is as follows: 
# Begin with the least restrictive consideration,
# 30 day SPI beginning June 1, ending August 31. 
# This computation requires complete data from May - August.
# This is the base line definition of "valid stations" after which
# conditions become more constrained for longer timescales (implemented in following scripts).
# For example a 60 day timescale requires data from April - August,
# and 90 day timescale requires data from March - August.

#######################################################################

#load required libraries
library(rnoaa)
library(tidyverse)
library(lubridate)
library(magrittr)
library(lmomco)
library(doParallel)
library(ggplot2)
library(foreach)
library(doParallel)
library(sf)

#define months of interest
months_of_interest = c(5,6,7,8)

#define not in function
`%notin%` = Negate(`%in%`)

#states for clipping
states = st_read('~/drought-year-sensitivity/data/shp/states.shp/states.shp') %>%
  filter(STATE_ABBR %notin% c('AK', 'HI', 'VI')) %>%
  st_geometry()

#import all ghcn station meta
stations = ghcnd_stations()

#convert to SF
stations_sf = stations %>%
  st_as_sf(., coords = c(x = 'longitude', y = 'latitude'))
#Define CRS (Coordinate Reference System)
st_crs(stations_sf) = st_crs(4326)

#initial filter for stations with potentially enough data (70 years) and within the spatial domain of interest
filtered_stations = stations_sf %>%
  filter(element == 'PRCP',
         first_year <= 1950,
         last_year == 2020) %>%
  st_intersection(., states)

#rev up the socket cluster for parallel processing
cl = makeCluster(6)
registerDoParallel(cl)

#in parallel, compute the number of valid years in sequence
nobs_list = foreach(s = 1:length(filtered_stations$id))%dopar%{
  #load required package for each R instance
  library(rnoaa)
  library(tidyverse)
  library(lubridate)
  
  #import raw data
  data_raw = ghcnd_search(
    filtered_stations$id[s],
    date_min = NULL,
    date_max = NULL,
    var = "PRCP"
  )
  
  #compute summer obs per year
  seasonal_obs = data_raw %$%
    prcp %>%
    mutate(year = year(date),
           month = month(date)) %>%
    filter(month %in% months_of_interest) %>%
    drop_na() %>%
    group_by(year) %>%
    summarize(n = length(prcp)) %>%
    #filter for complete data
    filter(n == 123)# May 1 - Aug 31
  
  #define out df id = station ID, nobs = number of complete years
  out = data.frame(id = filtered_stations$id[s], nobs = length(seasonal_obs$year))
  
  #export in foreach
  out
}

#stop cluster
stopCluster(cl)

#combine results to single df for analysis 
valid_stations_params = data.table::rbindlist(nobs_list) %>%
  filter(nobs >=70 ) 

#compute final sf object
final_valid = filtered_stations %>%
  filter(id %in% valid_stations_params$id)

#plot for visual confirmation
plot(final_valid$geometry, xlab = '', ylab = '', main = paste0('70+ Years Summeritme Precipitation (n = ', length(final_valid$id), ')'))
plot(states, add = T)

#save it out for later use
saveRDS(final_valid, file = '~/drought-year-sensitivity/data/valid_stations_70year_summer_baseline.RDS')
