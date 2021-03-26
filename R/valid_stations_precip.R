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

#months of interest (summer)
months_of_interest = c(5,6,7,8)

#define not in function
`%notin%` = Negate(`%in%`)

#states for clipping
states = st_read('/home/zhoylman/mesonet-dashboard/data/shp/states.shp') %>%
  filter(STATE_ABBR %notin% c('AK', 'HI', 'VI')) %>%
  st_geometry()

#import all ghcn station meta
stations = ghcnd_stations()

#convert to SF
stations_sf = stations %>%
  st_as_sf(., coords = c(x = 'longitude', y = 'latitude'))
st_crs(stations_sf) = st_crs(4326)

#initial filter for stations with potentially enough data and within the spatial domain of interest
filtered_stations = stations_sf %>%
  filter(element == 'PRCP',
         first_year <= 1950,
         last_year == 2020) %>%
  st_intersection(., states)

#rev up the socket cluster
cl = makeCluster(6)
registerDoParallel(cl)

#in parallel, compute the number of valid years in sequence
nobs_list = foreach(s = 1:length(filtered_stations$id))%dopar%{
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
    filter(n == 123)
  
  #define out df
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

#plot for funzies
plot(final_valid$geometry, xlab = '', ylab = '', main = paste0('70+ Years Summeritme Precipitation (n = ', length(final_valid$id), ')'))
plot(states, add = T)

#save it out for later use
saveRDS(final_valid, file = '/home/zhoylman/drought-year-sensitivity/data/valid_stations_70year_summer.RDS')
