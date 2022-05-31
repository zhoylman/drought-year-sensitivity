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
library(automap)
library(gstat)
options(dplyr.summarise.inform = FALSE)


time_scale_id = 1
time_scale = list(30,60,90)

months_of_interest = list(c(5,6,7,8),
                          c(4,5,6,7,8),
                          c(3,4,5,6,7,8))
n_minimum = list(123,153,184)

contemporary_climatology_length = 30


#read in dataframe of valid stations
valid_stations = readRDS('/home/zhoylman/drought-year-sensitivity/data/valid_stations_70year_summer_baseline.RDS')

#generate the list to contain results
time_comparison = list()

#rev up a cluster for parallel computing
cl = makeCluster(detectCores()-1)
#register the cluster for doPar
registerDoParallel(cl)

#time parallel run
tictoc::tic()
#process all data in parralel (30 cores ~ 2.5 hrs), parallel processing by station
time_comparison = foreach(s = 1:length(valid_stations$id),
                         .packages = c('rnoaa', 'tidyverse', 'lubridate', 'magrittr',
                                       'lmomco', 'sf')) %dopar% {
       tryCatch({
         #pull in raw GHCN data
         data_raw = ghcnd_search(
           valid_stations$id[s],
           date_min = NULL,
           date_max = NULL,
           var = "PRCP"
         )
         
         #compute years with complete data. the relative restrictions defining
         #what a complete year is changes depending on the time scale.
         #for example, 30 day time scales only require complete data between May 1 - Aug. 31
         #because the period of interest is June 1 - Aug. 31. However a 60 day time scale
         #requires data back to April 1 and 90 day back to March 1.
         complete_years = data_raw %$%
           #its a list, select precip variable
           prcp %>%
           #drop nas
           drop_na() %>%
           #add some time meta data for sorting/slicing/filtering
           mutate(year = year(date),
                  month = month(date)) %>%
           #select vars of interest
           dplyr::select(date, prcp, year, month)%>%
           #rename the variables for proper function ingestion
           rename(time = date, data = prcp)%>%
           #filter data for the time period of interest
           filter(month %in% months_of_interest[[time_scale_id]]) %>%
           #group by the year ID
           group_by(year) %>%
           #summarize the number of observations in the data
           summarise(n = length(data)) %>%
           #filter for complete datasets
           filter(n == n_minimum[[time_scale_id]]) 
         
         years = 2020:1991
         
         n_years_30 = list()
         for(i in 1:length(years)){
           max_year = years[i]
           min_year = years[i]-29
           
           n_years_30[[i]] = complete_years %>%
             filter(year >= min_year,
                    year <= max_year) %>%
             summarise(n = length(year))
         }
         export = data.table::rbindlist(n_years_30)
         export
       }, error = function(e){
         export = NA
       })
       
}                        
tictoc::toc()

stopCluster(cl)

summary = lapply(time_comparison, function(x) tryCatch({return(sum(x$n < 28))}, error = function(e){return(NA)})) %>%
  unlist()
hist(summary)

#stations lost entirely
sum(summary == 30, na.rm = T)/length(summary)

#full data ratio lost
(summary[summary<28] %>% sum(., na.rm = T))/sum(summary, na.rm = T)
