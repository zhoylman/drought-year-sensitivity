#######################################################################

# Script accomponying Hoylman et al., 2022, Nat. Comm.  
# Drought assessment has been outpaced by climate change: 
# empirical arguments for a paradigm shift
# https://doi.org/10.1038/s41467-022-30316-5

# Author: Dr. Zachary Hoylman
# Contact: zachary.hoylman@umontana.edu

#######################################################################

# This script is accomplishes two major objectives. 
# 1. compute the moving window gamma distrobutional parameters for 
# sites with 100 years of record for use in the monte carlo simulation.
# 2. plot the moving window gamma distrobution PDFs, rate and shape parameters
# along with the moving window mean and coeffitient of variation values for 
# afformentioned sites. 

# start with loading nessisary libs

library(tictoc)
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

# define base parameters 
# ID to define time scale, months of interest and minimum
# number of records, coorisponding to "complete data"
# time scale here is the only parameter that needs to be changed
# to run different timescales. 1 = 30 day, 2 = 60 day, 3 = 30 day. 
time_scale_id = 1
time_scale = list(30,60,90)

months_of_interest = list(c(5,6,7,8),
                          c(4,5,6,7,8),
                          c(3,4,5,6,7,8))
n_minimum = list(123,153,184)

# define nice special
`%notin%` = Negate(`%in%`)

# import states to filter and for plotting
states = st_read('~/drought-year-sensitivity/data/shp/conus_states.shp') %>%
  st_geometry()

#source spi fun
source('~/drought-year-sensitivity/R/funs/gamma_fit_spi.R')

# wrapper function for spi_fun that processes precip data and 
# computes spi for different time periods. This function will also output the
# assosiated gamma parameters for plotting.
# NOTE! moving_window boolean must be set to true to clip data to a 30 year window!
# this is done by default. 
gamma_params = function(data, time_scale, index_of_interest, moving_window = T){ 
  tryCatch({ 
    #define some date based variables
    data$day = yday(data$time)
    data$mday = day(data$time)
    data$year = year(data$time)
    data$month = month(data$time)
    
    #define some meta data on the index of interest 
    max_year = data$year[index_of_interest]
    mday_of_interest = data$mday[index_of_interest]
    month_of_interest = data$month[index_of_interest]
    
    #filter data to not be newer then the index of interest
    #this is important to emulate practial constraints on 
    #drought monitoring. It is impossible to use future
    #data to compute the probability distrobution. Although,
    #it should be noted that it is routine practice to 
    #back calculate drought metrics using reference periods that
    #include data that is in the future with respect to the time
    #period of interest. for example, SPI for 1/1/2015 could be computed
    #using the 1991 - 2020 reference period. 
    precip_data = data %>%
      filter(year <= max_year)
    
    #calcualte index vectors for defining breaks used in the summation procedure.
    #the important function here is to compute the indicies for the start 
    #period and end period (summation period). 
    first_date_breaks = which(precip_data$mday == mday_of_interest & precip_data$month == month_of_interest)
    second_date_breaks = first_date_breaks-(time_scale-1)
    
    #if there are negative indexes remove last year (incomplete data range)
    pos_index = which(second_date_breaks > 0)
    first_date_breaks = first_date_breaks[c(pos_index)]
    second_date_breaks = second_date_breaks[c(pos_index)]
    
    #create slice vectors and group by vectors, this will be used for the 
    #summation procedure below. importantly, this is used instead of year
    # identifiers (or other date based identifiers) because summation periods 
    #can span multiple years. 
    for(j in 1:length(first_date_breaks)){
      if(j == 1){
        slice_vec = seq(second_date_breaks[j],first_date_breaks[j], by = 1)
        group_by_vec = rep(j,(first_date_breaks[j] - second_date_breaks[j]+1))
      }
      else{
        slice_vec = append(slice_vec, seq(second_date_breaks[j],first_date_breaks[j], by = 1))
        group_by_vec = append(group_by_vec, rep(j,(first_date_breaks[j] - second_date_breaks[j]+1)))
      }
    }
    
    #preform the summation procedure. this has a few steps.
    data_time_filter = precip_data %>%
      #slice data for appropriate periods
      slice(slice_vec) %>%
      #add the group_by_vec(tor) to the tibble
      tibble::add_column(group_by_vec = group_by_vec)%>%
      #group by the group_by_vec
      group_by(group_by_vec)%>%
      #summation of data (and convert to mm)
      dplyr::summarise(sum = sum(data/10, na.rm = T),
                       #add the year identifier here,
                       #this can be dont for this special
                       #case because summation does not 
                       #span multiple years. the group by/
                       #slice by methodology is more abstractable
                       year = median(year))
    
    if(moving_window == T){
      data_time_filter = data_time_filter %>%
        #this is computed for the most current 30 year time period
        filter(year >= max_year - 29)
    }
    
    #compute date time for day/year of interest
    date_time = precip_data$time[first_date_breaks[length(first_date_breaks)]] %>% as.Date()
    
    #compute gamma distrobution params
    params = gamma_fit_spi(data_time_filter$sum, 'params')
    
    if(anyNA(params) == T){
      output.df = data.frame(time = NA,
                             shape = NA,
                             rate = NA,
                             mean_p = NA,
                             cv_p = NA,
                             n = NA)
    }
    
    if(anyNA(params) == F){
      #define output dataframe and conduct the SPI calculation. SPI is computed using the 
      #afor-defined spi_fun
      output.df = data.frame(time = date_time,
                             shape = params$para[1],
                             rate = 1/params$para[2],
                             mean_p = mean(data_time_filter$sum),
                             cv_p = sd(data_time_filter$sum)/mean(data_time_filter$sum),
                             n = length(data_time_filter$sum))
      
    }
    #basic error handling. generally for if a gamma fit cannot be obtained
    #there is additional error handling in the spi_fun above
  }, error=function(cond) {
    #if there is an error output an empty but equivelent dataframe (place holder)
    output.df = data.frame(time = NA,
                           shape = NA,
                           rate = NA,
                           mean_p = NA,
                           cv_p = NA,
                           n = NA)
  })
  
  return(output.df)
} 

valid_stations = readRDS('~/drought-year-sensitivity/data/valid_stations_70year_summer_baseline.RDS')

# if plotting is true, plotting functions below will be run. 
# we only plot a few station as examples - c('USW00024137', 'USC00111265', 'USC00143239', 'USC00381770')
plotting = F

selected_sites = which(valid_stations$id %in% c('USW00024137', 'USC00111265', 'USC00143239', 'USC00381770')) #for plotting
#selected_sites = 1:length(valid_stations$id) # for monte carlo
#rev up a cluster for parallel computing
cl = makeCluster(detectCores()-1)
#register the cluster for doPar
registerDoParallel(cl)

tic()
#Fig 2 is USC00381770 or index # 1770

foreach(s = selected_sites, .packages = c('rnoaa', 'tidyverse', 'lubridate', 'magrittr',
                               'lmomco', 'sf')) %dopar% {
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
    filter(n == n_minimum[[time_scale_id]]) # defines number of obs per time scale for complete data by year
  
  #filter it for variable, time, completeness, etc
  data_filtered = data_raw %$%
    #its a list, select precip variable
    prcp %>%
    #add some time meta data for sorting/slicing/filtering
    mutate(year = year(date),
           month = month(date),
           mday = mday(date)) %>%
    #select vars of interest
    dplyr::select(date, prcp, year, month, mday)%>%
    rename(time = date, data = prcp)%>%
    #filter for the time period of interest (months) and completeness
    filter(month %in% months_of_interest[[time_scale_id]],
           year %in% complete_years$year)
  
  #define the indicies of interest, August 1st for plotting and Monte Carlo
  indicies_of_interest = which(data_filtered$month %in% 8 & data_filtered$mday %in% 1 & data_filtered$year <= 2020)
  
  # map the spi calculation through the indicies of interest with moving window defined as T
  temp = indicies_of_interest %>%
    purrr::map(function(i){
      export = gamma_params(data_filtered, time_scale[[time_scale_id]], i, moving_window = T)
      return(export)
    })
  
  #compute the time integrated distrobution
  time_integrated = gamma_params(data_filtered, time_scale[[time_scale_id]], 
                                 indicies_of_interest[length(indicies_of_interest)], moving_window = F)
  
  #merge (rbind) the results and order them by time # bind_rows - dplyr
  params_merged = temp %>%
    #bind rows
    bind_rows() %>%
    #drop nas
    drop_na() %>%
    #order in time
    .[order(.$time),] %>% 
    #convert to tibble
    as_tibble() %>%
    #filter for a minimum of 25 years of data for a given moving window
    filter(n >= 25) %>%
    #clean up data and colnames
    mutate(year = year(time)) %>%
    rename(Shape = shape, Rate = rate, `Mean Precipitation (mm)` = mean_p,
           `CV Precipitation` = cv_p)
  
  # save out data if plotting is F
  if(plotting == F){
    #if its greater than 100 years of data, save it out for monte carlo simulations
    if(length(params_merged$time) >= 100){
      saveRDS(params_merged, paste0('~/drought-year-sensitivity/data/params/param_shift_'
                                    , valid_stations$id[s], '_',time_scale[[time_scale_id]], '_days.RDS'))
    }
    #save out clemson Univ. station for figure and monte carlo
    if(s == 1770 & time_scale[[time_scale_id]] == 30){
      saveRDS(params_merged, paste0('~/drought-year-sensitivity/data/params/param_shift_'
                                    , valid_stations$id[s], '_',time_scale[[time_scale_id]], '_days.RDS'))
    }
  }
  
  #if user would like plotting functions, set to true above
  if(plotting == T){
    #only plot if there is more than 60 moving window's of time
    if(length(params_merged$time) >= 60){
      # compute a sypthetic precipitation distrobution to generate the PDFs of precipitaiton over
      # the synthetic distrobution here is defined as a sequence of length 1,000 
      # ranging from 1 to 3x the max obsereved mean value
      synthetic_precip = seq(1, params_merged$`Mean Precipitation (mm)` %>% max *3, length.out = 1000)
      
      # define "out" data frame to store data
      out = data.frame(matrix(ncol = length(params_merged$year), nrow = length(synthetic_precip)))
      
      #compute PDFs using gamma params 
      for(i in 1:length(params_merged$year)){
        shape = params_merged$Shape[i]
        rate = 1/params_merged$Rate[i]
        temp_params = c(shape, rate) %>% as.numeric() %>% vec2par(., type="gam")
        out[,i] = pdfgam(synthetic_precip, temp_params['para'])
      }
      #define year at end of 30 year moving window for plotting
      colnames(out) = params_merged$year
      out$precip = synthetic_precip
      
      # gather the dataset
      out_tibble = out %>%
        gather(key= 'key', value = 'value',-precip)
      
      # define plotting fucntion for the PDFs 
      plot_dist = ggplot()+
        geom_line(data = out_tibble, aes(x = precip, y = value, color = key %>% as.numeric), alpha = 0.5)+
        scale_colour_gradientn(colours = viridis::turbo(n = 100))+
        theme_bw(base_size = 20)+
        ggtitle(' ')+
        geom_line(data = NULL, aes(x = synthetic_precip, y = pdfgam(synthetic_precip, 
                                                                    vec2par(c(time_integrated$shape, 
                                                                              1/time_integrated$rate), type = 'gam'))), color = 'black', size = 1.5)+
        geom_line(data = NULL, aes(x = synthetic_precip, y = pdfgam(synthetic_precip, 
                                                                    vec2par(c(time_integrated$shape, 
                                                                              1/time_integrated$rate), type = 'gam'))), linetype = 'dashed',
                  color = 'white', size = 1)+
        labs(x = 'Accumulated Precipitation', y = 'PDF')+
        theme(legend.position = 'bottom',
              legend.title = element_blank(),
              legend.key.width=unit(2,"cm"),
              plot.title = element_text(hjust = 0.5))
      
      # plot it
      plot_dist
      
      #compute long tibble for facet plotting of parameters
      long_tibble = params_merged %>% 
        select(-n) %>%
        pivot_longer(cols = -c(year,time)) 
      
      # rename for axes
      long_tibble$name = factor(long_tibble$name, levels = (c('Rate','Shape',
                                                              "Mean Precipitation (mm)",
                                                              "CV Precipitation")),ordered = TRUE)
      
      # define the time integrated dataset/ data frame for plotting
      time_integrated_long = time_integrated %>%
        select(-n) %>%
        rename(Rate = rate, Shape = shape, `CV Precipitation` = cv_p,
               `Mean Precipitation (mm)` = mean_p) %>%
        pivot_longer(cols = -c(time))
      
      # re name for axes and for consistancy
      time_integrated_long$name = factor(time_integrated_long$name, levels = (c('Rate','Shape',
                                                                                "Mean Precipitation (mm)",
                                                                                "CV Precipitation")),ordered = TRUE)
      
      # plot parameters using the facet grid 
      plot_param = ggplot(long_tibble, aes(x = year, y = value))+
        geom_smooth(method = 'loess')+
        geom_point()+
        geom_hline(data = time_integrated_long, aes(yintercept = value))+
        facet_wrap(~name, scales = 'free')+
        theme_bw(base_size = 20)+
        labs(x = NULL, y = 'Parameter Value')+
        theme(plot.title = element_text(hjust = 0.5),
              strip.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA))
      
      # arrange the two plots into a single grob
      plot_grid = ggpubr::ggarrange(plot_param, plot_dist)
      
      # add a title and define some other graphical parameters
      final = ggpubr::annotate_figure(plot_grid,
                                      top = ggpubr::text_grob(paste0(time_scale[[time_scale_id]], ' Day Timescale for August 1\nGHCN Site #: ',valid_stations$id[s], ' (', valid_stations$name[s], ', '
                                                                     , valid_stations$state[s], ')'),
                                                              color = "black",
                                                              face = "bold", size = 22))
      
      # save out the final figure
      ggsave(final, file = paste0('~/drought-year-sensitivity/figs/distrobution_shift/site_',
                                  valid_stations$id[s],'_',time_scale[[time_scale_id]],'_day','.png'),
             width = 15, height = 8, units = 'in', dpi = 300)
    } 
  }
}

toc()

stopCluster(cl)

