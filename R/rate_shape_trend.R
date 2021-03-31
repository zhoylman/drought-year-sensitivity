library(rnoaa)
library(tidyverse)
library(lubridate)
library(magrittr)
library(lmomco)
library(doParallel)
library(ggplot2)

time_scale = 90
months_of_interest = c(6,7,8)

spi_params = function(precip_data, time_scale){
  #define some date based variables
  precip_data$day = yday(precip_data$time)
  precip_data$year = year(precip_data$time)
  precip_data$month = month(precip_data$time)
  
  output.list = list()
  
  #Start SPI calculation
  for(t in 1:length(time_scale)){
    for(i in rev((length(precip_data$time)-364):length(precip_data$time))[1]){
      #calcualte index vectors of interest based on time
      first_date_breaks = which(precip_data$day == precip_data$day[i])
      second_date_breaks = first_date_breaks-(time_scale[t]-1)
      
      #if there are negative indexes remove last year (incomplete data range)
      #change this to remove all indexes from both vectors that are negative
      if(!all(second_date_breaks < 0)){
        pos_index = which(second_date_breaks > 0)
        first_date_breaks = first_date_breaks[c(pos_index)]
        second_date_breaks = second_date_breaks[c(pos_index)]
      }
      
      #create slice vectors and group by vectors
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
      
      #slice data for appropriate periods
      data_time_filter = precip_data %>%
        slice(slice_vec) %>%
        tibble::add_column(group_by_vec = group_by_vec)%>%
        group_by(group_by_vec)%>%
        dplyr::summarise(sum = sum(data, na.rm = T))
      
      if(length(data_time_filter$group_by_vec) >= 30){
        #remove zeros because they cause the gamma dist to blow up to Inf
        data_time_filter$sum[data_time_filter$sum == 0] = 0.01
        
        #compute date time for day/year of interest
        date_time = as.POSIXct(paste(precip_data$day[first_date_breaks], precip_data$year[first_date_breaks], sep = "-"), format = "%j-%Y")
        
        #fit gamma distrobution to data
        #fit.gamma = gamma_fit(data_time_filter$sum)
        #calcualte CDF values for the theoretical distrobution
        #fit.cdf = pgamma(data_time_filter$sum, shape = fit.gamma$shape, rate = fit.gamma$rate)
        
        #Unbiased Sample Probability-Weighted Moments (following Beguer ́ıa et al 2014)
        params = tryCatch(
          {
            pwm = pwm.ub(data_time_filter$sum)
            #Probability-Weighted Moments to L-moments
            lmoments_x = pwm2lmom(pwm)
            #fit gamma
            fit.parglo = pargam(lmoments_x)
            #compute probabilistic cdf
            params = fit.parglo$para %>% as.data.frame() %>% t() %>% as.data.frame()
          },
          error = function(e){
            params = rep(NA, 2)
          }
        )
      }
      else {
        params = rep(NA, 4)
      }
    }
    params$mean = mean(data_time_filter$sum)/10
    params$cv = (sd(data_time_filter$sum)/10)/(mean(data_time_filter$sum)/10)
    return(params)
  }
}

stations = ghcnd_stations()
valid_station_ids = readRDS('/home/zhoylman/drought-year-sensitivity/data/valid_stations.RDS')

#initial filter for stations with potentially enough data
valid_stations = stations %>%
  filter(element == 'PRCP',
         first_year <= 1940,
         last_year == 2020,
         latitude > 22 & latitude < 50,
         longitude < -50,
         id %in% valid_station_ids$id)

for(s in 1:length(valid_stations$id)){
  data_raw = ghcnd_search(
    valid_stations$id[s],
    date_min = NULL,
    date_max = NULL,
    var = "PRCP"
  ) 
  
  annual_nobs = data_raw %$%
    prcp %>%
    mutate(year = year(date),
           month = month(date)) %>%
    drop_na() %>%
    filter(month %in% months_of_interest) %>%
    group_by(year) %>%
    summarize(n = length(prcp)) %>%
    filter(n >= 90) 
  
  print(paste0('Years = ', length(annual_nobs$year)))
  
  data_filtered = data_raw %$%
    prcp %>%
    mutate(year = year(date)) %>%
    filter(year %in% annual_nobs$year) %>% # filter for full years
    select(date, prcp)%>%
    rename(time = date, data = prcp)%>%
    mutate(month = month(time)) %>%
    filter(month %in% months_of_interest)
  
  moving_window_index = data.frame(first = annual_nobs$year[1]:annual_nobs$year[length(annual_nobs$year)]) %>%
    mutate(last = first + 40) %>%
    filter(last <= max(annual_nobs$year)) 
  
  params_out = data.frame(matrix(ncol = 4, nrow = length(moving_window_index$first)))
  colnames(params_out) = c('Alpha (Shape)', 'Beta (Rate)', 
                           'Mean Precipitation (mm; June - Aug)', 'CV Precipitation (mm; June - Aug)')
  for(i in 1:length(moving_window_index$first)){
    params_out[i,] = spi_params(data_filtered %>%
                                  mutate(year = year(time)) %>%
                                  filter(year >= moving_window_index$first[i] &
                                           year <= moving_window_index$last[i]), 90)
  }
  
  params_out_tibble = params_out %>%
    mutate(time = moving_window_index$last) %>%
    gather('key', 'value', -time)
  
  plot = ggplot(params_out_tibble, aes(x = time, y = value))+
    geom_smooth(method = 'loess')+
    geom_point()+
    facet_wrap(~key, scales = 'free')+
    theme_bw(base_size = 16)+
    labs(x = NULL, y = 'Parameter Value')+
    ggtitle(paste0('GHCN Site #: ',valid_stations$id[s]))+
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(plot, file = paste0('/home/zhoylman/drought-year-sensitivity/figs/parameters/site_',
                             valid_stations$id[s],'_',time_scale,'_day','.png'),
         width = 11, height = 8, units = 'in', dpi = 300)
  print(s)
}

