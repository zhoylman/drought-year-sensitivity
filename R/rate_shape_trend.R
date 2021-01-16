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

time_scale = 60
months_of_interest = c(5,6,7,8)

spi_params = function(precip_data, time_scale){
  #define some date based variables
  precip_data$day = yday(precip_data$time)
  precip_data$year = year(precip_data$time)
  precip_data$month = month(precip_data$time)
  
  output.list = list()
  
  #Start SPI calculation
  for(t in 1:length(time_scale)){
    #for(i in rev((length(precip_data$time)-364):length(precip_data$time))[2]){
    for(i in which(precip_data$time == as.Date(paste0(max(precip_data$year),'-08-01')))){ #parameters computed for Aug 1
      #calcualte index vectors of interest based on time and timescale
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
        dplyr::summarise(sum = sum(data/10, na.rm = T))
      
      if(length(data_time_filter$group_by_vec) >= 29){
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
    params$mean = mean(data_time_filter$sum)
    params$cv = (sd(data_time_filter$sum))/(mean(data_time_filter$sum))
    return(params)
  }
} 

valid_stations = readRDS('/home/zhoylman/drought-year-sensitivity/data/valid_stations_70year_summer.RDS')

#for(s in 1:length(valid_stations$id)){
for(s in 1:50){
  data_raw = ghcnd_search(
    valid_stations$id[s],
    date_min = NULL,
    date_max = NULL,
    var = "PRCP"
  ) 
  
  data_filtered = data_raw %$%
    #its a list, select precip
    prcp %>%
    mutate(year = year(date),
           month = month(date)) %>%
    dplyr::select(date, prcp, year, month)%>%
    rename(time = date, data = prcp)%>%
    filter(month %in% months_of_interest)
  
  moving_window_index = data.frame(first = unique(data_filtered$year)[1]:max(unique(data_filtered$year))) %>%
    mutate(last = first + 39) %>%
    filter(last <= max(unique(data_filtered$year))) 
  
  params_out = data.frame(matrix(ncol = 4, nrow = length(moving_window_index$first)))
  params_out$year = moving_window_index$last
  colnames(params_out) = c('Alpha (Shape)', 'Beta (Rate)', 
                           'Mean Precipitation (mm; 60 day Aug. 1)', 'CV Precipitation (mm; 60 day Aug. 1)', 'year')
  
  for(i in 1:length(moving_window_index$first)){
    precip_data = data_filtered %>%
      filter(year >= moving_window_index$first[i],
             year <= moving_window_index$last[i])
    
    if(length(unique(precip_data$year)) == 40){
      params_out[i,1:4] = spi_params(precip_data, time_scale)
    }
  }
  
  params_all = spi_params(data_filtered, time_scale)
  colnames(params_all) = c('Alpha (Shape)', 'Beta (Rate)', 
                           'Mean Precipitation (mm; 60 day Aug. 1)', 'CV Precipitation (mm; 60 day Aug. 1)')
  params_all_wide = params_all %>%
    gather('key', 'value')
  params_all_vec = params_all[1,1:2] %>% as.numeric() %>% vec2par(., type="gam")
  
  params_out_tibble = params_out %>%
    drop_na() %>%
    gather('key', 'value', -year)
  
  plot_param = ggplot(params_out_tibble, aes(x = year, y = value))+
    geom_smooth(method = 'loess')+
    geom_point()+
    geom_hline(data = params_all_wide, aes(yintercept = value))+
    facet_wrap(~key, scales = 'free')+
    theme_bw(base_size = 14)+
    labs(x = NULL, y = 'Parameter Value')+
    theme(plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  
  params_out_wide = params_out_tibble %>%
    pivot_wider(names_from = key, values_from = value) %>%
    drop_na()
  
  synthetic_precip = seq(1, params_all$`Mean Precipitation (mm; 60 day Aug. 1)`*3, length.out = 1000)
  
  out = data.frame(matrix(ncol = length(params_out_wide$year), nrow = length(synthetic_precip)))
  
  for(i in 1:length(params_out_wide$year)){
    temp_params = params_out_wide[i,2:3] %>% as.numeric() %>% vec2par(., type="gam")
    out[,i] = pdfgam(synthetic_precip, temp_params['para'])
  }
  colnames(out) = params_out_wide$year
  out$precip = synthetic_precip
  
  out_tibble = out %>%
    gather(key= 'key', value = 'value',-precip)
  
  col = colorRampPalette(c("#8b0000", "#ff0000", "#ffff00", '#00FF00', "#0000ff", '#9932CC', '#4B0082'))
  
  plot_dist = ggplot()+
    geom_line(data = out_tibble, aes(x = precip, y = value, color = key %>% as.numeric), alpha = 0.3)+
    scale_colour_gradientn(colours = col(100))+
    theme_bw(base_size = 16)+
    geom_line(data = NULL, aes(x = synthetic_precip, y = pdfgam(synthetic_precip, params_all_vec['para'])), color = 'black', size = 1.5)+
    labs(x = 'Accumulated Precipitation', y = 'PDF')+
    theme(legend.position = 'bottom',
          legend.title = element_blank(),
          legend.key.width=unit(2,"cm"),
          plot.title = element_text(hjust = 0.5))+
    xlim(0,params_all$`Mean Precipitation (mm; 60 day Aug. 1)`*3)
    
  plot_grid = ggpubr::ggarrange(plot_param, plot_dist)
  
  final = ggpubr::annotate_figure(plot_grid,
                  top = ggpubr::text_grob(paste0('GHCN Site #: ',valid_stations$id[s], ' (', valid_stations$name[s], ', '
                                                 , valid_stations$state[s], ')'),
                                          color = "black",
                                          face = "bold", size = 16))
  
  ggsave(final, file = paste0('/home/zhoylman/drought-year-sensitivity/figs/parameters/site_',
                             valid_stations$id[s],'_',time_scale,'_day','.png'),
         width = 15, height = 8, units = 'in', dpi = 300)
  
  print(s)
}
