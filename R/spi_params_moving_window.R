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
time_scale_id = 1
time_scale = list(30,60,90)

months_of_interest = list(c(5,6,7,8),
                          c(4,5,6,7,8),
                          c(3,4,5,6,7,8))
n_minimum = list(123,153,184)

# define nice special
`%notin%` = Negate(`%in%`)

# import states to filter and for plotting
states = st_read('/home/zhoylman/mesonet-dashboard/data/shp/states.shp') %>%
  filter(STATE_ABBR %notin% c('AK', 'HI', 'VI')) %>%
  st_geometry()

# function to compute SPI
gamma_fit_spi = function(x, export_opts = 'SPI') {
  #load the package needed for these computations
  library(lmomco)
  #first try gamma
  tryCatch(
    {
      x = as.numeric(x)
      #if precip is 0, replace it with 0.01mm Really Dry
      if(any(x == 0, na.rm = T)){
        index = which(x == 0)
        x[index] = 0.01
      }
      #Unbiased Sample Probability-Weighted Moments (following Beguer ́ıa et al 2014)
      pwm = pwm.ub(x)
      #Probability-Weighted Moments to L-moments
      lmoments_x = pwm2lmom(pwm)
      #fit gamma
      fit.gam = pargam(lmoments_x)
      #compute probabilistic cdf 
      fit.cdf = cdfgam(x, fit.gam)
      #compute spi
      spi = qnorm(fit.cdf, mean = 0, sd = 1)
      if(export_opts == 'CDF'){
        return(fit.cdf) 
      }
      if(export_opts == 'params'){
        return(fit.gam) 
      }
      if(export_opts == 'SPI'){
        return(spi) 
      }
    },
    #else return NA
    error=function(cond) {
      return(NA)
    })
}

# wrapper function for spi_fun that processes precip data and 
# computes spi for different time periods. The pricipal purpose 
# of this function is to properly 
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
    if(!all(second_date_breaks < 0)){
      pos_index = which(second_date_breaks > 0)
      first_date_breaks = first_date_breaks[c(pos_index)]
      second_date_breaks = second_date_breaks[c(pos_index)]
    }
    
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
    
    #remove zeros because they cause the gamma dist to blow up to Inf
    data_time_filter$sum[data_time_filter$sum == 0] = 0.01
    
    #compute date time for day/year of interest
    date_time = precip_data$time[first_date_breaks[length(first_date_breaks)]] %>% as.Date()
    
    #compute gamma distrobution params
    params = gamma_fit_spi(data_time_filter$sum, 'params')
    
    if(is.na(params) == T){
      output.df = data.frame(time = NA,
                             shape = NA,
                             rate = NA,
                             mean_p = NA,
                             cv_p = NA,
                             n = NA)
    }
    
    if(is.na(params) == F){
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

valid_stations = readRDS('/home/zhoylman/drought-year-sensitivity/data/valid_stations_70year_summer_baseline.RDS')

plotting = T

selected_sites = which(valid_stations$id %in% c('USW00024137', 'USC00111265', 'USC00143239', 'USC00381770'))
#selected_sites = which(valid_stations$state == 'MD')
#rev up a cluster for parallel computing
cl = makeCluster(detectCores()-1)
#register the cluster for doPar
registerDoParallel(cl)

tic()

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
  
  #define the indicies of interest, June - August
  indicies_of_interest = which(data_filtered$month %in% 8 & data_filtered$mday %in% 1)
  
  # map teh spi calculation through the indicies of interest 
  temp = indicies_of_interest %>%
    purrr::map(function(i){
      export = gamma_params(data_filtered, time_scale[[time_scale_id]], i, moving_window = T)
      return(export)
    })
  
  time_integrated = gamma_params(data_filtered, time_scale[[time_scale_id]], 
                                 indicies_of_interest[length(indicies_of_interest)], moving_window = F)
  
  #merge (rbind) the results and order them by time # bind_rows - dplyr
  params_merged = temp %>%
    bind_rows() %>%
    drop_na() %>%
    .[order(.$time),] %>% 
    as_tibble() %>%
    filter(n >= 25) %>%
    mutate(year = year(time)) %>%
    rename(Shape = shape, Rate = rate, `Mean Precipitation (mm)` = mean_p,
           `CV Precipitation` = cv_p)
  
  if(length(params_merged$time) >= 100){
    saveRDS(params_merged, paste0('/home/zhoylman/drought-year-sensitivity/data/params/param_shift_'
                                  , valid_stations$id[s], '_',time_scale[[time_scale_id]], '_days.RDS'))
  }
  
  if(plotting == T){
    if(length(params_merged$time) >= 60){
      
      synthetic_precip = seq(1, params_merged$`Mean Precipitation (mm)` %>% max *3, length.out = 1000)
      
      out = data.frame(matrix(ncol = length(params_merged$year), nrow = length(synthetic_precip)))
      
      for(i in 1:length(params_merged$year)){
        shape = params_merged$Shape[i]
        rate = 1/params_merged$Rate[i]
        temp_params = c(shape, rate) %>% as.numeric() %>% vec2par(., type="gam")
        out[,i] = pdfgam(synthetic_precip, temp_params['para'])
      }
      colnames(out) = params_merged$year
      out$precip = synthetic_precip
      
      out_tibble = out %>%
        gather(key= 'key', value = 'value',-precip)
      
      col = colorRampPalette(rev(c("#8b0000", "#ff0000", "#ffff00", '#00FF00',  '#5aede1',"#0000ff",'#9932CC', '#4B0082')))
      
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
      
      plot_dist
      
      long_tibble = params_merged %>% 
        select(-n) %>%
        pivot_longer(cols = -c(year,time)) 
      
      long_tibble$name = factor(long_tibble$name, levels = (c('Rate','Shape',
                                                              "Mean Precipitation (mm)",
                                                              "CV Precipitation")),ordered = TRUE)
      
      time_integrated_long = time_integrated %>%
        select(-n) %>%
        rename(Rate = rate, Shape = shape, `CV Precipitation` = cv_p,
               `Mean Precipitation (mm)` = mean_p) %>%
        pivot_longer(cols = -c(time))
      
      time_integrated_long$name = factor(time_integrated_long$name, levels = (c('Rate','Shape',
                                                                                "Mean Precipitation (mm)",
                                                                                "CV Precipitation")),ordered = TRUE)
      
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
      
      plot_grid = ggpubr::ggarrange(plot_param, plot_dist)
      
      final = ggpubr::annotate_figure(plot_grid,
                                      top = ggpubr::text_grob(paste0(time_scale[[time_scale_id]], ' Day Timescale for August 1\nGHCN Site #: ',valid_stations$id[s], ' (', valid_stations$name[s], ', '
                                                                     , valid_stations$state[s], ')'),
                                                              color = "black",
                                                              face = "bold", size = 22))
      
      ggsave(final, file = paste0('/home/zhoylman/drought-year-sensitivity/figs/distrobution_shift/site_',
                                  valid_stations$id[s],'_',time_scale[[time_scale_id]],'_day','.png'),
             width = 15, height = 8, units = 'in', dpi = 300)
    } 
  }
}

toc()

stopCluster(cl)

