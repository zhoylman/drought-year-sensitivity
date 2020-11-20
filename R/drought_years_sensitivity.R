library(rnoaa)
library(tidyverse)
library(lubridate)
library(magrittr)
library(lmomco)
library(doParallel)
library(ggplot2)

#define some functions
#compute timeseries of SPI using the 
#Unbiased Sample Probability-Weighted Moments method (following Beguer ́ıa et al 2014)
spi_point_ghcn_moments = function(precip_data, time_scale){
  #define some date based variables
  precip_data$day = yday(precip_data$time)
  precip_data$year = year(precip_data$time)
  precip_data$month = month(precip_data$time)
  
  output.list = list()
  
  #Start SPI calculation
  for(t in 1:length(time_scale)){
    for(i in rev((length(precip_data$time)-364):length(precip_data$time))){
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
      
      #remove zeros because they cause the gamma dist to blow up to Inf
      data_time_filter$sum[data_time_filter$sum == 0] = 0.01
      
      #compute date time for day/year of interest
      date_time = as.POSIXct(paste(precip_data$day[first_date_breaks], precip_data$year[first_date_breaks], sep = "-"), format = "%j-%Y")
      
      #fit gamma distrobution to data
      #fit.gamma = gamma_fit(data_time_filter$sum)
      #calcualte CDF values for the theoretical distrobution
      #fit.cdf = pgamma(data_time_filter$sum, shape = fit.gamma$shape, rate = fit.gamma$rate)
      
      #Unbiased Sample Probability-Weighted Moments (following Beguer ́ıa et al 2014)
      fit.cdf = tryCatch(
        {
          pwm = pwm.ub(data_time_filter$sum)
          #Probability-Weighted Moments to L-moments
          lmoments_x = pwm2lmom(pwm)
          #fit gamma
          fit.parglo = pargam(lmoments_x)
          #compute probabilistic cdf
          fit.cdf = cdfgam(data_time_filter$sum, fit.parglo)
        },
        error = function(e){
          fit.cdf = rep(NA, length(data_time_filter$sum))
        }
      )
      
      #equaprobaility transformation for cdf quantiles
      if(i == length(precip_data$time)){
        output.df = data.frame(time = date_time, 
                               spi = qnorm(fit.cdf, mean = 0, sd = 1))
      }
      
      else{
        output.df = rbind(output.df, data.frame(time = date_time, 
                                                spi = qnorm(fit.cdf, mean = 0, sd = 1)))
      } 
    }
    output.df = output.df[order(output.df$time),]
    output.list[[t]] = output.df
  }
  #if there is only one timescale to calculate return a data frame
  if(length(time_scale) == 1){
    return(output.df)
  }
  # otherwise return a list
  else{
    return(output.list)
  }
}

`%notin%` = Negate(`%in%`)

error_comp = function(X){
  random_years = sample(unique(year(data_filtered$time)), X, replace = F)
  remaining_years = unique(year(data_filtered$time))[unique(year(data_filtered$time)) %notin% random_years]
  test_year = sample(remaining_years,1)
  
  temp_filtered = data_filtered %>%
    filter(year(time) %in% c(random_years, test_year)) 
  
  temp_spi = spi_point_ghcn_moments(temp_filtered, 30)
  
  temp_test = temp_spi %>% filter(year(time) == test_year)
  temp_real = full_spi %>% filter(year(time) == test_year)
  
  joined = left_join(temp_real, temp_test, by='time')
  
  error = joined$spi.x - joined$spi.y
  
  error[error == 'NaN'] = NA
  error_abs = abs(error)
  mae = mean(error_abs, na.rm = T)
  
  return(mae)
  print(X)
}

stations = ghcnd_stations()

filtered_stations = stations %>%
  filter(element == 'PRCP',
         first_year <= 1900,
         last_year == 2020,
         latitude > 22 & latitude < 50,
         longitude < -50) 

set.seed(8)
random_index = sample(1:length(filtered_stations$id), 200)
random_stations = filtered_stations[random_index,]

monte_carlo_out = list()
monte_carlo_summary = list()

for(s in 53:length(random_stations$id)){
#for(s in 1:2){
  time = Sys.time()
  data_raw = ghcnd_search(
    random_stations$id[s],
    date_min = NULL,
    date_max = NULL,
    var = "PRCP",
  ) 
  
  annual_nobs = data_raw %$%
    prcp %>%
    mutate(year = year(date)) %>%
    drop_na() %>%
    group_by(year) %>%
    summarize(n = length(prcp)) %>%
    filter(n >= 365) 
  
  print(paste0('Years = ', length(annual_nobs$year)))
  
  if(length(annual_nobs$year) > 80){data_filtered = data_raw %$%
    prcp %>%
    mutate(year = year(date)) %>%
    filter(year %in% annual_nobs$year) %>% # filter for full years
    select(date, prcp)%>%
    rename(time = date, data = prcp)
  
    full_spi = spi_point_ghcn_moments(data_filtered, 30)
    
    years = 1:(length(annual_nobs$year)-1)
    
    years_out = data.frame(years = years,
                           out = NA)
    
    ###################################
    cl = makeCluster(7)
    registerDoParallel(cl)
    
    years = 1:(length(annual_nobs$year)-1)
    n.simulations = 20
    set.seed(10)
    monte_carlo = list()
    
    monte_carlo = foreach(i = 1:n.simulations) %dopar% {
      library(tidyverse)
      library(lubridate)
      library(magrittr)
      library(lmomco)
      test = sapply(years, FUN = error_comp)
      test
    }
    
    stopCluster(cl)
    
    monte_carlo_out[[s]] = monte_carlo
    
    summary = do.call('cbind', monte_carlo) %>%
      as.matrix() %>%
      apply(., 1, function(x){quantile(x, c(0.1,0.5, 0.9), na.rm = T)}) %>%
      t() %>%
      as.data.frame()
    
    monte_carlo_summary[[s]] = summary
    
    plot = ggplot()+
      geom_ribbon(data = summary, aes(x = years, ymax = `90%`, ymin = `10%`), fill = "grey70")+
      geom_line(data = summary, aes(x = years, y = `50%`), size = 2)+
      theme_bw(base_size = 16)+
      labs(y = 'Mean Absolute Error', x = 'Years in Climatology')+
      ggtitle(paste0('GHCN Site #: ',random_stations$id[s], '\n', random_stations$latitude[s], ',', random_stations$longitude[s]))+
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(plot, file = paste0('/home/zhoylman/drought-year-sensitivity/figs/site_',random_stations$id[s],'.png'),)
    print(paste0('Site ', s, ' of ', 200, ' complete! Time = ', Sys.time()- time))
  }
  else{
    print(paste0('Less than 80 Years of Data for Site ', s, ' of ', 200))
  }
}

saveRDS(monte_carlo_out, '/home/zhoylman/drought-year-sensitivity/data/monte_carlo_out.RDS')
saveRDS(monte_carlo_summary, '/home/zhoylman/drought-year-sensitivity/data/monte_carlo_summary_out.RDS')

filtered_mc = Filter(length, monte_carlo_out)

final_list = list()

for(i in 1:12){
  final_list[[i]] = do.call('cbind', filtered_mc[[i]]) %>%
    as.data.frame()
  final_list[[i]]$row_id = 1:length(final_list[[i]][,1])
}

do.call('left_join', c(final_list, by = 'row_id'))

final_df = final_list %>% reduce(left_join, by = "row_id") %>%
  filter(row_id < 81) %>%
  select(-row_id)

final_df[final_df == 'NaN'] = NA

summary_mc = final_df %>%
  apply(., 1, function(x){quantile(x, c(0.1,0.5, 0.9), na.rm = T)}) %>%
  t() %>%
  as.data.frame()

plot = ggplot()+
  geom_ribbon(data = summary_mc, aes(x = years[1:80], ymax = `90%`, ymin = `10%`), fill = "grey70")+
  geom_line(data = summary_mc, aes(x = years[1:80], y = `50%`), size = 2)+
  theme_bw(base_size = 16)+
  labs(y = 'Mean Absolute Error', x = 'Years in Climatology')+
  ggtitle('Merged Results for 12 GHCN Sites')+
  theme(plot.title = element_text(hjust = 0.5))

ggsave(plot, file = paste0('/home/zhoylman/drought-year-sensitivity/figs/all_sites.png'),)
