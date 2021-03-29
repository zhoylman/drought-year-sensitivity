#################################################
########### Drought Metric Bias (SPEI) ###########
##########################  #######################

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
library(gstat)
library(raster)
library(stars)
library(automap)
options(dplyr.summarise.inform = FALSE)

#define base parameters
time_scale = 60
months_of_interest = c(4,5,6,7,8)
contemporary_climatology_length = 30

#define nice special  
`%notin%` = Negate(`%in%`)

#import PET algorithm
source('/home/zhoylman/drought-year-sensitivity/R/hargreaves_samani_eto.R')

#import states to filter and for plotting
states = st_read('/home/zhoylman/mesonet-dashboard/data/shp/states.shp') %>%
  filter(STATE_ABBR %notin% c('AK', 'HI', 'VI')) %>%
  st_geometry()

#define the SPEI function to compute SPEI based on defined params
spei_vals = function(diff_data, time_scale){ 
  #define some date based variables
  diff_data$day = yday(diff_data$time)
  diff_data$mday = day(diff_data$time)
  diff_data$year = year(diff_data$time)
  diff_data$month = month(diff_data$time)
  
  #Start SPEI calculation
  for(t in 1:length(time_scale)){ 
    for(i in rev((length(diff_data$time)-62):length(diff_data$time))){ # compute SPEI for 62 days (July 1 - August 31 is target range)
      #calcualte index vectors of interest based on time and timescale
      first_date_breaks = which(diff_data$mday == diff_data$mday[i] & diff_data$month == diff_data$month[i])
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
      data_time_filter = diff_data %>%
        slice(slice_vec) %>%
        tibble::add_column(group_by_vec = group_by_vec)%>%
        group_by(group_by_vec)%>%
        dplyr::summarise(sum = sum(diff, na.rm = T))
      
      #compute date time for day/year of interest
      date_time = diff_data$time[first_date_breaks] %>% as.Date()
      
      #Unbiased Sample Probability-Weighted Moments (following Beguer ́ıa et al 2014)
      pwm = pwm.ub(data_time_filter$sum)
      #Probability-Weighted Moments to L-moments
      lmoments_x = pwm2lmom(pwm)
      #fit generalized logistic
      fit.parglo = parglo(lmoments_x)
      #compute probabilistic cdf 
      fit.cdf = cdfglo(data_time_filter$sum, fit.parglo)
      
      if(i == length(diff_data$time)){
        output.df = data.frame(time = date_time,
                               spei = qnorm(fit.cdf, mean = 0, sd = 1),
                               n = length(data_time_filter$sum))
      }
      else{
        output.df = rbind(output.df, data.frame(time = date_time,
                                                spei = qnorm(fit.cdf, mean = 0, sd = 1),
                                                n = length(data_time_filter$sum)))
      }
    } 
  }
  output.df = output.df[order(output.df$time),]
  return(output.df)
}

#read in dataframe of valid stations
valid_stations = readRDS('/home/zhoylman/drought-year-sensitivity/data/valid_stations_70year_summer_precip_temp.RDS')

#generate the list to contain results
spei_comparison = list()

#rev up a cluster for parallel computing
cl = makeCluster(detectCores()-1)
#register the cluster for doPar
registerDoParallel(cl)

#time parallel run
tictoc::tic() 
spei_comparison = foreach(s = 1:length(valid_stations$id), .packages = c('rnoaa', 'tidyverse', 'lubridate', 'magrittr',
                                                                        'lmomco', 'sf')) %dopar% {
  tryCatch(
    {
      #import function so that all the cores have the script. 
      source('/home/zhoylman/drought-year-sensitivity/R/hargreaves_samani_eto.R')
      #pull in raw GHCN data
      data_raw = ghcnd_search(
        valid_stations$id[s],
        date_min = NULL,
        date_max = NULL,
        var = c('PRCP','TMAX','TMIN')
      )
      
      joined = left_join(data_raw %$% prcp, data_raw %$% tmax, by = c('date', 'id')) %>%
        left_join(data_raw %$% tmin, by = c('date', 'id')) %>%
        dplyr::select(id, date, prcp, tmax, tmin)
      
      #compute summer obs per year
      complete_years = joined %>%
        mutate(year = year(date),
               month = month(date)) %>%
        filter(month %in% months_of_interest) %>%
        drop_na() %>%
        dplyr::select(-month) %>%
        pivot_longer(cols = -c(id,date,year)) %>%
        group_by(year, name) %>%
        summarize(n = length(value)) %>%
        ungroup() %>%
        group_by(year)%>%
        summarize(n_combined = sum(n)) %>%
        #filter for complete data (3 vars for 153 days (5 months - Daily))
        filter(n_combined == 153*3)
      
      #filter it for variable, time etc
      data_filtered = joined %>%
        mutate(year = year(date),
               month = month(date)) %>%
        dplyr::select(date, prcp, tmax, tmin, year, month)%>%
        filter(month %in% months_of_interest,
               year %in% complete_years$year) %>%
        mutate(pet = hs_eto(Tmax = tmax/10, 
                            Tmin = tmin/10, 
                            Julian_day = yday(date), 
                            lat = st_coordinates(valid_stations[s,])[2]),
               diff = pet - prcp/10,
               year = year(date)) %>%
        dplyr::select(date,year,diff) %>%
        rename(time = date)
      
      #filter for the most current 30 years
      # data_contemporary = data_filtered %>%
      #   filter(year > max(.$year) - contemporary_climatology_length)
      
      complete_year_min_contemporary = complete_years$year[(length(complete_years$year) - (contemporary_climatology_length-1))]
      
      data_contemporary = data_filtered %>%
        filter(year >= complete_year_min_contemporary)
      
      #compute SPEI using all data
      spei_historic = spei_vals(data_filtered, time_scale) %>%
        rename(spei_historic = spei, n_historic = n)
      
      #compute SPEI using only the current 30 years
      spei_contemporary = spei_vals(data_contemporary, time_scale) %>%
        rename(spei_contemporary = spei, n_contemporary = n)
      
      #join datasets for the final comparison dataframe
      spei_merged = spei_historic %>%
        filter(time %in% spei_contemporary$time) %>%
        left_join(., spei_contemporary, by = 'time') %>%
        mutate(diff = spei_historic - spei_contemporary,
               complete_year_min_contemporary = complete_year_min_contemporary)
    },
    error = function(e){
      spei_merged = NA
    })  
  
  spei_merged                                             
}
tictoc::toc()

#stop cluster
stopCluster(cl)

#drought breaks to compute bias based on different classes
drought_breaks = c(-0.5, -0.7, -1.2, -1.5, -1.9, -Inf) %>% rev

#define function to compute bias for different drought classes
#determined by the longest period of record SPEI values and the
#UNL definitions of Dx classes (above)
# this analysis is for Daily Values between July 1 - Aug 31
drought_class_bias = function(x){
  #dummy data frame to ensure all levels are present. 
  dummy_df = x[1:5,]
  dummy_df[1:5,] = NA
  dummy_df$drought = as.factor(c('D0', 'D1', 'D2', 'D3', 'D4'))
  dummy_df$month = NA
  #summarise
  temp = x %>%
    mutate(month = month(time))%>%
    filter(month %in% c(7,8))%>%
    mutate(drought = .bincode(spei_historic, drought_breaks)) %>%
    tidyr::drop_na() %>%
    mutate(drought = drought %>% as.factor(), 
           drought = plyr::revalue(drought, c(`1` = 'D4',
                                              `2` = 'D3',
                                              `3` = 'D2',
                                              `4` = 'D1',
                                              `5` = 'D0'))) %>%
    as_tibble()%>%
    rbind(., dummy_df) %>%
    group_by(drought, .drop = FALSE) %>%
    summarise(bias = median(diff, na.rm = T)) 
  
  temp = with(temp, temp[order(drought %>% as.character()),])$bias
  export = data.frame(temp) %>% t %>% as.data.frame()
  colnames(export) = c('D0', 'D1', 'D2', 'D3', 'D4')
  rownames(export) = NULL
  return(export)
}

drought_bias_all = function(x){
  #summarise
  temp = x %>%
    mutate(month = month(time))%>%
    filter(month %in% c(7,8),
           spei_historic <= -0.5)%>%
    tidyr::drop_na() %>%
    as_tibble()%>%
    summarise(bias = median(diff, na.rm = T)) 
  
  return(temp$bias)
}

# compute average bias for all data together (wet and dry, all D classes together)
bias = lapply(spei_comparison, drought_bias_all) %>%
  unlist()

#drought classes broken out
drought_class = lapply(spei_comparison, drought_class_bias) %>%
  data.table::rbindlist(.) %>%
  as_tibble()

#compute the min number of years for the full record
min_clim = lapply(spei_comparison, FUN = function(x){median(x$n_historic)}) %>%
  unlist()

#compute the min number of years for the 30 year record
contemp_clim = lapply(spei_comparison, FUN = function(x){median(x$n_contemporary)}) %>%
  unlist()

contemp_min_year = lapply(spei_comparison, FUN = function(x){median(x$complete_year_min_contemporary)}) %>%
  unlist()

#merge into single dataframe that summarizes all results  
valid_stations_joined = valid_stations %>%
  mutate(`Average Bias` = bias,
         min_clim = min_clim,
         contemp_min_year = contemp_min_year,
         contemp_clim = contemp_clim,
         D0 = drought_class$D0,
         D1 = drought_class$D1,
         D2 = drought_class$D2,
         D3 = drought_class$D3,
         D4 = drought_class$D4)

#filter for climatology with greater than 100 years for the long "historic" analysis 
#and for locations with 30 years in the "contempary record"
# valid_stations_filtered = valid_stations_joined %>%
#   filter(min_clim >= 70,
#          contemp_clim == 30)

valid_stations_filtered = valid_stations_joined %>%
  filter(min_clim >= 70,
         contemp_clim == 30,
         contemp_min_year > 1985)

#################################################
################ Plot the Results ###############
#################################################

#define plotting
col = colorRampPalette((c('darkred', 'red', 'white', 'blue', 'darkblue')))

classes = c('Average Bias','D0', 'D1', 'D2', 'D3', 'D4')

for(c in 1:length(classes)){
  
  pts_plot = ggplot(valid_stations_filtered)+
    geom_sf(data = states)+
    geom_sf(aes(color = get(classes[c])))+
    scale_color_gradientn(colours = col(100), breaks = c(-0.5, 0, 0.5), limits = c(-0.5, 0.5),
                          labels = c('-0.5 (Dry Bias)', '0 (No Bias)', '0.5 (Wet Bias)'), name = "",
                          oob = scales::squish)+
    theme_bw()+
    ggtitle(paste0('Average Difference in Daily Summer SPEI Values (', classes[c], ')\n', time_scale, ' Day SPEI (July 1 - August 31)'))+
    theme(legend.position = 'bottom',
          legend.key.width=unit(2,"cm"),
          plot.title = element_text(hjust = 0.5))
  
  pts_plot
  
  # Krige
  
  template = raster::raster(resolution=c(1/3,1/3),
                            crs = sp::CRS("+init=epsg:4326")) %>%
    raster::crop(., as(states, 'Spatial')) %>%
    raster::rasterToPoints() %>%
    as.data.frame() %>%
    st_as_sf(., coords = c('x', 'y')) 
  
  st_crs(template) = st_crs(4326)
  
  temp_stations = valid_stations_filtered %>%
    drop_na(classes[c]) %>%
    st_as_sf
  
  vgm = autofitVariogram(temp_stations[classes[c]] %>% data.frame() %>% .[,1] ~ 1, as(temp_stations, 'Spatial'))
  krig = krige(temp_stations[classes[c]] %>% data.frame() %>% .[,1] ~ 1, as(temp_stations, 'Spatial'), template, model=vgm$var_model) %>%
    st_intersection(states)
  krig_pts = st_coordinates(krig) %>%
    as_tibble() %>%
    mutate(val = krig$var1.pred)
  
  krig_plot = ggplot(krig)+
    geom_tile(data = krig_pts, aes(x = X, y = Y, fill = val))+
    geom_sf(data = states, fill = 'transparent', color = 'black')+
    labs(x = "", y = "")+
    scale_fill_gradientn(colours = col(100), breaks = c(-0.5, 0, 0.5), limits = c(-0.5, 0.5),
                         labels = c('-0.5 (Dry Bias)', '0 (No Bias)', '0.5 (Wet Bias)'), name = "",
                         oob = scales::squish)+
    theme_bw()+
    theme(legend.position = 'bottom',
          legend.key.width=unit(2,"cm"))
  
  final = cowplot::plot_grid(pts_plot, krig_plot, ncol = 1)

  # vgm1 = variogram(temp_stations[classes[c]] %>% data.frame() %>% .[,1] ~ 1, data = temp_stations)
  # fit1 = fit.variogram(vgm1, model = gstat::vgm("Gau")) # fit model
  # krig = krige(temp_stations[classes[c]] %>% data.frame() %>% .[,1]~1, 
  #                      temp_stations, template, model=fit1) %>%
  #   st_intersection(states)
  # 
  # krig_pts = st_coordinates(krig) %>%
  #   as_tibble() %>%
  #   mutate(val = krig$var1.pred)
  # 
  # krig_plot = ggplot(krig)+
  #   geom_tile(data = krig_pts, aes(x = X, y = Y, fill = val))+
  #   geom_sf(data = states, fill = 'transparent', color = 'black')+
  #   labs(x = "", y = "")+
  #   scale_fill_gradientn(colours = col(100), breaks = c(-0.5, 0, 0.5), limits = c(-0.5, 0.5),
  #                        labels = c('-0.5 (Dry Bias)', '0 (No Bias)', '0.5 (Wet Bias)'), name = "",
  #                        oob = scales::squish)+
  #   theme_bw()+
  #   theme(legend.position = 'bottom',
  #         legend.key.width=unit(2,"cm"))
  # 
  # final = cowplot::plot_grid(pts_plot, krig_plot, ncol = 1)
  # 
  ggsave(final, file = paste0('/home/zhoylman/drought-year-sensitivity/figs/spei_bias_maps_',classes[c],'_',time_scale,'day_timescale_July1-Aug31.png'), width = 7, height = 10, units = 'in')
  #plot
}

