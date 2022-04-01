#################################################
########### Drought Metric Bias (SPI) ###########
#################################################

# moving window based analysis.
# this analysis computes SPI for "valid" GHCN sites
# as defined in ~/drought-year-sensitivity/R/valid_stations_precip.R
# for the longest period of record possible and for  
# the previous 30 year moving window for years with complete data.

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
library(spdplyr)
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
states = st_read('/home/zhoylman/drought-year-sensitivity/data/shp/conus_states.shp') %>%
  st_geometry()

#read in dataframe of valid stations
valid_stations = readRDS('/home/zhoylman/drought-year-sensitivity/data/valid_stations_70year_summer_baseline.RDS')

# function to compute SPI
source('~/drought-year-sensitivity/R/funs/gamma_fit_spi.R')

# wrapper function for gamma_fit_spi that processes precip data and 
# computes spi for different time periods. The pricipal purpose 
# of this function is to properly 
daily_spi = function(data, time_scale, index_of_interest){ 
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
    #this is important to emulate practical constraints on 
    #drought monitoring. It is impossible to use future
    #data to compute the probability distribution in practice. Although,
    #it should be noted that it is common to 
    #back calculate drought metrics using reference periods that
    #include data that is in the future with respect to the time
    #period of interest. for example, SPI for 1/1/2015 could be computed
    #using the 1991 - 2020 reference period. 
    precip_data = data %>%
      filter(year <= max_year)
    
    #calculate index vectors for defining breaks used in the summation procedure.
    #the important function here is to compute the indices for the start 
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
    #can span multiple years (not in this case but indeed in other cases). 
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
                       #this can be done for this special
                       #case because summation does not 
                       #span multiple years. the group by/
                       #slice by methodology is more abstractable
                       year = median(year))
    
    #compute date time for day/year of interest
    date_time = precip_data$time[first_date_breaks[length(first_date_breaks)]] %>% as.Date()
    
    #30 year moving window filter to compute "contempary SPI values"
    contempary_data = data_time_filter %>%
      #this is computed for the most current 30 year time period
      filter(year >= max_year - 29)
    
    #store params
    params_contempary = gamma_fit_spi(contempary_data$sum, 'params')
    params_historical = gamma_fit_spi(data_time_filter$sum, 'params')
  
    #generate storage data frame
    if(is.na(params_contempary) == T | is.na(params_historical) == T){
      output.df = data.frame(time = NA,
                             spi_historical = NA,
                             shape_historical = NA,
                             rate_historical = NA,
                             n_historical = NA,
                             spi_contemporary = NA,
                             shape_contemporary = NA,
                             rate_contemporary = NA,
                             n_contemporary = NA) %>%
        `rownames<-`(1)
    }
    
    if(is.na(params_contempary) == F & is.na(params_historical) == F){
      #define output dataframe and conduct the SPI calculation. SPI is computed using the 
      #afor-defined gamma_fit_spi
      output.df = data.frame(time = date_time,
                             #compute spi using l-moments and the gamma distribution
                             #for the historical time period (longest period of record)
                             spi_historical = gamma_fit_spi(data_time_filter$sum, 'SPI'),
                             #store parameters
                             shape_historical = params_historical$para[1],
                             rate_historical = 1/params_historical$para[2],
                             #report the number of years in this SPI calculation
                             n_historical = length(data_time_filter$sum),
                             #compute SPI using the most current 30 years of data
                             spi_contemporary = gamma_fit_spi(contempary_data$sum, 'SPI'),
                             #store parameters
                             shape_contemporary = params_contempary$para[1],
                             rate_contemporary = 1/params_contempary$para[2],
                             #report number of years in the SPI calculation
                             n_contemporary = length(contempary_data$sum))%>%
        `rownames<-`(1)
    }
    
    #basic error handling. generally for if a gamma fit cannot be obtained
    #there is additional error handling in the gamma_fit_spi above
  }, error=function(cond) {
    #if there is an error output an empty but equivelent dataframe (place holder)
    output.df = data.frame(time = NA,
                           spi_historical = NA,
                           shape_historical = NA,
                           rate_historical = NA,
                           n_historical = NA,
                           spi_contemporary = NA,
                           shape_contemporary = NA,
                           rate_contemporary = NA,
                           n_contemporary = NA)%>%
      `rownames<-`(1)
  })
  
  return(output.df)
} 

#generate the list to contain results
spi_comparison = list()

#rev up a cluster for parallel computing
cl = makeCluster(detectCores()-1)
#register the cluster for doPar
registerDoParallel(cl)

#time parallel run
tictoc::tic()
#process all data in parralel (31 cores ~ 2.5 hrs), parallel processing by station
spi_comparison = foreach(s = 1:length(valid_stations$id),
                         .packages = c('rnoaa', 'tidyverse', 'lubridate', 'magrittr',
                                                                        'lmomco', 'sf')) %dopar% {
    #basic error handling
    tryCatch(
      {
        #pull in raw GHCN data
        data_raw = ghcnd_search(
          valid_stations$id[s],
          date_min = NULL,
          date_max = NULL,
          var = "PRCP"
        ) 
        
        #compute years with complete data. the relative restrictions defining
        #what a complete year changes depending on the time scale. 
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
                 month = month(date)) %>%
          #select vars of interest
          dplyr::select(date, prcp, year, month)%>%
          rename(time = date, data = prcp)%>%
          #filter for the time period of interest (months) and completeness
          filter(month %in% months_of_interest[[time_scale_id]],
                 year %in% complete_years$year)
        
        #define the indicies of interest, June - August
        indicies_of_interest = which(data_filtered$month %in% c(6,7,8) & data_filtered$year > 1990 & data_filtered$year <= 2020)

        # map the spi calculation through the indicies of interest 
        temp = indicies_of_interest %>%
          purrr::map(function(i){
            export = daily_spi(data_filtered, time_scale[[time_scale_id]], i)
            return(export)
          })
        #merge (rbind) the results and order them by time # bind_rows - dplyr
        spi_merged = temp %>%
          bind_rows() %>%
          .[order(.$time),] %>% 
          as_tibble() %>%
          filter_all(any_vars(!is.na(.)))

      },
      #basic error handling
      error = function(e){
        spi_merged = data.frame(time = NA,
                                spi_historical = NA,
                                n_historical = NA,
                                spi_contemporary = NA,
                                n_contemporary = NA)
      })  
    #return the data
    spi_merged                                          
  }
tictoc::toc()

#stop cluster
stopCluster(cl)

#save out big list
saveRDS(spi_comparison, paste0('/home/zhoylman/temp', '/spi_comparision_moving_window_with_params_30year_', time_scale[[time_scale_id]], '_days.RDS'))

#read in big list if already processed
spi_comparison = readRDS(paste0('/home/zhoylman/temp', '/spi_comparision_moving_window_with_params_30year_', time_scale[[time_scale_id]], '_days.RDS'))

#drought breaks to compute bias based on different classes
drought_breaks = c(-0.5, -0.8, -1.3, -1.6, -2, -Inf) %>% rev

#define function to compute bias for different drought classes
#determined by the longest period of record SPI values and the
#UNL definitions of Dx classes (above)
# this analysis is for Daily Values between June 1 - Aug 31
drought_class_bias = function(x){
  #dummy data frame to ensure all levels are present.
  #we will bind this to the final summary data to ensure continuity
  dummy_df = x[1:5,]
  dummy_df[1:5,] = NA
  dummy_df = dummy_df%>%
    mutate(drought = as.factor(c('D0', 'D1', 'D2', 'D3', 'D4')),
           month = NA,
           year = NA,
           diff = NA)
  #summarise
  temp = x %>%
    #define some time meta data
    mutate(month = month(time),
           year = year(time),
           #compute the difference in historical (longest clim)
           #and the contempary data (last 30 years of data)
           diff = spi_historical - spi_contemporary)%>%
    #filter for time period of intest
    filter(month %in% c(6,7,8),
           #filter for a 25 year minimum climatology in the contemporary data
           n_contemporary >= 25,
           #filter for a 70 year minimum climatology in the historical data
           n_historical >= 70)%>%
    #compute the grouping bins
    mutate(drought = .bincode(spi_historical, drought_breaks)) %>%
    #drop data that doesn't fit the classes
    tidyr::drop_na() %>%
    #revalue the data to be human interperatble
    mutate(drought = drought %>% as.factor(),
           drought = plyr::revalue(drought, c(`1` = 'D4',
                                              `2` = 'D3',
                                              `3` = 'D2',
                                              `4` = 'D1',
                                              `5` = 'D0'))) %>%
    as_tibble()%>%
    #rbind our dummy dataframe to ensure all levels are present
    rbind(., dummy_df) %>%
    #refine the group_by structure by drought classes, don't drop missing classes
    group_by(drought, .drop = FALSE) %>%
    #compute the median difference
    summarise(bias = median(diff, na.rm = T))

  #reorganize the final data frame and rename data
  temp = with(temp, temp[order(drought %>% as.character()),])$bias[1:5]
  export = data.frame(temp) %>% t %>% as.data.frame()
  colnames(export) = c('D0', 'D1', 'D2', 'D3', 'D4')
  rownames(export) = NULL
  return(export)
}

#drought breaks to compute bias based on different classes
generalized_dryness = c(Inf, 2, 1, -1, -2, -Inf) %>% rev

#define function to compute bias for different drought classes
#determined by the longest period of record SPI values and the
#UNL definitions of Dx classes (above)
# this analysis is for Daily Values between June 1 - Aug 31
dryness_class_bias = function(x){
  #dummy data frame to ensure all levels are present.
  #we will bind this to the final summary data to ensure continuity
  dummy_df = x[1:5,]
  dummy_df[1:5,] = NA
  dummy_df = dummy_df%>%
    mutate(drought = as.factor(c('SPI > 2', '1 > SPI > 2', '1 > SPI > -1', '-1 > SPI > -2', '-2 > SPI')),
           month = NA,
           year = NA,
           diff = NA)
  #summarise
  temp = x %>%
    #define some time meta data
    mutate(month = month(time),
           year = year(time),
           #compute the difference in historical (longest clim)
           #and the contemporary data (last 30 years of data)
           diff = spi_historical - spi_contemporary)%>%
    #filter for time period of intest
    filter(month %in% c(6,7,8),
           #filter for a 25 year minimum climatology in the contemporary data
           n_contemporary >= 25,
           #filter for a 70 year minimum climatology in the historical data
           n_historical >= 70)%>%
    #compute the grouping bins
    mutate(drought = .bincode(spi_historical, generalized_dryness)) %>%
    #drop data that doesnt fit the classes
    tidyr::drop_na() %>%
    #revalue the data to be human interperatble
    mutate(drought = drought %>% as.factor(),
           drought = plyr::revalue(drought, c(`1` = '-2 > SPI',
                                              `2` = '-1 > SPI > -2',
                                              `3` = '1 > SPI > -1',
                                              `4` = '1 > SPI > 2',
                                              `5` = 'SPI > 2'))) %>%
    as_tibble()%>%
    #rbind our dummy dataframe to ensure all levels are present
    rbind(., dummy_df) %>%
    #define the group_by structure by drought classes, don't drop missing classes
    group_by(drought, .drop = FALSE) %>%
    #compute the median difference
    summarise(bias = median(diff, na.rm = T))

  #reorganize the final data frame and rename data
  temp = temp$bias[1:5]
  export = data.frame(temp) %>% t %>% as.data.frame()
  colnames(export) = c('SPI > 2', '1 > SPI > 2', '1 > SPI > -1', '-1 > SPI > -2', '-2 > SPI') %>% rev
  rownames(export) = NULL
  return(export)
}

#define function to compute bias for all time periods
bias_all = function(x){
  #summarise data
  temp = x %>%
    #define some time meta data
    mutate(month = month(time),
           year = year(time),
           #compute the difference in the historical and contemporary data
           diff = spi_historical - spi_contemporary)%>%
    #filter the data for months of interest
    filter(month %in% c(6,7,8),
           #filter for a 25 year minimum climatology in the contemporary data
           n_contemporary >= 25,
           #filter for a 70 year minimum climatology in the historical data
           n_historical >= 70)%>%
    #drop nas
    tidyr::drop_na() %>%
    #convert to tibble
    as_tibble()%>%
    #compute the median difference in values
    summarise(bias = median(diff, na.rm = T))

  return(temp$bias)
}

# compute average bias for all data together (wet and dry, all D classes together)
bias = lapply(spi_comparison, bias_all) %>%
  unlist()

print(paste0('Stations Lost (%): ', sum(is.na(bias))/length(bias)*100))

#drought classes broken out
drought_class = lapply(spi_comparison, drought_class_bias) %>%
  data.table::rbindlist(.) %>%
  as_tibble()

apply(drought_class, 2, FUN = function(x){sum(is.na(x))/length(x)*100})

#dryness classes broken out
dryness_class = lapply(spi_comparison, dryness_class_bias) %>%
  data.table::rbindlist(.) %>%
  as_tibble()

apply(dryness_class, 2, FUN = function(x){sum(is.na(x))/length(x)*100})

#merge into single dataframe that summarizes all results
valid_stations_filtered = valid_stations %>%
  mutate(`Average Bias` = bias,
         D0 = drought_class$D0,
         D1 = drought_class$D1,
         D2 = drought_class$D2,
         D3 = drought_class$D3,
         D4 = drought_class$D4,
         `-2 > SPI` = dryness_class$`-2 > SPI`,
         `-1 > SPI > -2` = dryness_class$`-1 > SPI > -2`,
         `1 > SPI > -1` = dryness_class$`1 > SPI > -1`,
         `1 > SPI > 2` = dryness_class$`1 > SPI > 2`,
         `SPI > 2` = dryness_class$`SPI > 2`)

print(paste0('Negative Bias during SPI < -2: ', 
             (sum(valid_stations_filtered$`-2 > SPI` < 0, na.rm = T)/ 
                sum(!is.na(valid_stations_filtered$`-2 > SPI`)))*100))

print(paste0('Positive Bias during SPI > 2: ', 
             (sum(valid_stations_filtered$`SPI > 2` > 0, na.rm = T)/ 
                sum(!is.na(valid_stations_filtered$`SPI > 2`)))*100))

#################################################
################ Plot the Results ###############
#################################################

#define plotting parameters
#color ramp
col = colorRampPalette((c('darkred', 'red', 'white', 'blue', 'darkblue')))

#classes to loop through for plotting and kriging
classes = c('Average Bias','D0', 'D1', 'D2', 'D3', 'D4',
            '-2 > SPI', '-1 > SPI > -2', '1 > SPI > -1',
            '1 > SPI > 2', 'SPI > 2')

layman = c('Average Bias','D0', 'D1', 'D2', 'D3', 'D4',
           'Very Dry Conditions', 'Dry Conditions', 'Normal Conditions',
           'Wet Conditions', 'Very Wet Conditions')

for(c in 1:length(classes)){
  #define the temp stations assosiated with the class of interest
  temp_stations = valid_stations_filtered %>%
    tidyr::drop_na(classes[c]) %>%
    st_as_sf

  #first lets plot the point data itself
  pts_plot = ggplot(temp_stations)+
    geom_sf(data = states)+
    geom_sf(aes(color = get(classes[c])))+
    scale_color_gradientn(colours = col(100), breaks = c(-0.5, 0, 0.5), limits = c(-0.5, 0.5),
                          labels = c('-0.5 (Dry Bias)', '0 (No Bias)', '0.5 (Wet Bias)'), name = "",
                          oob = scales::squish, guide = F)+
    theme_bw(base_size = 15)+
    ggtitle(paste0('Daily Summer Bias (',layman[c], ', ', classes[c], ')\n', time_scale[[time_scale_id]], ' Day SPI (June 1 - August 31, 1991-2020)'))+
    theme(legend.position = 'none',
          legend.key.width=unit(2,"cm"),
          plot.title = element_text(hjust = 0.5))

  #krige the reults using autofitting the variogram (package automap)
  #template raster to interpolate over
  template = raster::raster(resolution=c(1/3,1/3),
                            crs = sp::CRS("+init=epsg:4326")) %>%
    raster::crop(., as(states, 'Spatial')) %>%
    raster::rasterToPoints() %>%
    as.data.frame() %>%
    st_as_sf(., coords = c('x', 'y'))
  st_crs(template) = st_crs(4326)

  #fit the variogram
  vgm = autofitVariogram(temp_stations[classes[c]] %>% data.frame() %>% .[,1] ~ 1, as(temp_stations, 'Spatial'))
  
  #krige the variogram results
  krig = krige(temp_stations[classes[c]] %>% data.frame() %>% .[,1] ~ 1, temp_stations, template, model=vgm$var_model) %>%
    st_intersection(states)
  #convert to a point system for tile plotting
  krig_pts = st_coordinates(krig) %>%
    as_tibble() %>%
    mutate(val = krig$var1.pred)
  #define the kriged map plot
  krig_plot = ggplot(krig)+
    geom_tile(data = krig_pts, aes(x = X, y = Y, fill = val))+
    geom_sf(data = states, fill = 'transparent', color = 'black')+
    labs(x = "", y = "")+
    ggtitle(NULL)+
    scale_fill_gradientn(colours = col(100), breaks = c(-0.5, 0, 0.5), limits = c(-0.5, 0.5),
                         labels = c('-0.5 (Dry Bias)', '0 (No Bias)', '0.5 (Wet Bias)'), name = "",
                         oob = scales::squish)+
    theme_bw(base_size = 15)+
    theme(legend.position = 'bottom',
          legend.key.width=unit(2,"cm"),
          legend.text=element_text(size=15))
  #generate the final plot by doing a plot_grid call, first
  #align the plots and shrink the space between them by using an empty plot
  #and reducing the relative height of the middle plot to a negative value
  final = cowplot::plot_grid(pts_plot, NULL, krig_plot, ncol = 1, rel_heights = c(1,-0.17,1), align = 'v')
  #save it out
  ggsave(final, file = paste0('/home/zhoylman/drought-year-sensitivity/figs/moving_window_test/spi_bias_maps_',classes[c],'_',time_scale[[time_scale_id]],'day_timescale_June1-Aug31.png'), width = 7, height = 10, units = 'in')

  #fin
}
    
