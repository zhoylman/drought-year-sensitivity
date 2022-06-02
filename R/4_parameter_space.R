#######################################################################

# Script accomponying Hoylman et al., 2022, Nat. Comm.  
# Drought assessment has been outpaced by climate change: 
# empirical arguments for a paradigm shift
# https://doi.org/10.1038/s41467-022-30316-5

# Author: Dr. Zachary Hoylman
# Contact: zachary.hoylman@umontana.edu

#######################################################################

# This script serves two purposes. Firstly, it is meant to compute the 
# parameter space from which the SPI calculations in this analysis were 
# based upon. This script then plots a random million parameters to show
# the parameter space. Second, we select a random 100 parameters 
# to preform the Monte Carlo simulation for multiple parameters &
# stationary assumptions. These parameters are then saved out for later use. 

#######################################################################

#load some packages
library(tidyverse)
library(ggplot2)

#load in parameter space 
#can be downloaded from zenodo repo
spi_comparison_30 = readRDS(paste0('~/temp', '/spi_comparision_moving_window_with_params_30year_30_days.RDS'))
spi_comparison_60 = readRDS(paste0('~/temp', '/spi_comparision_moving_window_with_params_30year_60_days.RDS'))
spi_comparison_90 = readRDS(paste0('~/temp', '/spi_comparision_moving_window_with_params_30year_90_days.RDS'))
  
# function to process data, filtering for >70 years data for "period-of-record"
# data and > 25 years of data for contempary 30 year moving window data.
process_data = function(x, time_scale){
  temp = x %>% 
    bind_rows() %>%
    filter( n_contemporary >= 25,
            n_historical >= 70) %>%
    mutate(`Timescale` = time_scale)
  return(temp)
}

# process and bind
full = bind_rows(process_data(spi_comparison_30, '30 Day'),
                 process_data(spi_comparison_60, '60 Day'),
                 process_data(spi_comparison_90, '90 Day'))

# set random seed for reporduceability
set.seed(100220)

# randomly generate 1 million indicies to slice data
random_index = runif(1000000, 1, length(full$time)) %>% as.integer()

# slice contempary data based on random indicies
data_contemporary = full[random_index,]%>%
  select(shape_contemporary, rate_contemporary, Timescale) %>%
  rename(Shape = shape_contemporary, Rate = rate_contemporary) %>%
  mutate(ID = "Contemporary")

# slice historical data based on random indicies
data_historical =full[random_index,]%>%
  select(shape_historical, rate_historical, Timescale) %>%
  rename(Shape = shape_historical, Rate = rate_historical) %>%
  mutate(ID = "Period of Record")

# combine datasets
final = bind_rows(data_contemporary, data_historical)

# select random data for the monte carlo simulations
random_montecarlo = final[runif(100,1,length(final$Shape)) %>% as.integer() ,]

# save it out for later use
saveRDS(random_montecarlo, '~/drought-year-sensitivity/data/random_parameters_for_monte_carlo.RDS')

# plot the distrobution space as a 2d density plot seperated by 
# contemapry/ historical and timescale
distrobution_plot = ggplot(final, aes(x = Shape, y = Rate))+
  geom_density_2d_filled(
    aes(fill = ..level..),
    contour_var = "ndensity", 
    breaks = seq(0, 1.0, length.out = 10))+ 
  facet_grid(Timescale~ID)+
  scale_fill_viridis_d(guide = F)+
  theme_bw(base_size = 16)+
  labs(x = 'Shape', y = 'Rate')+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('Gamma Distribution Parameter Space')+
  xlim(0,15)+
  ylim(0.01,0.05)

# save
ggsave(distrobution_plot, file = '~/drought-year-sensitivity/figs/parameter_space/parameter_space_facets.png', height = 7, width = 7, units = 'in')

# plot them all together, with the random points for the monte carlo simulations as well
all_together = ggplot(final, aes(x = Shape, y = Rate))+
  geom_density_2d_filled(
    aes(fill = ..level..),
    contour_var = "ndensity", # normalize to each QBs total passes
    breaks = seq(0, 1.0, length.out = 10))+ # drop the lowest passes
  scale_fill_viridis_d(guide = F)+
  theme_bw(base_size = 16)+
  labs(x = 'Shape', y = 'Rate')+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('Gamma Distribution Parameter Space')+
  xlim(0,15)+
  ylim(0.01,0.05)+
  geom_point(data = random_montecarlo, aes(x = Shape, y = Rate), color = 'black', fill = 'white', shape = 21)

# save
ggsave(all_together, file = '~/drought-year-sensitivity/figs/parameter_space/parameter_space_merged.png', height = 7, width = 7, units = 'in')
        
# fin
