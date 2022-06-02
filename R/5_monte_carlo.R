#######################################################################

# Script accomponying Hoylman et al., 2022, Nat. Comm.  
# Drought assessment has been outpaced by climate change: 
# empirical arguments for a paradigm shift
# https://doi.org/10.1038/s41467-022-30316-5

# Author: Dr. Zachary Hoylman
# Contact: zachary.hoylman@umontana.edu

#######################################################################

# This script preforms the monte carlo simulations presented in the manuscript.
# The monte carlo simulations are conducted by generating a synthetic precipitation distrobution
# of known parameters to evalaute how much data it takes to apporximate the
# known distrobtion. We also use this script to analyze error in the resulting 
# approximated distrobution, CDF and SPI values. This is conducted 4 times 
# using different assumptions, 1. stationary distrobution with known parameters, 
# 2. stationary assumption with many (100) known parameters selected from the observed 
# parameter space (see script 4), 3. single non-stationary distrobution for known 30 year 
# moving window parameters from Clemson Univ., SC and finally 4. non-stationary distrobutions
# for an additional 10 sites using observed 30 year moving window parameters. 

#######################################################################
# load required packages

library(tidyverse)
library(lmomco)
library(cowplot)
library(foreach)
library(doParallel)
library(sf)

#source functions to compute probabilistic CDF and probabilistic parameters
source('~/drought-year-sensitivity/R/funs/gamma_fit_spi.R')

#######################################################################
####################### STATIONARY DISTROBUTION #######################
#######################################################################   

#define proabaility distrobution paramters
shape = 2.5 #alpha
rate = 0.03 #beta = 1/rate

#set seed for reproducability
set.seed(98)

#define monte carlo paramters
n_samples = seq(1,100,1)
n_simulation = 1000

#set up export data frames
export_df_mae = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_df_spi_mae = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_df_shape = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_df_rate = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))

#run monte carlo
#first for loop is per simulation 1:1000
for(i_n_simulations in 1:n_simulation){
  #for each simulation first:
  #generate synthetic data from which samples will be randomly sampled
  data = rgamma(1000, shape, rate)
  #second for loop is per number of samples, 1:100
  for(i_n_samples in n_samples){
    #randomly sample from the synthetic data above
    #number of samples, 1:100 is defined by the for itterator
    temp_data = sample(data, n_samples[i_n_samples])
    
    #CDF comparison
    #compute the true vals with known parameters
    true_vals = cdfgam(temp_data, vec2par(c(shape, 1/rate), 'gam'))
    #compute probabilistic CDF
    probabilistic_vals = gamma_fit_spi(temp_data, 'CDF', return_latest = T)
    #compute the Mean Absolute Error (Prob - True) for the most recent observation
    mae = abs(probabilistic_vals - true_vals[length(true_vals)])

    #SPI comparison
    true_spi_vals = qnorm(true_vals)
    #comupute probabilistic vals from limited climatology
    probabilistic_spi_vals = gamma_fit_spi(temp_data, 'SPI', return_latest = T)
    #compute MAE of the spi vals for the most recent observation
    mae_spi = abs(probabilistic_spi_vals - true_spi_vals[length(true_spi_vals)])

    #compute parameters of probability dist given the limited sample
    params = gamma_fit_spi(temp_data, 'params')
    
    #populate the export dataframes
    export_df_mae[i_n_samples, i_n_simulations] = mae
    export_df_spi_mae[i_n_samples, i_n_simulations] = mae_spi
    
    #if the L-moments couldn't fit a distrobution (AKA only 1 observation)
    #than return an NA, else store the paramters. 
    if(anyNA(params) == T){
      export_df_shape[i_n_samples, i_n_simulations] = NA
      export_df_rate[i_n_samples, i_n_simulations] = NA
    }
    if(anyNA(params) == F){
      export_df_shape[i_n_samples, i_n_simulations] = params$para[1]
      export_df_rate[i_n_samples, i_n_simulations] = 1/params$para[2]
    }
  }
  #print simulation itteration
  print(i_n_simulations)
}

#define summary function to process results
summarize_fun = function(x){
  export = x %>%
    # save number of samples for grouping
    mutate(n_obs = n_samples) %>%
    # pivot longer for grouping
    tidyr::pivot_longer(cols = -c(n_obs)) %>%
    # clean up columns selected
    dplyr::select(-name) %>%
    # group by the number of observations
    group_by(n_obs) %>%
    # sumamrize, median and IQR
    summarise(median = median(value, na.rm = T),
              upper =  quantile(value, 0.75, na.rm = T),
              lower = quantile(value, 0.25, na.rm = T))
  return(export)
}

#summarize monte carlo results
summary_mae = summarize_fun(export_df_mae)
summary_spi_mae = summarize_fun(export_df_spi_mae)
summary_rate = summarize_fun(export_df_rate)
summary_shape = summarize_fun(export_df_shape)

#plot the results (MAE)
plot_mae = ggplot(data = summary_mae, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = 'Number of Observations in Climatology', y = 'CDF Absolute Error')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_segment(data = NULL, aes(x = 30, y = 0, xend = 30, yend = summary_mae[30,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 30, y = summary_mae[30,]$median, xend = 45, yend = summary_mae[30,]$median + .15), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = 0, xend = 60, yend = summary_mae[60,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = summary_mae[60,]$median, xend = 70, yend = summary_mae[60,]$median + .11), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = 0, xend = 90, yend = summary_mae[90,]$median  + .04), linetype = 'dashed', color = 'red')+
  geom_text(data = NULL, aes(x = 46, y = summary_mae[30,]$median + .15,  
                             label = paste0(summary_mae[30,]$median %>% round(., 2), ' [', (summary_mae[30,]$upper - summary_mae[30,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 71, y = summary_mae[60,]$median + .11, 
                             label = paste0(summary_mae[60,]$median %>% round(., 2), ' [', (summary_mae[60,]$upper - summary_mae[60,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 79.5, y = summary_mae[90,]$median + .06, 
                             label = paste0(summary_mae[90,]$median %>% round(., 2), ' [', (summary_mae[90,]$upper - summary_mae[90,]$lower) %>% round(., 2), ']')), hjust = 0)+
  scale_x_continuous(breaks = c(0,30,60,90))+
  geom_line()
  
#plot the results (SPI)
plot_spi_mae = ggplot(data = summary_spi_mae, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = 'Number of Observations in Climatology', y = 'SPI Absolute Error')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_segment(data = NULL, aes(x = 30, y = 0, xend = 30, yend = summary_spi_mae[30,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 30, y = summary_spi_mae[30,]$median, xend = 45, yend = summary_spi_mae[30,]$median + .45), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = 0, xend = 60, yend = summary_spi_mae[60,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = summary_spi_mae[60,]$median, xend = 70, yend = summary_spi_mae[60,]$median + .35), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = 0, xend = 90, yend = summary_spi_mae[90,]$median  + .13), linetype = 'dashed', color = 'red')+
  geom_text(data = NULL, aes(x = 46, y = summary_spi_mae[30,]$median + .45, 
                             label = paste0(summary_spi_mae[30,]$median %>% round(., 2), ' [', (summary_spi_mae[30,]$upper - summary_spi_mae[30,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 71, y = summary_spi_mae[60,]$median + .35, 
                             label = paste0(summary_spi_mae[60,]$median %>% round(., 2), ' [', (summary_spi_mae[60,]$upper - summary_spi_mae[60,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 79.5, y = summary_spi_mae[90,]$median + .18, 
                             label = paste0(summary_spi_mae[90,]$median %>% round(., 2), ' [', (summary_spi_mae[90,]$upper - summary_spi_mae[90,]$lower) %>% round(., 2), ']')), hjust = 0)+
  scale_x_continuous(breaks = c(0,30,60,90))+
  geom_line()

#plot the results (rate parameter)
plot_rate = ggplot(data = summary_rate, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = NULL, y = 'Rate Parameter')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_text(data = NULL, aes(x = 50, y = .17, label = paste0('True Rate Parameter = ', rate)), size = 6)+
  geom_hline(yintercept=rate, linetype="dashed", color = "red")+
  geom_segment(data = NULL, aes(x = 30, y = 0, xend = 30, yend = summary_rate[30,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 30, y = summary_rate[30,]$median, xend = 45, yend = summary_rate[30,]$median + .075), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = 0, xend = 60, yend = summary_rate[60,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = summary_rate[60,]$median, xend = 70, yend = summary_rate[60,]$median + 0.05), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = 0, xend = 90, yend = summary_rate[90,]$median  + .015), linetype = 'dashed', color = 'red')+
  geom_text(data = NULL, aes(x = 46, y = summary_rate[30,]$median + .075, 
                             label = paste0(summary_rate[30,]$median %>% round(., 2), ' [',
                                            (summary_rate[30,]$upper - summary_rate[30,]$lower) %>% round(., 3), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 71, y = summary_rate[60,]$median + 0.055, 
                             label = paste0(summary_rate[30,]$median %>% round(., 2), ' [',
                                             (summary_rate[60,]$upper - summary_rate[60,]$lower) %>% round(., 3), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 79.5, y = summary_rate[90,]$median + .025, 
                             label = paste0(summary_rate[90,]$median %>% round(., 2), ' [', (summary_rate[90,]$upper - summary_rate[90,]$lower) %>% round(., 3), ']')), hjust = 0)+
  scale_x_continuous(breaks = c(0,30,60,90))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#plot the results (shape parameter)
plot_shape = ggplot(data = summary_shape, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = NULL, y = 'Shape Parameter')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_text(data = NULL, aes(x = 50, y = 12.5, label = paste0('True Shape Parameter = ', shape)), size = 6)+
  geom_hline(yintercept=shape, linetype="dashed", color = "red")+
  geom_segment(data = NULL, aes(x = 30, y = 0, xend = 30, yend = summary_shape[30,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 30, y = summary_shape[30,]$median, xend = 45, yend = summary_shape[30,]$median + 5.5), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = 0, xend = 60, yend = summary_shape[60,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = summary_shape[60,]$median, xend = 70, yend = summary_shape[60,]$median + 3.5), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = 0, xend = 90, yend = summary_shape[90,]$median  + 1), linetype = 'dashed', color = 'red')+
  geom_text(data = NULL, aes(x = 46, y = summary_shape[30,]$median + 5.5, 
                             label = paste0(summary_shape[30,]$median %>% round(., 2), ' [',
                                            (summary_shape[30,]$upper - summary_shape[30,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 71, y = summary_shape[60,]$median + 3.8, 
                             label = paste0(summary_shape[30,]$median %>% round(., 2), ' [',
                                            (summary_shape[60,]$upper - summary_shape[60,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 79.5, y = summary_shape[90,]$median + 1.5, 
                             label = paste0(summary_shape[90,]$median %>% round(., 2), ' [', (summary_shape[90,]$upper - summary_shape[90,]$lower) %>% round(., 2), ']')), hjust = 0)+
  scale_x_continuous(breaks = c(0,30,60,90))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#combine to single plot
plot_final = cowplot::plot_grid(plot_rate, plot_shape, NULL, NULL, plot_mae, plot_spi_mae,
                                nrow = 3, ncol = 2, align = 'hv', rel_heights = c(1,-0.18,1))
title = ggdraw() + 
  draw_label("                Stationary Climate\n(Single Parameter Pair, 1000 Iterations)",x = 0,hjust = 0,size = 20) +
  theme(plot.margin = margin(0, 0, 0, 200))

final = plot_grid(title, plot_final,ncol = 1,rel_heights = c(0.1, 1))

final
#save
ggsave(final, file = '~/drought-year-sensitivity/figs/monte_carlo/year_sensitivity_monte_carlo.png', width = 10, height = 7, units = 'in')
#fin single parameter stationary dist

#######################################################################
######## STATIONARY DISTROBUTION Multiple Parameters ##################
#######################################################################

#import randomly selected parameters from script 4
params_space = readRDS('~/drought-year-sensitivity/data/random_parameters_for_monte_carlo.RDS')

#set seed for reproducability
set.seed(98)

#define monte carlo parameters
n_samples = seq(1,100,1)
n_simulation = 1000

#rev up a cluster for parallel computing
cl = makeCluster(detectCores()-2)
#register the cluster for doPar
registerDoParallel(cl)

#define out list for saving paralell output
out = list()

#run simulations, very similar to last nested for loop, but in this case
#add a foreach loop outside for each parameter pair. In this case we are 
#focusing on SPI mean absolute error
out = foreach(i = 1:length(params_space$Shape), .packages = c('lmomco')) %dopar% {
  #for each paraemter pair 
  #define the output dataframe
  export_df_spi_mae = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
  #run monte carlo
  for(i_n_simulations in 1:n_simulation){
    #for each simulation
    #generate random data following a gamma distrobution using foreach parameter pair
    data = rgamma(1000, params_space$Shape[i], params_space$Rate[i])
    for(i_n_samples in n_samples){
      #randomly sample from the synthetic data above
      temp_data = sample(data, n_samples[i_n_samples])
      
      #CDF comparison
      #compute the true vals with known parameters
      true_vals = cdfgam(temp_data, vec2par(c(params_space$Shape[i], 1/params_space$Rate[i]), 'gam'))
      
      #SPI comparison
      true_spi_vals = qnorm(true_vals)
      #comupute probabilistic vals from limited climatology
      probabilistic_spi_vals = gamma_fit_spi(temp_data, 'SPI', return_latest = T)
      #compute MAE of the most recent spi vals
      mae_spi = abs(probabilistic_spi_vals - true_spi_vals[length(true_spi_vals)])
      
      #populate the export dataframes
      export_df_spi_mae[i_n_samples, i_n_simulations] = mae_spi
    }
  }
  export_df_spi_mae
}

stopCluster(cl) 

#save out siulation reults if you wish
#saveRDS(out, '~/temp-drought/stationary_monte_carlo_100_params.RDS')
#out = readRDS('~/temp-drought/stationary_monte_carlo_100_params.RDS')
#add parameter pair ID for grouping to each dataframe in list
for(i in 1:length(out)){
  out[[i]]$param_pair = i
}

#define summary functions for all simulation together
summaries = lapply(out, function(x){return(x %>% mutate(n_obs = n_samples))}) %>%
  bind_rows() %>%
  tidyr::pivot_longer(cols = -c(n_obs, param_pair)) %>%
  dplyr::select(-c(name, param_pair)) %>%
  group_by(n_obs) %>%
  summarise(median = median(value, na.rm = T),
            upper =  quantile(value, 0.75, na.rm = T),
            lower = quantile(value, 0.25, na.rm = T))

#define summary function for each paraemter pair seperated
summaries_each = lapply(out, function(x){return(x %>% mutate(n_obs = n_samples))}) %>%
  bind_rows() %>%
  tidyr::pivot_longer(cols = -c(n_obs, param_pair)) %>%
  dplyr::select(-c(name)) %>%
  group_by(n_obs, param_pair) %>%
  summarise(median = median(value, na.rm = T),
            upper =  quantile(value, 0.75, na.rm = T),
            lower = quantile(value, 0.25, na.rm = T))

#plot the results for the SPI MAE
plot_spi_mae_param_space = ggplot()+
  geom_ribbon(data = summaries, aes(x = n_obs, y = median, ymax = upper, ymin = lower), fill = 'grey70')+
  theme_bw(base_size = 16)+
  geom_line(data = summaries_each, aes(x = n_obs, y = median, colour = as.factor(param_pair)), alpha = 0.3)+
  geom_line(data = summaries, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  scale_colour_discrete(guide = F, limits = rep(1,100))+
  labs(x = 'Number of Observations in Climatology', y = 'SPI Absolute Error')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_segment(data = NULL, aes(x = 30, y = 0, xend = 30, yend = summaries[30,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 30, y = summaries[30,]$median, xend = 45, yend = summaries[30,]$median + .25), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = 0, xend = 60, yend = summaries[60,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = summaries[60,]$median, xend = 70, yend = summaries[60,]$median + .15), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = 0, xend = 90, yend = summaries[90,]$median  + .08), linetype = 'dashed', color = 'red')+
  geom_text(data = NULL, aes(x = 46, y = summaries[30,]$median + .25, 
                             label = paste0(summaries[30,]$median %>% round(., 2), ' [', (summaries[30,]$upper - summaries[30,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 71, y = summaries[60,]$median + .15, 
                             label = paste0(summaries[60,]$median %>% round(., 2), ' [', (summaries[60,]$upper - summaries[60,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 88, y = summaries[90,]$median + .11, 
                             label = paste0(summaries[90,]$median %>% round(., 2), ' [', (summaries[90,]$upper - summaries[90,]$lower) %>% round(., 2), ']')), hjust = 0)+
  scale_x_continuous(breaks = c(0,30,60,90))+
  geom_line()+
  ggtitle('Stationary Climate\n(100 Parameter Pairs, 1000 Iterations per Pair)')+
  theme(plot.title = element_text(hjust = 0.5))

#save them out
ggsave(plot_spi_mae_param_space, file = '~/drought-year-sensitivity/figs/monte_carlo/year_sensitivity_monte_carlo_multiple_paras.png', width = 7, height = 8*.8, units = 'in')

#######################################################################
##################### NON-STATIONARY DISTROBUTION #####################
#######################################################################
#define function to infill missing years rate and shape parameters.
#this is only done if there is not an "observed value".
#the majority of the data is observed, but to ensure that the parameters
#are changing on a realistic time step, (each year) we infill missing data.
#this way paraemters wont jump from say, 1920 to 1930 without accounting for the
#10 years in between. 
spline_fill = function(x){
  export = data.frame(year = seq(min(x$year), max(x$year), by = 1)) %>%
    left_join(., x, by = 'year') %>%
    mutate(spline_rate = predict(smooth.spline(x$year, x$Rate), year)$y,
           spline_shape = predict(smooth.spline(x$year, x$Shape), year)$y,
           Rate = ifelse(is.na(Rate), spline_rate, Rate),
           Shape = ifelse(is.na(Shape), spline_shape, Shape))
  return(export)
}
#set seed for reporducability
set.seed(100)
#import observed data from Clemson Univ. SC to begin with
non_stationary_example = readRDS('~/drought-year-sensitivity/data/params/param_shift_USC00381770_30_days.RDS') %>%
  spline_fill(.)
  
#non-stationary distrobution
#define probability distrobution with shifting parameters (30 year moving window)
#define monte carlo simulation parameters
n_samples = seq(1,length(non_stationary_example$Shape),1)
n_ = length(non_stationary_example$Shape)

#define input matrix
input_matrix = data.frame(x = n_samples, shape = non_stationary_example$Shape, rate = non_stationary_example$Rate)

#define number of simulations in the monte carlo simulation
n_simulation = 1000

#set up export data frames
export_non_stationary_df_mae_spi = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_non_stationary_df_mae_cdf = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))

#non-stationary data
for(i_n_simulations in 1:n_simulation){
  #for each simulation
  print(i_n_simulations)
  #randomly generate the data for the time-varying distrobution 
  data_non_stationary = input_matrix %>%
    mutate(data =  rgamma(x, shape, rate)) %>%
    mutate(time = seq(1:n_))
  
  for(i_n_samples in n_samples){
    #pull data from the end of the distrobution backwards
    #simulates what we would do in the case of a moving window analysis
    #i.e. 2020-1991
    temp_data = data_non_stationary[length(data_non_stationary$data):(length(data_non_stationary$data)-(i_n_samples-1)),]
    #stract the most recent value to compare against 
    recent_slice = temp_data[1,]
    #compute the true value using the known distrobution
    true_vals = cdfgam(recent_slice$data, vec2par(c(recent_slice$shape, 1/recent_slice$rate), 'gam'))
    #CDF comparison
    #compute probabilistic CDF using the random sample of random data
    #a little confusing here but we need to reverse the vector again to return the 
    #latest value (2020). - - This is here and below for the probabilistic SPI vals
    probabilistic_vals = gamma_fit_spi(temp_data$data %>% rev, 'CDF', return_latest = T)
    #compute the Absolute Error (AE) Prob - True
    mae = abs(probabilistic_vals- true_vals)
    #SPI comparison
    true_spi_vals = qnorm(true_vals)
    #comupute probabilistic vals from limited climatology
    probabilistic_spi_vals = gamma_fit_spi(temp_data$data %>% rev, 'SPI', return_latest = T)
    #compute AE of the spi vals
    mae_spi = mean(abs(probabilistic_spi_vals - true_spi_vals))
    #populate the export dataframes
    export_non_stationary_df_mae_cdf[i_n_samples, i_n_simulations] = mae
    export_non_stationary_df_mae_spi[i_n_samples, i_n_simulations] = mae_spi
  }
}

# summarize the results using the priorly defined function (median and IQR)
summary_non_stationary_mae_spi = summarize_fun(export_non_stationary_df_mae_spi)
summary_non_stationary_mae_cdf = summarize_fun(export_non_stationary_df_mae_cdf)

#plot the results
plot_mae_spi = ggplot(data = summary_non_stationary_mae_spi, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = 'Number of Observations in Climatology', y = 'SPI Absolute Error')+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle('Non-Stationary Climate (CLEMSON UNIV, SC)\n(100 Parameter Pairs, 1000 Iterations)')+
  scale_x_continuous(breaks = c(0,30,60,90))+
  geom_segment(data = NULL, aes(x = 30, y = 0, xend = 30, yend = summary_non_stationary_mae_spi[30,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 30, y = summary_non_stationary_mae_spi[30,]$median, xend = 25, yend = summary_non_stationary_mae_spi[30,]$median + .25), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = 0, xend = 60, yend = summary_non_stationary_mae_spi[60,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = summary_non_stationary_mae_spi[60,]$median, xend = 55, yend = summary_non_stationary_mae_spi[60,]$median + .25), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = 0, xend = 90, yend = summary_non_stationary_mae_spi[90,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = summary_non_stationary_mae_spi[90,]$median, xend = 83, yend = summary_non_stationary_mae_spi[90,]$median + .25), linetype = 'dashed', color = 'red')+
  
  geom_text(data = NULL, aes(x = 18, y = summary_non_stationary_mae_spi[30,]$median + .27, 
                             label = paste0(summary_non_stationary_mae_spi[30,]$median %>% round(., 2), ' [', (summary_non_stationary_mae_spi[30,]$upper - summary_non_stationary_mae_spi[30,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 48, y = summary_non_stationary_mae_spi[60,]$median + .27, 
                             label = paste0(summary_non_stationary_mae_spi[60,]$median %>% round(., 2), ' [', (summary_non_stationary_mae_spi[60,]$upper - summary_non_stationary_mae_spi[60,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 76, y = summary_non_stationary_mae_spi[90,]$median + .27, 
                             label = paste0(summary_non_stationary_mae_spi[90,]$median %>% round(., 2), ' [', (summary_non_stationary_mae_spi[90,]$upper - summary_non_stationary_mae_spi[90,]$lower) %>% round(., 2), ']')), hjust = 0)+
  scale_x_continuous(breaks = c(0,30,60,90))+
  geom_line()
  
#visualize
plot_mae_spi

#save the plot
ggsave(plot_mae_spi, file = '~/drought-year-sensitivity/figs/monte_carlo/year_sensitivity_monte_carlo_non_stationary.png', width = 12, height = 7*.8, units = 'in')

#######################################################################
########### NON-STATIONARY DISTROBUTION Multiple Sites ################
#######################################################################
#set seed for reproduceability
set.seed(100)
#non-stationary distrobution for multiple sites
#define probability distrobution with shifting parameters (observed from script 3)
#sites from names
sites = list.files('~/drought-year-sensitivity/data/params') %>%
  substr(., start = 13, stop = 23)
#timescales from names
time_scales = list.files('~/drought-year-sensitivity/data/params') %>%
  substr(., start = 25, stop = 26)
#data
data = list.files('~/drought-year-sensitivity/data/params', full.names = T) %>%
  lapply(., readRDS)
#define monte carlo simulation parameters
n_simulation = 1000

#rev up a cluster for parallel computing
cl = makeCluster(detectCores()-1)
#register the cluster for doPar
registerDoParallel(cl)
#define output list for simulation results
out_non_stationary = list()
#foreach loop for each site 
out_non_stationary = foreach(s = 1:length(sites), .packages = c('lmomco', 'tidyverse')) %dopar% {
  #foreach site, define input data
  input_data = data[[s]] %>% spline_fill(.)
  #define sample information
  n_samples = seq(1,length(input_data$time),1)
  n_ = length(input_data$time)
  #difine data input matrix per site
  input_matrix = data.frame(x = n_samples, shape = input_data$Shape, rate = input_data$Rate)
  #set up export data frames
  export_non_stationary_df_mae_spi = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
  #non-stationary data
  for(i_n_simulations in 1:n_simulation){
    #for each simulation 
    print(i_n_simulations)
    #randomly generate the data for the moving window parameters per simulation 
    data_non_stationary = input_matrix %>%
      mutate(data =  rgamma(x, shape, rate)) %>%
      mutate(time = seq(1:n_))
    for(i_n_samples in n_samples){
      #for each sample itteration
      #pull data from the end of the distrobution backwards
      #simulates what we would do in the case of a moving window analysis
      #i.e. 2020-1991
      temp_data = data_non_stationary[length(data_non_stationary$data):(length(data_non_stationary$data)-(i_n_samples-1)),]
      #extract the most recent obs (2020) for comparison
      recent_slice = temp_data[1,]
      #CDF estimation
      #compute the true CDF values using the know distorbution
      true_vals = cdfgam(recent_slice$data, vec2par(c(recent_slice$shape, 1/recent_slice$rate), 'gam'))
      #SPI comparison
      # compute the known SPI value from the know CDF value
      true_spi_vals = qnorm(true_vals)
      #comupute probabilistic vals from limited climatology
      #a little confusing here but we need to reverse the vector again to return the 
      #latest value (2020).
      probabilistic_spi_vals = gamma_fit_spi(temp_data$data %>% rev, 'SPI', return_latest = T)
      #compute AE of the spi vals
      mae_spi = abs(probabilistic_spi_vals - true_spi_vals[length(true_spi_vals)])
      #populate the export dataframes
      export_non_stationary_df_mae_spi[i_n_samples, i_n_simulations] = mae_spi
    }
  }
  #add some meta data about each site
  export = export_non_stationary_df_mae_spi %>%
    mutate(n_obs = n_samples,
           site = paste0(sites[s]),
           time_scale = time_scales[s])
  #return final data
  export
}

#end cluster
stopCluster(cl)

#define nice special
`%notin%` = Negate(`%in%`)

#define summary function
summaries_non_stationary = out_non_stationary %>%
  bind_rows() %>%
  tidyr::pivot_longer(cols = -c(n_obs, site, time_scale)) %>%
  dplyr::select(-name) %>%
  #remove Clemson Univ. as it is already in previous fig. 
  dplyr::filter(site %notin% c('USC00381770')) %>%
  group_by(n_obs, time_scale) %>%
  summarise(median = median(value, na.rm = T),
            upper =  quantile(value, 0.75, na.rm = T),
            lower = quantile(value, 0.25, na.rm = T)) %>%
  filter(n_obs <= 100) %>%
  mutate(time_scale = paste0(time_scale, ' Days')) %>%
  rename(Timescale = time_scale)
  
#import valid station information for plotting meta
valid_stations = readRDS('~/drought-year-sensitivity/data/valid_stations_70year_summer_baseline.RDS') %>%
  select(id, state, name) %>%
  as_tibble() %>%
  rename(Site = id)

#add in the site information
summaries_non_stationary_sites = out_non_stationary %>%
  bind_rows() %>%
  mutate(time_scale = paste0(time_scale, ' Days')) %>%
  tidyr::pivot_longer(cols = -c(n_obs, site, time_scale)) %>%
  dplyr::select(-name) %>%
  dplyr::filter(site %notin% c('USC00381770')) %>%
  group_by(n_obs, site, time_scale) %>%
  summarise(median = median(value, na.rm = T),
            upper =  quantile(value, 0.75, na.rm = T),
            lower = quantile(value, 0.25, na.rm = T)) %>%
  filter(n_obs <= 100) %>%
  rename(Site = site, Timescale = time_scale) %>% 
  left_join(., valid_stations, by = 'Site') %>%
  mutate(Station = paste0(name, ', ', state))

#save out simulation results if wanted
#saveRDS(summaries_non_stationary_sites, '~/drought-year-sensitivity/data/param_shift_summary_sites.RDS')
#saveRDS(summaries_non_stationary, '~/drought-year-sensitivity/data/param_shift_summary.RDS')
  
#summaries_non_stationary_sites = readRDS('~/drought-year-sensitivity/data/param_shift_summary_sites.RDS')
#summaries_non_stationary = readRDS('~/drought-year-sensitivity/data/param_shift_summary.RDS')

#plot final monte carlo figure
param_shift = ggplot(data = summaries_non_stationary_sites, aes(x = n_obs, y = median, ymax = upper, ymin = lower, color = Station))+
  geom_ribbon(data = summaries_non_stationary, aes(x = n_obs, y = median, ymax = upper, ymin = lower, color = NULL), fill = 'grey70')+
  geom_line(size = 1)+
  #geom_line(data = summaries_non_stationary, aes(x = n_obs, y = median, ymax = upper, ymin = lower, color = 'All', linetype = 'All'))+
  theme_bw(base_size = 20)+
  scale_colour_manual(values = viridis::turbo(10) %>% as.vector(), name = NULL)+
  labs(x = 'Number of Observations in Climatology', y = 'SPI Absolute Error')+
  theme(plot.title = element_text(hjust = 0.5))+
  #ggtitle('Non-Stationary Climate\n(88 Parameter Pairs, 1000 Iterations)')+
  scale_x_continuous(breaks = c(0,30,60,90))+
  theme(legend.position="bottom", 
        legend.box = "vertical", 
        legend.justification = c(0,0),
        legend.margin = margin(t = 10, r = 0, b = 0, l = -5, unit = "pt"))+
  facet_wrap(~Timescale)+ 
  ylim(0,1)+
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.text=element_text(size=12))+
  guides(colour=guide_legend(nrow=2,byrow=TRUE))

#save it out
ggsave(param_shift, file = '~/drought-year-sensitivity/figs/monte_carlo/year_sensitivity_monte_carlo_non_stationary_multiple_sites.png',
       width = 15, height = 6.5, units = 'in', dpi = 300)
#fin
