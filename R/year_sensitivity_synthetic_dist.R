# generate a synthetic distrobution of known parameters
# to evalaute how much data it takes to apporximate the
# known distrobtion
library(tidyverse)
library(lmomco)
library(cowplot)
library(foreach)
library(doParallel)
library(sf)

#source functions to compute probabilistic CDF and probabilistic parameters
source('~/drought-year-sensitivity/R/gamma_fit_spi.R')

#######################################################################
####################### STATIONARY DISTROBUTION #######################
#######################################################################   

#define proabaility distrobution paramters
shape = 2.5 #alpha
rate = 0.03 #beta = 1/rate

#set seed for reproducability
set.seed(98)

#define monte carlo information
n_samples = seq(1,100,1)
n_simulation = 1000

#set up export data frames
export_df_mae = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_df_spi_mae = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_df_shape = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_df_rate = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))

#run monte carlo
for(i_n_simulations in 1:n_simulation){
  #generate synthetic data (synthetic precipitation distrobution based values from Wu et al., 2005)
  data = rgamma(1000, shape, rate)
  for(i_n_samples in n_samples){
    #Stationary data first
    #randomly sample from the synthetic data above
    temp_data = sample(data, n_samples[i_n_samples])
    
    #compute the emperical CDF - This is a sub optimal approach
    #replaced the ECDF with the true vals (below)
    #emperical_model = ecdf(temp_data) 
    #emperical_vals = emperical_model(temp_data)
    
    #CDF comparison
    #compute the true vals with known parameters
    true_vals = cdfgam(temp_data, vec2par(c(shape, 1/rate), 'gam'))
    #compute probabilistic CDF
    probabilistic_vals = gamma_fit_spi(temp_data, 'CDF', return_latest = T)
    #compute the Mean Absolute Error Prob - True
    #mae = mean(abs(probabilistic_vals - true_vals))
    #if you want to only compute error on latest data use:
    mae = abs(probabilistic_vals - true_vals[length(true_vals)])

    #SPI comparison
    true_spi_vals = qnorm(true_vals)
    #comupute probabilistic vals from limited climatology
    probabilistic_spi_vals = gamma_fit_spi(temp_data, 'SPI', return_latest = T)
    #compute MAE of the spi vals
    #mae_spi = mean(abs(probabilistic_spi_vals - true_spi_vals))
    #if you want to only compute error on latest data use:
    mae_spi = abs(probabilistic_spi_vals - true_spi_vals[length(true_spi_vals)])

    #compute parameters of probability dist given the limited sample
    params = gamma_fit_spi(temp_data, 'params')
    
    #populate the export dataframes
    export_df_mae[i_n_samples, i_n_simulations] = mae
    export_df_spi_mae[i_n_samples, i_n_simulations] = mae_spi
    
    if(is.na(params) == T){
      export_df_shape[i_n_samples, i_n_simulations] = NA
      export_df_rate[i_n_samples, i_n_simulations] = NA
    }
    if(is.na(params) == F){
      export_df_shape[i_n_samples, i_n_simulations] = params$para[1]
      export_df_rate[i_n_samples, i_n_simulations] = 1/params$para[2]
    }
  }
  #print simulation itteration
  print(i_n_simulations)
}

#define summary function
summarize_fun = function(x){
  export = x %>%
    mutate(n_obs = n_samples) %>%
    tidyr::pivot_longer(cols = -c(n_obs)) %>%
    dplyr::select(-name) %>%
    group_by(n_obs) %>%
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

#plot the results
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
ggsave(final, file = '/home/zhoylman/drought-year-sensitivity/figs/year_sensitivity_monte_carlo.png', width = 10, height = 7, units = 'in')
#fin single parameter stationary dist

#######################################################################
######## STATIONARY DISTROBUTION Multiple Parameters ##################
#######################################################################

params_space = readRDS('/home/zhoylman/drought-year-sensitivity/data/random_parameters_for_monte_carlo.RDS')

#set seed for reproducability
set.seed(98)

#define monte carlo information
n_samples = seq(1,100,1)
n_simulation = 1000

#rev up a cluster for parallel computing
cl = makeCluster(detectCores()-2)
#register the cluster for doPar
registerDoParallel(cl)

out = list()

out = foreach(i = 1:length(params_space$Shape), .packages = c('lmomco')) %dopar% {
  export_df_spi_mae = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
  #run monte carlo
  for(i_n_simulations in 1:n_simulation){
    #set up export data frames
    data = rgamma(1000, params_space$Shape[i], params_space$Rate[i])
    for(i_n_samples in n_samples){
      #Stationary data first
      #randomly sample from the synthetic data above
      temp_data = sample(data, n_samples[i_n_samples])
      
      #CDF comparison
      #compute the true vals with known parameters
      true_vals = cdfgam(temp_data, vec2par(c(params_space$Shape[i], 1/params_space$Rate[i]), 'gam'))
      
      #SPI comparison
      true_spi_vals = qnorm(true_vals)
      #comupute probabilistic vals from limited climatology
      probabilistic_spi_vals = gamma_fit_spi(temp_data, 'SPI', return_latest = T)
      #compute MAE of the spi vals
      #mae_spi = mean(abs(probabilistic_spi_vals - true_spi_vals))
      #if you want to only compute error on latest data use:
      mae_spi = abs(probabilistic_spi_vals - true_spi_vals[length(true_spi_vals)])
      
      #populate the export dataframes
      export_df_spi_mae[i_n_samples, i_n_simulations] = mae_spi
    }
  }
  export_df_spi_mae
}

stopCluster(cl) 

#saveRDS(out, '/home/zhoylman/temp/stationary_monte_carlo_100_params.RDS')
#out = readRDS('/home/zhoylman/temp/stationary_monte_carlo_100_params.RDS')
for(i in 1:length(out)){
  out[[i]]$param_pair = i
}

summaries = lapply(out, function(x){return(x %>% mutate(n_obs = n_samples))}) %>%
  bind_rows() %>%
  tidyr::pivot_longer(cols = -c(n_obs, param_pair)) %>%
  dplyr::select(-c(name, param_pair)) %>%
  group_by(n_obs) %>%
  summarise(median = median(value, na.rm = T),
            upper =  quantile(value, 0.75, na.rm = T),
            lower = quantile(value, 0.25, na.rm = T))

summaries_each = lapply(out, function(x){return(x %>% mutate(n_obs = n_samples))}) %>%
  bind_rows() %>%
  tidyr::pivot_longer(cols = -c(n_obs, param_pair)) %>%
  dplyr::select(-c(name)) %>%
  group_by(n_obs, param_pair) %>%
  summarise(median = median(value, na.rm = T),
            upper =  quantile(value, 0.75, na.rm = T),
            lower = quantile(value, 0.25, na.rm = T))
  

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

ggsave(plot_spi_mae_param_space, file = '/home/zhoylman/drought-year-sensitivity/figs/year_sensitivity_monte_carlo_multiple_paras.png', width = 7, height = 7*.8, units = 'in')

#######################################################################
##################### NON-STATIONARY DISTROBUTION #####################
#######################################################################

spline_fill = function(x){
  export = data.frame(year = seq(min(x$year), max(x$year), by = 1)) %>%
    left_join(., x, by = 'year') %>%
    mutate(spline_rate = predict(smooth.spline(x$year, x$Rate), year)$y,
           spline_shape = predict(smooth.spline(x$year, x$Shape), year)$y,
           Rate = ifelse(is.na(Rate), spline_rate, Rate),
           Shape = ifelse(is.na(Shape), spline_shape, Shape))
  return(export)
}

set.seed(100)
non_stationary_example = readRDS('/home/zhoylman/drought-year-sensitivity/data/params/param_shift_USC00381770_30_days.RDS') %>%
  spline_fill(.)
#non-stationary distrobution
#define probability distrobution with shifting parameters (10 decade chunks)
n_samples = seq(1,length(non_stationary_example$Shape),1)
n_ = length(non_stationary_example$Shape)

input_matrix = data.frame(x = n_samples, shape = non_stationary_example$Shape, rate = non_stationary_example$Rate)

n_simulation = 1000

#set up export data frames
export_non_stationary_df_mae_spi = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_non_stationary_df_mae_cdf = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))

#non-stationary data
for(i_n_simulations in 1:n_simulation){
  print(i_n_simulations)
  #randomly generate the distrobtuion each simulation
  data_non_stationary = input_matrix %>%
    mutate(data =  rgamma(x, shape, rate)) %>%
    mutate(time = seq(1:n_))
  
  for(i_n_samples in n_samples){
    #pull data from the end of the distrobution backwards
    #simulates what we would do in the case of a moving window analysis
    #to test parameter stability 
    temp_data = data_non_stationary[length(data_non_stationary$data):(length(data_non_stationary$data)-(i_n_samples-1)),]

    recent_slice = temp_data[1,]
    true_vals = cdfgam(recent_slice$data, vec2par(c(recent_slice$shape, 1/recent_slice$rate), 'gam'))
    #compute probabilistic CDF
    probabilistic_vals = gamma_fit_spi(temp_data$data %>% rev, 'CDF', return_latest = T)
    #compute the Mean Absolute Error Prob - True
    mae = abs(probabilistic_vals- true_vals)
    #if you want to only compute error on latest data use:
    #mae = mean(abs(probabilistic_vals - true_vals[length(true_vals)]))
    
    #SPI comparison
    true_spi_vals = qnorm(true_vals)
    #comupute probabilistic vals from limited climatology
    probabilistic_spi_vals = gamma_fit_spi(temp_data$data %>% rev, 'SPI', return_latest = T)
    #compute MAE of the spi vals
    mae_spi = mean(abs(probabilistic_spi_vals - true_spi_vals))
    #if you want to only compute error on latest data use:
    #mae_spi = mean(abs(probabilistic_spi_vals - true_spi_vals[length(true_spi_vals)]))
    
    #populate the export dataframes
    export_non_stationary_df_mae_cdf[i_n_samples, i_n_simulations] = mae
    export_non_stationary_df_mae_spi[i_n_samples, i_n_simulations] = mae_spi
  }
}

summary_non_stationary_mae_spi = summarize_fun(export_non_stationary_df_mae_spi)
summary_non_stationary_mae_cdf = summarize_fun(export_non_stationary_df_mae_cdf)

#plot the results
plot_mae_spi = ggplot(data = summary_non_stationary_mae_spi, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = 'Number of Observations in Climatology', y = 'SPI Absolute Error')+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle('Non-Stationary Climate (CLEMSON UNIV, SC)\n(99 Parameter Pairs, 1000 Iterations)')+
  scale_x_continuous(breaks = c(0,30,60,90))+
  geom_segment(data = NULL, aes(x = 30, y = 0, xend = 30, yend = summary_non_stationary_mae_spi[30,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 30, y = summary_non_stationary_mae_spi[30,]$median, xend = 25, yend = summary_non_stationary_mae_spi[30,]$median + .25), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = 0, xend = 60, yend = summary_non_stationary_mae_spi[60,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = summary_non_stationary_mae_spi[60,]$median, xend = 55, yend = summary_non_stationary_mae_spi[60,]$median + .25), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = 0, xend = 90, yend = summary_non_stationary_mae_spi[90,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = summary_non_stationary_mae_spi[90,]$median, xend = 83, yend = summary_non_stationary_mae_spi[90,]$median + .25), linetype = 'dashed', color = 'red')+
  
  geom_text(data = NULL, aes(x = 12, y = summary_non_stationary_mae_spi[30,]$median + .27, 
                             label = paste0(summary_non_stationary_mae_spi[30,]$median %>% round(., 2), ' [', (summary_non_stationary_mae_spi[30,]$upper - summary_non_stationary_mae_spi[30,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 42, y = summary_non_stationary_mae_spi[60,]$median + .27, 
                             label = paste0(summary_non_stationary_mae_spi[60,]$median %>% round(., 2), ' [', (summary_non_stationary_mae_spi[60,]$upper - summary_non_stationary_mae_spi[60,]$lower) %>% round(., 2), ']')), hjust = 0)+
  geom_text(data = NULL, aes(x = 70, y = summary_non_stationary_mae_spi[90,]$median + .27, 
                             label = paste0(summary_non_stationary_mae_spi[90,]$median %>% round(., 2), ' [', (summary_non_stationary_mae_spi[90,]$upper - summary_non_stationary_mae_spi[90,]$lower) %>% round(., 2), ']')), hjust = 0)+
  scale_x_continuous(breaks = c(0,30,60,90))+
  geom_line()
  
plot_mae_spi

ggsave(plot_mae_spi, file = '/home/zhoylman/drought-year-sensitivity/figs/year_sensitivity_monte_carlo_non_stationary.png', width = 7, height = 7*.8, units = 'in')


#######################################################################
########### NON-STATIONARY DISTROBUTION Multiple Sites ################
#######################################################################
set.seed(100)
#non-stationary distrobution
#define probability distrobution with shifting parameters (10 decade chunks)
sites = list.files('/home/zhoylman/drought-year-sensitivity/data/params') %>%
  substr(., start = 13, stop = 23)

time_scales = list.files('/home/zhoylman/drought-year-sensitivity/data/params') %>%
  substr(., start = 25, stop = 26)

data = list.files('/home/zhoylman/drought-year-sensitivity/data/params', full.names = T) %>%
  lapply(., readRDS)

n_simulation = 1000

#rev up a cluster for parallel computing
cl = makeCluster(detectCores()-2)
#register the cluster for doPar
registerDoParallel(cl)

out_non_stationary = list()

out_non_stationary = foreach(s = 1:length(sites), .packages = c('lmomco', 'tidyverse')) %dopar% {
  input_data = data[[s]] %>% spline_fill(.)
  
  n_samples = seq(1,length(input_data$time),1)
  n_ = length(input_data$time)
  
  input_matrix = data.frame(x = n_samples, shape = input_data$Shape, rate = input_data$Rate)

  #set up export data frames
  export_non_stationary_df_mae_spi = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))

  #non-stationary data
  for(i_n_simulations in 1:n_simulation){
    print(i_n_simulations)
    #randomly generate the distrobtuion each simulation 
    data_non_stationary = input_matrix %>%
      mutate(data =  rgamma(x, shape, rate)) %>%
      mutate(time = seq(1:n_))
    
    for(i_n_samples in n_samples){
      #pull data from the end of the distrobution backwards
      #simulates what we would do in the case of a moving window analysis
      #to test parameter stability 
      temp_data = data_non_stationary[length(data_non_stationary$data):(length(data_non_stationary$data)-(i_n_samples-1)),]
      
      recent_slice = temp_data[1,]
      true_vals = cdfgam(recent_slice$data, vec2par(c(recent_slice$shape, 1/recent_slice$rate), 'gam'))
      
      #SPI comparison
      true_spi_vals = qnorm(true_vals)
      #comupute probabilistic vals from limited climatology
      probabilistic_spi_vals = gamma_fit_spi(temp_data$data %>% rev, 'SPI', return_latest = T)
      #compute MAE of the spi vals
      #mae_spi = mean(abs(probabilistic_spi_vals - true_spi_vals))
      #if you want to only compute error on latest data use:
      mae_spi = abs(probabilistic_spi_vals - true_spi_vals[length(true_spi_vals)])
      
      #populate the export dataframes
      export_non_stationary_df_mae_spi[i_n_samples, i_n_simulations] = mae_spi
    }
  }
  export = export_non_stationary_df_mae_spi %>%
    mutate(n_obs = n_samples,
           site = paste0(sites[s]),
           time_scale = time_scales[s])
  export
}

stopCluster(cl)

summaries_non_stationary = out_non_stationary %>%
  bind_rows() %>%
  tidyr::pivot_longer(cols = -c(n_obs, site, time_scale)) %>%
  dplyr::select(-name) %>%
  group_by(n_obs, time_scale) %>%
  summarise(median = median(value, na.rm = T),
            upper =  quantile(value, 0.75, na.rm = T),
            lower = quantile(value, 0.25, na.rm = T)) %>%
  filter(n_obs <= 100) %>%
  mutate(time_scale = paste0(time_scale, ' Days')) %>%
  rename(Timescale = time_scale)
  

valid_stations = readRDS('/home/zhoylman/drought-year-sensitivity/data/valid_stations_70year_summer_baseline.RDS') %>%
  select(id, state, name) %>%
  as_tibble() %>%
  rename(Site = id)

summaries_non_stationary_sites = out_non_stationary %>%
  bind_rows() %>%
  mutate(time_scale = paste0(time_scale, ' Days')) %>%
  tidyr::pivot_longer(cols = -c(n_obs, site, time_scale)) %>%
  dplyr::select(-name) %>%
  group_by(n_obs, site, time_scale) %>%
  summarise(median = median(value, na.rm = T),
            upper =  quantile(value, 0.75, na.rm = T),
            lower = quantile(value, 0.25, na.rm = T)) %>%
  filter(n_obs <= 100) %>%
  rename(Site = site, Timescale = time_scale) %>% 
  left_join(., valid_stations, by = 'Site') %>%
  mutate(Station = paste0(name, ', ', state))

saveRDS(summaries_non_stationary_sites, '/home/zhoylman/drought-year-sensitivity/data/param_shift_summary_sites.RDS')
saveRDS(summaries_non_stationary, '/home/zhoylman/drought-year-sensitivity/data/param_shift_summary.RDS')
  
summaries_non_stationary_sites = readRDS('/home/zhoylman/drought-year-sensitivity/data/param_shift_summary_sites.RDS')
summaries_non_stationary = readRDS('/home/zhoylman/drought-year-sensitivity/data/param_shift_summary.RDS')

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
  theme(legend.position="top", 
        legend.box = "vertical", 
        legend.justification = c(0,0),
        legend.margin = margin(t = 10, r = 0, b = 0, l = -5, unit = "pt"))+
  facet_wrap(~Timescale)+ 
  ylim(0,1)+
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.text=element_text(size=10))+
  guides(colour=guide_legend(nrow=2,byrow=TRUE))

ggsave(param_shift, file = '/home/zhoylman/drought-year-sensitivity/figs/year_sensitivity_monte_carlo_non_stationary_multiple_sites.png',
       width = 13, height = 6.5, units = 'in', dpi = 300)
  