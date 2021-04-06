# generate a synthetic distrobution of known parameters
# to evalaute how much data it takes to apporximate the
# known distrobtion
library(tidyverse)
library(lmomco)
library(magrittr)

#source functions to compute probabilistic CDF and probabilistic parameters
source('R/gamma_fit_spi.R')

# Plan: 1000 replicates
# - 100 years, stationary data
# - 100 years, evolving alpha
# - 100 years, evolving beta
# - 100 years, co-evolution
# Assumptions:
# - No autocorrelation

calc_theoretical_spi <- 
  Vectorize(
    function(value, alpha, beta){
      qnorm(cdfgam(value, vec2par(c(alpha, beta), 'gam')))
    }
  )

calc_retrospective_spi <- 
    function(x){
      x %>%
        magrittr::set_names(1:length(x)) %>%
        purrr::imap_dbl(~gamma_fit_spi(x[1:.y], 
                                       export_opts = 'SPI',
                                       return_latest = T))
    }

calc_climatological_spi <- 
  function(x){
    x %>%
      gamma_fit_spi(export_opts = 'SPI',
                                     return_latest = F)
  }

test <-
  list(
    Stationary = 
      1:99 %>%
      purrr::map_dfr(
        ~tibble::tibble(Run = .x,
                        Timestep = 1:100,
                        Value = rgamma(n = 100, 
                                       shape = 40,
                                       rate = 0.7),
                        alpha = 40,
                        beta = 1/0.7)
      ),
    `Changing Alpha` = 
      1:99 %>%
      purrr::map_dfr(
        ~tibble::tibble(Run = .x,
                        Timestep = 1:100,
                        Value = rgamma(n = 100, 
                                       shape = 40 + ((1:100) - 1),
                                       rate = 0.7),
                        alpha = 40 + (1:100),
                        beta = 1/0.7)
      )
  ) %>%
  dplyr::bind_rows(.id = "Run Type") %>%
  dplyr::group_by(`Run Type`, Run) %>%
  dplyr::mutate(`Actual SPI` = calc_theoretical_spi(value = Value, 
                                                    alpha = alpha, 
                                                    beta = beta),
                `Retrospective SPI` = calc_retrospective_spi(Value),
                `Climatological SPI` = calc_climatological_spi(Value)) %>%
  dplyr::ungroup()


test %>%
  tidyr::pivot_longer(`Retrospective SPI`:`Climatological SPI`, names_to = "Empirical SPI Type",
                      values_to = "Empirical SPI") %>%
  dplyr::group_by(`Run Type`,
                  Timestep,
                  `Empirical SPI Type`) %>%
  # dplyr::mutate(Error = abs(`Empirical SPI` - `Actual SPI` )) %>%
  dplyr::summarise(Error = (`Empirical SPI` - `Actual SPI`) %>%
                     abs() %>%
                     mean(na.rm = TRUE)) %>%
  ggplot(aes(x = Timestep,
             y = Error)) +
  # geom_smooth() +
  geom_line() +
  facet_grid(`Run Type` ~ `Empirical SPI Type`)
  
  





#######################################################################
####################### STATIONARY DISTRIBUTION #######################
#######################################################################

#define probability distribution parameters
shape = 40 #alpha
rate = 0.7 #beta = 1/rate

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
    probabilistic_vals = gamma_fit_spi(temp_data, 'CDF')
    #compute the Mean Absolute Error Prob - True
    mae = mean(abs(probabilistic_vals - true_vals))
    
    #SPI comparison
    true_spi_vals = qnorm(true_vals)
    #comupute probabilistic vals from limited climatology
    probabilistic_spi_vals = gamma_fit_spi(temp_data, 'SPI')
    #compute MAE of the spi vals
    mae_spi = mean(abs(probabilistic_spi_vals - true_spi_vals))
    
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
  labs(x = 'Number of Observations in "Climatology"', y = 'CDF Mean Absolute Error')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_segment(data = NULL, aes(x = 30, y = 0, xend = 30, yend = summary_mae[30,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 30, y = summary_mae[30,]$median, xend = 45, yend = summary_mae[30,]$median + .08), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = 0, xend = 60, yend = summary_mae[60,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = summary_mae[60,]$median, xend = 70, yend = summary_mae[60,]$median + .05), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = 0, xend = 90, yend = summary_mae[90,]$median  + .02), linetype = 'dashed', color = 'red')+
  geom_text(data = NULL, aes(x = 46, y = summary_mae[30,]$median + .08,  
                             label = paste0(summary_mae[30,]$median %>% round(., 2), ' ± ', (summary_mae[30,]$upper - summary_mae[30,]$lower) %>% round(., 2))), hjust = 0)+
  geom_text(data = NULL, aes(x = 71, y = summary_mae[60,]$median + .05, 
                             label = paste0(summary_mae[60,]$median %>% round(., 2), ' ± ', (summary_mae[60,]$upper - summary_mae[60,]$lower) %>% round(., 2))), hjust = 0)+
  geom_text(data = NULL, aes(x = 88, y = summary_mae[90,]$median + .027, 
                             label = paste0(summary_mae[90,]$median %>% round(., 2), ' ± ', (summary_mae[90,]$upper - summary_mae[90,]$lower) %>% round(., 2))), hjust = 0)+
  scale_x_continuous(breaks = c(0,30,60,90))+
  geom_line()

plot_spi_mae = ggplot(data = summary_spi_mae, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = 'Number of Observations in "Climatology"', y = 'SPI Mean Absolute Error')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_segment(data = NULL, aes(x = 30, y = 0, xend = 30, yend = summary_spi_mae[30,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 30, y = summary_spi_mae[30,]$median, xend = 45, yend = summary_spi_mae[30,]$median + .25), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = 0, xend = 60, yend = summary_spi_mae[60,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = summary_spi_mae[60,]$median, xend = 70, yend = summary_spi_mae[60,]$median + .15), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = 0, xend = 90, yend = summary_spi_mae[90,]$median  + .05), linetype = 'dashed', color = 'red')+
  geom_text(data = NULL, aes(x = 46, y = summary_spi_mae[30,]$median + .25, 
                             label = paste0(summary_spi_mae[30,]$median %>% round(., 2), ' ± ', (summary_spi_mae[30,]$upper - summary_spi_mae[30,]$lower) %>% round(., 2))), hjust = 0)+
  geom_text(data = NULL, aes(x = 71, y = summary_spi_mae[60,]$median + .15, 
                             label = paste0(summary_spi_mae[60,]$median %>% round(., 2), ' ± ', (summary_spi_mae[60,]$upper - summary_spi_mae[60,]$lower) %>% round(., 2))), hjust = 0)+
  geom_text(data = NULL, aes(x = 88, y = summary_spi_mae[90,]$median + .08, 
                             label = paste0(summary_spi_mae[90,]$median %>% round(., 2), ' ± ', (summary_spi_mae[90,]$upper - summary_spi_mae[90,]$lower) %>% round(., 2))), hjust = 0)+
  scale_x_continuous(breaks = c(0,30,60,90))+
  geom_line()

plot_rate = ggplot(data = summary_rate, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = NULL, y = 'Rate Parameter')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_text(data = NULL, aes(x = 50, y = 3.8, label = 'True Rate Parameter = 0.7'), size = 6)+
  geom_hline(yintercept=0.7, linetype="dashed", color = "red")+
  geom_segment(data = NULL, aes(x = 30, y = 0, xend = 30, yend = summary_rate[30,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 30, y = summary_rate[30,]$median, xend = 45, yend = summary_rate[30,]$median + 1.25), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = 0, xend = 60, yend = summary_rate[60,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = summary_rate[60,]$median, xend = 70, yend = summary_rate[60,]$median + 0.75), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = 0, xend = 90, yend = summary_rate[90,]$median  + .15), linetype = 'dashed', color = 'red')+
  geom_text(data = NULL, aes(x = 46, y = summary_rate[30,]$median + 1.25, 
                             label = paste0(summary_rate[30,]$median %>% round(., 2), ' ± ',
                                            (summary_rate[30,]$upper - summary_rate[30,]$lower) %>% round(., 2))), hjust = 0)+
  geom_text(data = NULL, aes(x = 71, y = summary_rate[60,]$median + 0.75, 
                             label = paste0(summary_rate[30,]$median %>% round(., 2), ' ± ',
                                            (summary_rate[60,]$upper - summary_rate[60,]$lower) %>% round(., 2))), hjust = 0)+
  geom_text(data = NULL, aes(x = 88, y = summary_rate[90,]$median + .28, 
                             label = paste0(summary_rate[90,]$median %>% round(., 2), ' ± ', (summary_rate[90,]$upper - summary_rate[90,]$lower) %>% round(., 2))), hjust = 0)+
  scale_x_continuous(breaks = c(0,30,60,90))

plot_shape = ggplot(data = summary_shape, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = NULL, y = 'Shape Parameter')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_text(data = NULL, aes(x = 50, y = 225, label = 'True Shape Parameter = 40'), size = 6)+
  geom_hline(yintercept=40, linetype="dashed", color = "red")+
  geom_segment(data = NULL, aes(x = 30, y = 0, xend = 30, yend = summary_shape[30,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 30, y = summary_shape[30,]$median, xend = 45, yend = summary_shape[30,]$median + 75), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = 0, xend = 60, yend = summary_shape[60,]$median), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 60, y = summary_shape[60,]$median, xend = 70, yend = summary_shape[60,]$median + 50), linetype = 'dashed', color = 'red')+
  geom_segment(data = NULL, aes(x = 90, y = 0, xend = 90, yend = summary_shape[90,]$median  + 10), linetype = 'dashed', color = 'red')+
  geom_text(data = NULL, aes(x = 46, y = summary_shape[30,]$median + 75, 
                             label = paste0(summary_shape[30,]$median %>% round(., 2), ' ± ',
                                            (summary_shape[30,]$upper - summary_shape[30,]$lower) %>% round(., 2))), hjust = 0)+
  geom_text(data = NULL, aes(x = 71, y = summary_shape[60,]$median + 50, 
                             label = paste0(summary_shape[30,]$median %>% round(., 2), ' ± ',
                                            (summary_shape[60,]$upper - summary_shape[60,]$lower) %>% round(., 2))), hjust = 0)+
  geom_text(data = NULL, aes(x = 87, y = summary_shape[90,]$median + 20, 
                             label = paste0(summary_shape[90,]$median %>% round(., 2), ' ± ', (summary_shape[90,]$upper - summary_shape[90,]$lower) %>% round(., 2))), hjust = 0)+
  scale_x_continuous(breaks = c(0,30,60,90))

#combine to single plot
plot_final = cowplot::plot_grid(plot_rate, plot_shape,plot_mae, plot_spi_mae,
                                nrow = 2, ncol = 2, align = 'hv')
title = ggdraw() + 
  draw_label("Stationary Monte Carlo Simulation (1000 iterations)", fontface = 'bold',x = 0,hjust = 0,size = 20) +
  theme(plot.margin = margin(0, 0, 0, 300))

final = plot_grid(title, plot_final,ncol = 1,rel_heights = c(0.1, 1))

final
#save
ggsave(final, file = '/home/zhoylman/drought-year-sensitivity/figs/year_sensitivity_monte_carlo.png', width = 13, height = 9, units = 'in')
#fin stationary dist




# coga -> for defining gamma distrobutions

# likely wont do the below analysis.... saving for now. 
######################################################################################################################
######################################################################################################################
######################################################################################################################

#######################################################################
##################### NON-STATIONARY DISTROBUTION #####################
#######################################################################

set.seed(100)
#non-stationary distrobution
#define probability distrobution with shifting parameters (10 decade chunks)
n_ = rep(1,100)
n_samples = seq(1,100,1)
shape_non_stationary = seq(40,49, length.out = 100) # alpha
rate_non_stationary = seq(0.6,3, length.out = 100) # rate

input_matrix = data.frame(x = n_, shape = shape_non_stationary, rate = rate_non_stationary)

#n_simulation = 100

#set up export data frames
export_non_stationary_df_mae_spi = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_non_stationary_df_mae_cdf = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_non_stationary_df_shape = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_non_stationary_df_rate = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_non_stationary_rate_shape_full = data.frame(matrix(nrow = n_simulation, ncol = 2))

#non-stationary data
for(i_n_simulations in 1:n_simulation){
  print(i_n_simulations)
  #randomly generate the distrobtuion each simulation
  data_non_stationary = input_matrix %>%
    mutate(data =  rgamma(x, shape, rate)) %>%
    mutate(time = seq(1:100))
  export_non_stationary_rate_shape_full[i_n_simulations,1] = gamma_fit_spi(data_non_stationary$data, 'params')$para[1]
  export_non_stationary_rate_shape_full[i_n_simulations,2] = gamma_fit_spi(data_non_stationary$data, 'params')$para[2]
  
  for(i_n_samples in n_samples){
    #pull data from the end of teh distrobution backwards
    #simulates what we would do in the case of a moving window analysis
    #to test parameter stability 
    temp_data = data_non_stationary[length(data_non_stationary$data):(length(data_non_stationary$data)-(i_n_samples-1)),] %>%
      mutate(contemporary_vals_spi = gamma_fit_spi(data, 'SPI'),
             contemporary_vals_cdf = gamma_fit_spi(data, 'CDF'))
    
    #here we will simulate the assumptions of Wu and Guttman
    historical_values = data_non_stationary %>%
      mutate(historical_values_spi = gamma_fit_spi(data, 'SPI'),
             historical_vals_cdf = gamma_fit_spi(data, 'CDF'))
    #join data
    joined = left_join(historical_values, temp_data, by = 'time') %>%
      drop_na() %>%
      mutate(diff_spi = abs(contemporary_vals_spi - historical_values_spi),
             diff_cdf = abs(contemporary_vals_cdf - historical_vals_cdf))
    
    mae_spi = mean(joined$diff_spi)
    mae_cdf = mean(joined$diff_cdf)
    
    #compute parameters of probability dist given the limited sample
    params = gamma_fit_spi(temp_data$data, 'params')
    #populate the export dataframes
    export_non_stationary_df_mae_spi[i_n_samples, i_n_simulations] = mae_spi
    export_non_stationary_df_mae_cdf[i_n_samples, i_n_simulations] = mae_cdf
    if(is.na(params) == T){
      export_non_stationary_df_shape[i_n_samples, i_n_simulations] = NA
      export_non_stationary_df_rate[i_n_samples, i_n_simulations] = NA
    }
    if(is.na(params) == F){
      export_non_stationary_df_shape[i_n_samples, i_n_simulations] = params$para[1]
      export_non_stationary_df_rate[i_n_samples, i_n_simulations] = 1/params$para[2]
    }
  }
}

summary_non_stationary_mae_spi = summarize_fun(export_non_stationary_df_mae_spi)
summary_non_stationary_mae_cdf = summarize_fun(export_non_stationary_df_mae_cdf)
summary_non_stationary_rate = summarize_fun(export_non_stationary_df_rate)
summary_non_stationary_shape = summarize_fun(export_non_stationary_df_shape)

#plot the results
plot_mae_spi = ggplot(data = summary_non_stationary_mae_spi, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(data = summary_spi_mae, aes(x = n_obs, y = median, ymax = upper, ymin = lower), fill = 'red', alpha = .20)+
  geom_line(data = summary_spi_mae, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = 'Number of Observations in "Climatology"', y = 'SPI Mean Absolute Error')+
  theme(plot.title = element_text(hjust = 0.5))

plot_mae_spi

#plot the results
plot_mae_cdf = ggplot(data = summary_non_stationary_mae_cdf, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(data = export_df_mae, aes(x = n_obs, y = median, ymax = upper, ymin = lower), fill = 'red', alpha = .20)+
  geom_line(data = export_df_mae, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = 'Number of Observations in "Climatology"', y = 'CDF Mean Absolute Error')+
  theme(plot.title = element_text(hjust = 0.5))

plot_mae_cdf

plot_rate = ggplot(data = summary_non_stationary_rate, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = NULL, y = 'Rate Parameter')+
  ggtitle('Monte Carlo Simulation (1000 simulations)')+
  theme(plot.title = element_text(hjust = 0.5))
#geom_text(data = NULL, aes(x = 80, y = 5, label = 'True Rate Parameter = 0.7'))+
#geom_hline(yintercept=0.7, linetype="dashed", color = "red")

plot_shape = ggplot(data = summary_non_stationary_shape, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = NULL, y = 'Shape Parameter')+
  theme(plot.title = element_text(hjust = 0.5))
#geom_text(data = NULL, aes(x = 80, y = 280, label = 'True Shape Parameter = 40'))+
#geom_hline(yintercept=40, linetype="dashed", color = "red")

#combine to single plot
plot_final = cowplot::plot_grid(plot_rate, plot_shape, plot_mae_cdf, plot_mae_spi,
                                nrow = 2, ncol = 2, align = 'hv')
title = ggdraw() + 
  draw_label("Non-Stationary Monte Carlo Simulation (1000 iterations)", fontface = 'bold',x = 0,hjust = 0,size = 20) +
  theme(plot.margin = margin(0, 0, 0, 300))

final = plot_grid(title, plot_final,ncol = 1,rel_heights = c(0.1, 1))

final
#save
ggsave(plot_final, file = '/home/zhoylman/drought-year-sensitivity/figs/year_sensitivity_monte_carlo_non_stationary.png', width = 6, height = 8, units = 'in')
#fin stationary dist




library(coga)
library(lmomco)
#define some params for the gamma function
shape = c(40,41,42)
rate = c(0.5,1,2)

#using coga calculate the convoluted distrobution
pcoga(1:200, shape=shape, rate=rate) %>% plot(1:200,.)
pcoga(1:200, shape=shape, rate=rate) %>% points(1:200/3,., col = 'orange')


#plot each independent function
for(i in 1:3){
  #compute the dist emulating prob dist fitting
  if(i == 1){
    data_new = rgamma(10000, shape[i], rate[i])
  }
  if(i!=1){
    data_new = append(data_new, rgamma(10000, shape[i], rate[i]))
  }
  if(i == 3){
    points(data_new, gamma_fit_spi(data_new, 'CDF'), col = 'blue')
  }
  
  cdfgam(1:200, vec2par(c(shape[i], 1/rate[i]), 'gam')) %>% lines(.)
  pcoga(1:200, shape=shape[i], rate=rate[i]) %>% points(., col = 'red')
}
