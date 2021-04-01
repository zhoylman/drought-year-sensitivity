# generate a synthetic distrobution of known parameters
# to evalaute how much data it takes to apporximate the
# known distrobtion
library(tidyverse)
library(lmomco)
library(cowplot)

#define functions to compute probabilistic CDF and probabilistic parameters
gamma_cdf = function(x) {
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
      return(fit.cdf) 
    },
    #else return NA
    error=function(cond) {
      return(NA)
    })
}
gamma_params = function(x) {
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
      #compute standard normal equivelant
      standard_norm = qnorm(fit.cdf, mean = 0, sd = 1)
      return(fit.gam) 
    },
    #else return NA
    error=function(cond) {
      return(NA)
    })
}

#######################################################################
####################### STATIONARY DISTROBUTION #######################
#######################################################################

#define proabaility distrobution paramters
shape = 40 #alpha
rate = 0.7 #beta = 1/rate

#set seed for reproducability
set.seed(99)

#generate synthetic data (synthetic precipitation distrobution based values from Wu et al., 2005)
data = rgamma(100, shape, rate)

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
    probabilistic_vals = gamma_cdf(temp_data)
    #compute the Mean Absolute Error Prob - True
    mae = mean(abs(probabilistic_vals - true_vals))
    
    #SPI comparison
    true_spi_vals = qnorm(true_vals)
    #comupute probabilistic vals from limited climatology
    probabilistic_spi_vals = qnorm(probabilistic_vals)
    #compute MAE of the spi vals
    mae_spi = mean(abs(probabilistic_spi_vals - true_spi_vals))
    
    #compute parameters of probability dist given the limited sample
    params = gamma_params(temp_data)
    
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

which(summary_mae$upper < 0.05) %>% first()
which(summary_spi_mae$upper < 0.2) %>% first()


#plot the results
plot_mae = ggplot(data = summary_mae, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = 'Number of Observations in "Climatology"', y = 'CDF Mean Absolute Error\n(Limited Climatology - True Value)')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept=0.05, linetype="dashed", color = "red")+
  geom_vline(xintercept = which(summary_mae$upper < 0.05) %>% first(), linetype="dashed", color = "red")+
  geom_text(data = NULL, aes(x =65,
                             y = 0.08, label = 'CDF MAE for 75% of Simulations:\n< 0.05 at 33 years'), hjust = 0.5)+
  geom_text(data = NULL, aes(x =65, y = 0, label = 'Median CDF MAE = 0.0323 at 33 years'), hjust = 0.5)+
  geom_segment(data = NULL, aes(x = 34, y = 0.028, xend = 37.5, yend = 0.005))

plot_spi_mae = ggplot(data = summary_spi_mae, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = 'Number of Observations in "Climatology"', y = 'SPI Mean Absolute Error\n(Limited Climatology - True Value)')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_hline(yintercept=0.2, linetype="dashed", color = "red")+
  geom_vline(xintercept = which(summary_spi_mae$upper < 0.2) %>% first(), linetype="dashed", color = "red")+
  geom_text(data = NULL, aes(x =65,
                             y = 0.3, label = 'SPI MAE for 75% of Simulations:\n< 0.2 at 31 years'), hjust = 0.5)+
  geom_text(data = NULL, aes(x =65,
                             y = 0, label = 'Median SPI MAE = 0.139 at 31 years'), hjust = 0.5)+
  geom_segment(data = NULL, aes(x = 32, y = 0.12, xend = 39, yend = 0.02))

plot_spi_mae

plot_rate = ggplot(data = summary_rate, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = NULL, y = 'Rate Parameter')+
  #ggtitle('Monte Carlo Simulation (1000 simulations)')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_text(data = NULL, aes(x = 65, y = 1.5, label = 'True Rate Parameter = 0.7'), size = 5)+
  geom_hline(yintercept=0.7, linetype="dashed", color = "red")

plot_shape = ggplot(data = summary_shape, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = NULL, y = 'Shape Parameter')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_text(data = NULL, aes(x = 65, y = 80, label = 'True Shape Parameter = 40'), size = 5)+
  geom_hline(yintercept=40, linetype="dashed", color = "red")

#combine to single plot
plot_final = cowplot::plot_grid(plot_rate, plot_shape,plot_mae, plot_spi_mae,
                                nrow = 2, ncol = 2, align = 'hv')

title = ggdraw() + 
  draw_label(
    "Monte Carlo Simulation (1000 iterations)",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 20
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 300)
  )

final = plot_grid(
  title, plot_final,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

#save
ggsave(final, file = '/home/zhoylman/drought-year-sensitivity/figs/year_sensitivity_monte_carlo.png', width = 13, height = 9, units = 'in')
#fin stationary dist





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
shape_non_stationary = seq(40,49, length.out = 100) # alpha
rate_non_stationary = seq(0.6,3, length.out = 100) # rate

input_matrix = data.frame(x = n_, shape = shape_non_stationary, rate = rate_non_stationary)

n_simulation = 100

#set up export data frames
export_non_stationary_df_mae = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_non_stationary_df_shape = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_non_stationary_df_rate = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_non_stationary_rate_shape_full = data.frame(matrix(nrow = n_simulation, ncol = 2))

#non-stationary data
for(i_n_simulations in 1:n_simulation){
  print(i_n_simulations)
  #randomly generate the distrobtuion each simulation
  data_non_stationary = input_matrix %>%
    mutate(data =  rgamma(x, shape, rate)) %>%
    select(data)%>%
    mutate(time = seq(1:100))
  export_non_stationary_rate_shape_full[i_n_simulations,1] = gamma_params(data_non_stationary$data)$para[1]
  export_non_stationary_rate_shape_full[i_n_simulations,2] = gamma_params(data_non_stationary$data)$para[2]
  
  for(i_n_samples in n_samples){
    #pull data from the end of teh distrobution backwards
    #simulates what we would do in the case of a moving window analysis
    #to test parameter stability 
    temp_data = data_non_stationary[length(data_non_stationary$data):(length(data_non_stationary$data)-(i_n_samples-1)),] %>%
      mutate(contemporary_vals = gamma_cdf(data))
    #here we will simulate the assumptions of Wu and Guttman
    historical_values = data_non_stationary %>%
      mutate(historical_values = gamma_cdf(data))
    #join data
    joined = left_join(historical_values, temp_data, by = 'time') %>%
      drop_na() %>%
      mutate(diff = abs(contemporary_vals - historical_values))
    
    mae = mean(joined$diff)
    #compute parameters of probability dist given the limited sample
    params = gamma_params(temp_data$contemporary_vals)
    #populate the export dataframes
    export_non_stationary_df_mae[i_n_samples, i_n_simulations] = mae
    if(is.na(params) == T){
      export_non_stationary_df_shape[i_n_samples, i_n_simulations] = NA
      export_non_stationary_df_rate[i_n_samples, i_n_simulations] = NA
    }
    if(is.na(params) == F){
      export_non_stationary_df_shape[i_n_samples, i_n_simulations] = params$para[1]
      export_non_stationary_df_rate[i_n_samples, i_n_simulations] = params$para[2]
    }
  }
}

summary_non_stationary_mae = summarize_fun(export_non_stationary_df_mae)
summary_non_stationary_rate = summarize_fun(export_non_stationary_df_rate)
summary_non_stationary_shape = summarize_fun(export_non_stationary_df_shape)

#plot the results
plot_mae = ggplot(data = summary_non_stationary_mae, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = 'Number of Observations in "Climatology"', y = 'Mean Absolute Error\n(Contempary - Historical)')+
  theme(plot.title = element_text(hjust = 0.5))

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
plot_final = cowplot::plot_grid(plot_rate, plot_shape, plot_mae, nrow = 3)
#save
ggsave(plot_final, file = '/home/zhoylman/drought-year-sensitivity/figs/year_sensitivity_monte_carlo_non_stationary.png', width = 6, height = 8, units = 'in')
#fin stationary dist
