# generate a synthetic distrobution of known parameters
# to evalaute how much data it takes to apporximate the
# known distrobtion
library(tidyverse)

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
      #compute standard normal equivelant
      standard_norm = qnorm(fit.cdf, mean = 0, sd = 1)
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

#define proabaility distrobution paramters
shape = 40 #alpha
rate = 0.7 #beat = 1/rate

#set seed for reproducability
set.seed(99)

#generate synthetic data (synthetic precipitation distrobution based on Wu et al., 2005)
data = rgamma(100, shape, rate)

#define monte carlo information
n_samples = seq(1,100,1)
n_simulation = 1000

#set up export data frames
export_df_mae = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_df_shape = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))
export_df_rate = data.frame(matrix(nrow = length(n_samples), ncol = n_simulation))

#run monte carlo
for(i_n_simulations in 1:n_simulation){
  for(i_n_samples in n_samples){
    #randomly sample from the synthetic data above
    temp_data = sample(data, n_samples[i_n_samples])
    #compute the emperical CDF
    emperical_model = ecdf(temp_data) 
    emperical_vals = emperical_model(temp_data)
    #compute probabilistic CDF
    probabilistic_vals = gamma_cdf(temp_data)
    #compute the Mean Absolute Error Prob - Emperical
    mae = mean(abs(probabilistic_vals - emperical_vals))
    #compute parameters of probability dist given the limited sample
    params = gamma_params(temp_data)
    #populate the export dataframes
    export_df_mae[i_n_samples, i_n_simulations] = mae
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
    pivot_longer(cols = -c(n_obs)) %>%
    dplyr::select(-name) %>%
    group_by(n_obs) %>%
    summarise(median = median(value, na.rm = T),
              upper =  quantile(value, 0.75, na.rm = T),
              lower = quantile(value, 0.25, na.rm = T))
  return(export)
}

#summarize monte carlo results
summary_mae = summarize_fun(export_df_mae)
summary_rate = summarize_fun(export_df_rate)
summary_shape = summarize_fun(export_df_shape)

#plot the results
plot_mae = ggplot(data = summary_mae, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = 'Number of Observations in "Climatology"', y = 'Mean Absolute Error\n(Probabilistic - Emperical)')+
  theme(plot.title = element_text(hjust = 0.5))

plot_rate = ggplot(data = summary_rate, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = NULL, y = 'Rate Parameter')+
  ggtitle('Monte Carlo Simulation (1000 simulations)')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_text(data = NULL, aes(x = 80, y = 5, label = 'True Rate Parameter = 0.7'))+
  geom_hline(yintercept=0.7, linetype="dashed", color = "red")

plot_shape = ggplot(data = summary_shape, aes(x = n_obs, y = median, ymax = upper, ymin = lower))+
  geom_ribbon(fill = 'grey70')+
  geom_line()+
  theme_bw(base_size = 16)+
  labs(x = NULL, y = 'Shape Parameter')+
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_text(data = NULL, aes(x = 80, y = 280, label = 'True Shape Parameter = 40'))+
  geom_hline(yintercept=40, linetype="dashed", color = "red")

#combine to single plot
plot_final = cowplot::plot_grid(plot_rate, plot_shape, plot_mae, nrow = 3)
#save
ggsave(plot_final, file = '/home/zhoylman/drought-year-sensitivity/figs/year_sensitivity_monte_carlo.png', width = 6, height = 8, units = 'in')
#fin