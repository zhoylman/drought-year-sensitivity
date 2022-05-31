library(tidyverse)

data = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/drought-metrics/spi-data-wide-10s.csv')

vwc_data = read_csv('/home/zhoylman/soil-moisture-validation-data/processed/merged-soil-moisture/soil-moisture-data-long.csv')

sites = unique(vwc$site_id)

drought = data %>%
  filter(site_id == sites[1]) %>%
  pivot_longer(cols = -c(site_id, time))

vwc = vwc_data %>%
  filter(site_id == sites[1]) 

drought_correlation = function(drought, vwc){
  # remove frozen timestamps
  test = vwc %>%
    pivot_wider(names_from = name, values_from = value)
    
}