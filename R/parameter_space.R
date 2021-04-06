library(tidyverse)

time_scale_id = 1
time_scale = list(30,60,90)

spi_comparison = readRDS(paste0('/home/zhoylman/temp', '/spi_comparision_moving_window_with_params_', time_scale[[time_scale_id]], '_days.RDS'))

full = bind_rows(spi_comparison) %>%
  filter( n_contemporary >= 25,
          n_historical >= 70,
          shape_contemporary > quantile(shape_contemporary, 0.01),
          shape_contemporary < quantile(shape_contemporary, 0.99),
          rate_contemporary > quantile(rate_contemporary, 0.01),
          rate_contemporary < quantile(rate_contemporary, 0.99))

random_index = runif(100000, 1, length(full$time)) %>% as.integer()

data_contemporary = full[random_index, ] %>%
  select(shape_contemporary, rate_contemporary) %>%
  rename(Shape = shape_contemporary, Rate = rate_contemporary) %>%
  mutate(ID = "Contemporary")

data_historical =full[random_index, ] %>%
  select(shape_historical, rate_historical) %>%
  rename(Shape = shape_historical, Rate = rate_historical) %>%
  mutate(ID = "Historical")

final = bind_rows(data_contemporary, data_historical)

ggplot(final, aes(x = Shape, y = Rate))+
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis", name = 'Count')+
  facet_grid(~ID)

ggplot(final, aes(x = Shape, y = Rate))+
  stat_density_2d(geom = "polygon", 
                  aes(fill = as.factor(..level..)),
                  bins = 10,
                  alpha = 0.5) +
  facet_grid(~ID)+
  scale_fill_viridis_d(guide = F)+
  theme_bw(base_size = 16)+
  labs(x = 'Shape', y = 'Rate')+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('Gamma Distrobution Parameter Space (30 Day Timescale)')
