library(tidyverse)

data = read_csv('/home/zhoylman/drought-year-sensitivity/data/ncei_economic_summary/events-US-1980-2020.csv', skip = 1) %>%
  filter(Disaster == 'Drought')

test = data %>%
  mutate(year = substr(as.character(`Begin Date`), 1,4) %>% as.numeric(),
         Decade = .bincode(year, breaks = c(1979,1989,1999,2009,2019))) %>%
  group_by(Decade) %>%
  summarise(`Number of Billion Dollar Events` = length(`Total CPI-Adjusted Cost (Millions of Dollars)`)) %>%
  drop_na() %>%
  mutate(Decade = c('1980s', '1990s', '2000s', '2010s'))

plot = ggplot()+
  geom_bar(data = test, aes(x = Decade, y = `Number of Billion Dollar Events`), stat = 'identity')+
  theme_bw(base_size = 16)

plot

ggsave(plot, file = '/home/zhoylman/drought-year-sensitivity/figs/number_of_billion_dollar_events.png')
