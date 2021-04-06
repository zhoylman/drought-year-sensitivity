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
  