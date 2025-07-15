## joint CI visualization
# snr_eps snr_b    In   BRR
# <dbl> <dbl> <dbl> <dbl>
#   0.1   0.1   250 0.870
# 5     0.1   250 0.930
#  0.1   0.5   250 0.75


# visualize the data
library(future)
library(furrr)
library(tidyverse)
library(here)
library(devtools)
library(fastFMM)
library(haven)
library(dplyr)
library(survey)
library(progress)
library(lme4)
library(paletteer)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(tidyfun)
library(mvtnorm)
library(refund)
library(svrep)
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "00_data_gen_function_twostage.R"))
source(here::here("R", "utils.R"))
library(patchwork)


nsim = 250
B = 500
force = TRUE

ncores = parallelly::availableCores() - 1
fit_types = c('weighted')
boot_types = c('BRR')

options(survey.lonely.psu = "adjust")

# just use PSU RE = 0.5

temp = tibble(
  psu_re = 0.5,
  len = 50,
  snr_b = 0.5,
  snr_eps = 0.1,
  scen = 1,
  family = "gaussian",
  In = 250
)

strata_sigma = 0.05
psu_sigma = sqrt(strata_sigma ^ 2 * temp$psu_re)

lst = generate_superpopulation(
  scenario = temp$scen,
  family = temp$family,
  I = 10e6,
  L = temp$len,
  psu_sigma = psu_sigma,
  snr_b = temp$snr_b,
  snr_eps = temp$snr_eps
)
x = lst$stratum_assignments %>% table() %>% as_tibble() %>% arrange(n)
x
x %>% tail


x = lst$psu_assignments %>% table() %>% as_tibble() %>% arrange(n)
x
x %>% tail

rm(x)
rm(lst)

lst = generate_superpopulation(
  scenario = temp$scen,
  family = temp$family,
  I = 5 * 10e6,
  L = temp$len,
  psu_sigma = psu_sigma,
  snr_b = temp$snr_b,
  snr_eps = temp$snr_eps
)

x = lst$psu_assignments %>% table() %>% as_tibble() %>% arrange(n)
x
x %>% tail
rm(x)
rm(lst)

lst = generate_superpopulation(
  scenario = temp$scen,
  family = temp$family,
  I = 5 * 10e7,
  L = temp$len,
  psu_sigma = psu_sigma,
  snr_b = temp$snr_b,
  snr_eps = temp$snr_eps
)

x = lst$psu_assignments %>% table() %>% as_tibble() %>% arrange(n)
x
x %>% tail
