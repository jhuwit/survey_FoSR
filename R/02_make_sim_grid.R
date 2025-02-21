library(tidyverse)

I_n = c(30, 50, 100, 250, 500)
sigma = c(0.5, 1, 1.5)
sampled_strata = c(10, 15, 30)
num_boots = c(100, 500, 1000)
families = c("gaussian", "binomial", "poisson")
scenarios = c(1, 2)


sim_settings = expand_grid(
  I_n = I_n,
  scenario = scenarios,
  sigma = sigma,
  sampled_strata = sampled_strata,
  num_boots = num_boots,
  family = families
) %>%
  mutate(fold = row_number())

# we don't need noise for poisson or binomial

extra = sim_settings %>%
  filter(family != "gaussian") %>%
  group_by(I_n, scenario, sampled_strata, family, num_boots) %>%
  slice(2:3) %>%
  pull(fold)

sim_settings =
  sim_settings %>%
  filter(!(fold %in% extra)) %>%
  mutate(fold = row_number())

write_rds(sim_settings, here::here("sim_data", "sim_settings.rds"))

sim_settings = expand_grid(
  I_n = I_n,
  scenario = scenarios,
  sigma = sigma,
  sampled_strata = sampled_strata,
  num_boots = num_boots,
  family = "gaussian"
) %>%
  mutate(fold = row_number())

write_rds(sim_settings, here::here("sim_data", "sim_settings_inf.rds"))
