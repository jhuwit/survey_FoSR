library(tidyverse)

n = c(500, 1000, 5000, 10000)
sigma = c(0.5, 1, 5)
families = c("gaussian", "binomial", "poisson")
scenarios = c(1, 2)
num_boots = 500
L = c(50, 100, 500)
sim_settings = expand_grid(
  n = n,
  scenario = scenarios,
  sigma = sigma,
  num_boots = num_boots,
  family = families,
  L = L
) %>%
  mutate(fold = row_number())

# repetitive bc don't need noise for poisson or binomial
extra = sim_settings %>%
  filter(family != "gaussian") %>%
  filter(sigma != 0.5) %>%
  pull(fold)



sim_settings =
  sim_settings %>%
  filter(!(fold %in% extra)) %>%
  arrange(L, n) %>%
  mutate(fold = row_number())
