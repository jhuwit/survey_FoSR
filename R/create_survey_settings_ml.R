### survey settings final

library(tidyr)
library(dplyr)


settings1 = tidyr::expand_grid(strata_scale = c(0, 0.125),
                               strata_sigma = c(0, 0.05),
                               snr_b = c(0.5, 1, 2),
                               snr_eps = c(0.5, 1, 2),
                               snr_u = c(0.5, 1, 2),
                               inf_level = c(0, 10),
                               family = "gaussian",
                               len = c(50, 100, 200),
                               In = c(100, 500))

settings1 =
  settings1 %>%
  filter((strata_scale == 0 & strata_sigma == 0) | (strata_scale > 0 & strata_sigma > 0)) |>
  mutate(snr_b = if_else(strata_scale == 0 & strata_sigma == 0, NA_real_, snr_b)) %>%
  distinct()

settings2 =
  expand_grid(strata_scale = c(0, 0.125),
              strata_sigma = c(0, 0.05),
              snr_b = c(0.5, 1, 2),
              snr_eps = NA_real_,
              inf_level = c(0, 10),
              family = c("poisson", "binomial"),
              len = c(50, 100, 200),
              In = c(100, 500))

settings2 =
  settings2 %>%
  filter((strata_scale == 0 & strata_sigma == 0) | (strata_scale > 0 & strata_sigma > 0)) |>
  mutate(snr_b = if_else(strata_scale == 0 & strata_sigma == 0, NA_real_, snr_b)) %>%
  distinct()

settings = settings1 %>%
  bind_rows(settings2) %>%
  mutate(fold = row_number())

