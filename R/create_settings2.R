library(tidyverse)


settings1 = expand_grid(strata_scale = c(0, 0.125),
                        strata_sigma = c(0, 0.05),
                        snr_b = c(0.5, 1, 5),
                        snr_eps = c(0.5, 1, 5),
                        inf_level = c(0, 10, 15),
                        family = "gaussian",
                        len = c(50, 100),
                        In = c(100, 500),
                        scenario = 1)

settings1 =
  settings1 %>%
  filter(!(strata_scale == 0 & strata_sigma == 0 & snr_b %in% c(1, 5)))


settings2 =
  expand_grid(strata_scale = c(0, 0.125),
              strata_sigma = c(0, 0.05),
              snr_b = c(0.5, 1, 5),
              snr_eps = 1,
              inf_level = c(0, 10, 15),
              family = c("poisson", "binomial"),
              len = c(50, 100),
              In = c(100, 500),
              scenario = 1)

settings2 =
  settings2 %>%
  filter(!(strata_scale == 0 & strata_sigma == 0 & snr_b %in% c(1, 5)))

settings_new = settings1 %>%
  bind_rows(settings2) %>%
  mutate(fold = row_number())

rm(settings1, settings2)


