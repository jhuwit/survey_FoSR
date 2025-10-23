library(tidyverse)
# getremote survey_FoSR/results/mfpca_sim3/*.rds results/mfpca_sim3/
source(here::here("R", "create_survey_settings.R"))
files = list.files(here::here("results", "mfpca_sim3"), full.names = TRUE)


res = map_dfr(.x = files, .f = \(x) read_rds(x) %>%
                mutate(fold =
                         as.numeric(sub(".*fold\\_(.+).rds.*", "\\1", basename(x))),
                       iter = row_number()))


res = res %>%
  left_join(settings_new, by = "fold")

res = res %>%
  mutate(type = case_when(
    strata_scale == 0 & strata_sigma == 0 ~ "No random effects",
    strata_scale > 0 & strata_sigma > 0 ~ "Strata scaling and noise",
    strata_scale > 0 & strata_sigma == 0 ~ "Strata scaling only",
    strata_scale == 0 & strata_sigma > 0 ~ "Strata noise only"
  ))

write_rds(res, here::here("results", "mfpca_sim.rds"))
