library(tidyverse)
source(here::here("R_cp", "create_survey_settings_final.R"))
files = list.files(here::here("results", "mfpca_sim"), full.names = TRUE,
                   pattern = "*rds")


res = map_dfr(.x = files, .f = \(x) read_rds(x) %>%
                mutate(fold =
                       as.numeric(sub(".*fold\\_(.+).rds.*", "\\1", basename(x)))))

res =
  res %>%
  group_by(fold) %>%
  summarize(across(c(pve_b, pve_w, tv_bw, tv_win), mean), .groups = "drop")

res = res %>%
  left_join(settings, by = "fold")

res = res %>%
  mutate(type = case_when(
    strata_scale == 0 & strata_sigma == 0 ~ "No random effects",
    strata_scale > 0 & strata_sigma > 0 ~ "Strata scaling and noise",
    strata_scale > 0 & strata_sigma == 0 ~ "Strata scaling only",
    strata_scale == 0 & strata_sigma > 0 ~ "Strata noise only"
  ))

write_rds(res, here::here("results", "mfpca_sim.rds"))


######
