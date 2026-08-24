library(future)
library(furrr)
library(tidyverse)
library(here)
library(devtools)
library(fastFMM)
# library(haven)
library(survey)
library(progress)
library(lme4)
library(paletteer)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(tidyfun)
library(mvtnorm)
library(refund)
library(svrep)
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "00_data_gen_function_ff.R"))
source(here::here("R", "00_data_gen_function_fast.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "create_survey_settings.R"))


force = FALSE
force_iter = FALSE
B = 500
# ncores = parallelly::availableCores()
ncores = future::availableCores()

fit_types = c('weighted', 'unweighted')
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'weighted', 'unweighted')
options(survey.lonely.psu = "adjust")


ifold = get_fold()

temp = settings[ifold, ]

nsim = switch(
  as.character(temp$len),
  "50" = 200,
  "100" = 200,
  "1440" = 100
)

outfile = here::here("results",
                     "simulations",
                     "survey_sim",
                     paste0("fold_", sprintf("%03d", ifold), ".rds"))
if (file.exists(outfile)) {
  x = read_rds(outfile)
  if (nrow(x) >= nsim * length(boot_types) * 2) {
    message("Output file already has results for all iterations. Exiting.")
    quit(save = "no", status = 0)
  } else if (force) {
    message("Output file exists but force = TRUE, so re-running all iterations.")
  } else {
    message(
      "Output file exists but does not have results for all iterations. Continuing to run missing iterations."
    )
  }
}

if (!dir.exists(dirname(outfile)))
  dir.create(dirname(outfile), recursive = TRUE)
partial_dir = here::here("results",
                         "simulations",
                         "survey_sim",
                         paste0("fold_", sprintf("%03d", ifold), "_partials"))
if (!dir.exists(partial_dir))
  dir.create(partial_dir)
partial_files = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$", full.names = TRUE)

if (length(partial_files) >= nsim) {
  sim_res = map(partial_files, \(x) read_rds(x) %>% mutate(id = sub(
    ".*iter_(.+)\\.rds.*", "\\1", basename(x)
  ))) %>%
    keep(., is.data.frame) %>%
    list_rbind()

  write_rds(sim_res, file = outfile)

  unlink(partial_dir, recursive = TRUE)
  message("All iterations complete, writing and exiting.")
  quit(save = "no", status = 0)
}


if (temp$len != 1440 && !file.exists(here::here(partial_dir, "bh_data.rds"))) {
  lst = generate_superpopulation(
    family = temp$family,
    I = 10e6,
    L = temp$len,
    snr_b = temp$snr_b,
    snr_eps = temp$snr_eps,
    strata_sigma = temp$strata_sigma,
    strata_scale = temp$strata_scale,
    seed = 111
  )

  df = data.frame(X = lst$X_des[, 2], Y = lst$Y_obs)
  rm(lst); gc()
  bt_data = get_betatilde(family = temp$family,
                          type = "unweighted",
                          data = df)
  rm(df); gc()
  bh_data  = get_betahat(bt_data, L = temp$len)

  write_rds(bh_data, here::here(partial_dir, "bh_data.rds"))
}
