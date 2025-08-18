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
source(here::here("R", "00_data_gen_function_ff.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "create_settings.R"))

force = FALSE
force_iter = FALSE
nsim = 200
B = 500
ncores = parallelly::availableCores() - 1
fit_types = c('weighted', 'unweighted')
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'weighted', 'unweighted')
options(survey.lonely.psu = "adjust")


ifold = get_fold()

temp = settings %>%
  filter(fold == ifold)

outfile = here::here("results", "simulations", "survey_sim2", paste0("fold_", sprintf("%03d", ifold), ".rds"))
outfile2 = here::here("results", "simulations", "survey_sim2", paste0("fold_", sprintf("%03d", ifold), "_res.rds"))
# create partial dir to store ea. iteration

if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)

partial_dir = here::here("results", "simulations", "survey_sim2", paste0("fold_", sprintf("%03d", ifold), "_partials"))
if (!dir.exists(partial_dir)) dir.create(partial_dir)

partial_dir2 = here::here("results", "simulations", "survey_sim2", paste0("fold_", sprintf("%03d", ifold), "_partials2"))
if (!dir.exists(partial_dir2)) dir.create(partial_dir2)

completed_iters = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$") %>%
  str_extract("\\d+") %>%
  as.integer()

if(!file.exists(outfile) || force) {
  temp = settings[ifold, ] # set settings for this simulation

  lst = generate_superpopulation(
    scenario = temp$scenario,
    family = temp$family,
    I = 10e6,
    L = temp$len,
    snr_b = temp$snr_b,
    snr_eps = temp$snr_eps,
    strata_sigma = temp$strata_sigma,
    strata_scale = temp$strata_scale,
    seed = 111
  )


  plan(multisession, workers = ncores)

  for (iter in 1:nsim) {
    if (iter %in% completed_iters && !force_iter) {
      message("Skipping completed iteration: ", iter)
      next
    }

    set.seed(iter)
    x = try({
      data <- sample_from_population_wor(
        X_des = lst$X_des,
        Y_obs = lst$Y_obs,
        L = temp$len,
        I_n = temp$In,
        num_strata = 30,
        stratum_assignments = lst$stratum_assignments,
        psu_assignments = lst$psu_assignments,
        dirichlet_probs = lst$dirichlet_probs,
        seed = iter,
        inf_level = temp$inf_level * 4,
        compression = 2,
      )
      print(summary(data$weight))

      beta_true = lst$beta_true
      betaTilde = map(
        .x = fit_types,
        .f = get_betatilde,
        data = data,
        family = temp$family
      )

      betaHat = map(.x = betaTilde, get_betahat, L = temp$len)


      res = future_map(
        .x = boot_types,
        .f = run_boots_fast,
        betaHat = betaHat[[1]], # don't actually use betahat in calculation so we can just use the weighted one
        data = data,
        family = temp$family,
        num_boots = B,
        seed = 2025,
        L = temp$len,
        samp_stages = c("PPSWOR", "Poisson"),
        .options = furrr_options(seed = TRUE)
      )

      fit_list = list(betaHat[[1]], betaHat[[1]], betaHat[[1]], betaHat[[2]])

      cis = future_map2(
        .x = res,
        .y = fit_list,
        .f = get_cis,
        L = temp$len,
        smooth_for_ci = TRUE,
        .options = furrr_options(seed = TRUE)
      )

      stats = future_map2(
        .x = cis,
        .y = boot_types,
        .f = get_coverage_stats,
        beta_true = beta_true,
        L = temp$len
      ) %>%
        list_rbind() %>%
        mutate(n = nrow(data), n_boot = B)


      write_rds(stats, file = file.path(partial_dir, sprintf("iter_%03d.rds", iter)))
      write_rds(res, file = file.path(partial_dir2, sprintf("iter_%03d.rds", iter)))

    })
    rm(x)
    print(iter)
  }

  plan(sequential)

  partial_files = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$", full.names = TRUE)

  sim_res = map(partial_files, read_rds) %>%
    keep(., is.data.frame) %>%
    list_rbind(names_to = "id")

  write_rds(sim_res, file = outfile)

  partial_files2 = list.files(partial_dir2, pattern = "^iter_\\d+\\.rds$", full.names = TRUE)

  sim_res = map(partial_files, read_rds)

  write_rds(sim_res, file = outfile2)

  unlink(partial_dir, recursive = TRUE)
  unlink(partial_dir2, recursive = TRUE)


}


