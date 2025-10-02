library(future)
library(furrr)
library(tidyverse)
library(here)
library(devtools)
library(fastFMM)
library(haven)
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
source(here::here("R_cp", "01_sim_functions.R"))
source(here::here("R_cp", "00_data_gen_function_ff.R"))
source(here::here("R", "utils.R"))


settings1 = expand_grid(strata_scale = c(0, 0.125),
                        strata_sigma = c(0, 0.05),
                        snr_b = c(0.5, 1, 5),
                        snr_eps = c(0.5, 1, 5),
                        inf_level = c(0, 10, 15),
                        family = "gaussian",
                        len = 1440,
                        In = 100,
                        scenario = 1)

settings1 =
  settings1 %>%
  filter(!(strata_scale == 0 & strata_sigma == 0 & snr_b %in% c(1, 5))) %>% # snr_b is not applicable when there are no random effects
  filter(!(strata_scale == 0 & strata_sigma > 0)) %>%
  filter(!(strata_scale > 0 & strata_sigma == 0))


settings2 =
  expand_grid(strata_scale = c(0, 0.125),
              strata_sigma = c(0, 0.05),
              snr_b = c(0.5, 1, 5),
              snr_eps = 1,
              inf_level = c(0, 10, 15),
              family = c("poisson", "binomial"),
              len = 1440,
              In = 100,
              scenario = 1)

settings2 =
  settings2 %>%
  filter(!(strata_scale == 0 & strata_sigma == 0 & snr_b %in% c(1, 5))) %>%
  filter(!(strata_scale == 0 & strata_sigma > 0)) %>%
  filter(!(strata_scale > 0 & strata_sigma == 0))

settings_new = settings1 %>%
  bind_rows(settings2) %>%
  mutate(fold = row_number())

rm(settings1, settings2)



force = FALSE
force_iter = FALSE
nsim = 100
B = 500
# ncores = parallelly::availableCores()
ncores = future::availableCores()

fit_types = c('weighted', 'unweighted')
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'weighted', 'unweighted')
options(survey.lonely.psu = "adjust")


ifold = get_fold()

outfile = here::here("results", "simulations", "survey_sim_1440", paste0("fold_", sprintf("%03d", ifold), ".rds"))

if(!file.exists(outfile) || force) { # if the file doesn't exist or we want to force re-do it
  if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)
  partial_dir = here::here("results", "simulations", "survey_sim_1440", paste0("fold_", sprintf("%03d", ifold), "_partials"))
  if (!dir.exists(partial_dir)) dir.create(partial_dir)

  temp = settings_new[ifold,]


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

  # check for completed iters
  # completed_iters = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$") %>%
  #   str_extract("\\d+") %>%
  #   as.integer()

  for (iter in 1:nsim) {
    completed_iters = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$") %>%
      str_extract("\\d+") %>%
      as.integer()

    if (iter %in% completed_iters && !force_iter) {
      message("Skipping completed iteration: ", iter)
      next
    }

    set.seed(iter)
    x = try({
      data = sample_from_population_wor(
        X_des = lst$X_des,
        Y_obs = lst$Y_obs,
        L = temp$len,
        I_n = temp$In,
        num_strata = 30,
        stratum_assignments = lst$stratum_assignments,
        psu_assignments = lst$psu_assignments,
        dirichlet_probs = lst$dirichlet_probs,
        seed = iter,
        inf_level = temp$inf_level,
        compression = 2,
        family = temp$family
      )

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
        .f = run_boots_xfast,
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
        .f = get_cis_bam,
        L = temp$len,
        nknots_min = 35,
        nknots_min_cov = 35,
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

    })
    rm(x)
    print(iter)
  }

  plan(sequential)

  partial_files = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$", full.names = TRUE)

  sim_res = map(partial_files, \(x) read_rds(x) %>% mutate(id = sub(".*iter_(.+)\\.rds.*", "\\1", basename(x)))) %>%
    keep(., is.data.frame) %>%
    list_rbind()

  write_rds(sim_res, file = outfile)


  unlink(partial_dir, recursive = TRUE)
}
