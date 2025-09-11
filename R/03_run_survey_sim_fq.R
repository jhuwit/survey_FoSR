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
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "00_data_gen_function_ff.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "create_settings.R"))
source(here::here("R", "create_settings2.R"))


force = FALSE
force_iter = FALSE
nsim = 200
B = 500
# ncores = parallelly::availableCores()
ncores = future::availableCores()

fit_types = c('weighted', 'unweighted')
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'weighted', 'unweighted')
options(survey.lonely.psu = "adjust")


ifold = get_fold()

outfile = here::here("results", "simulations", "survey_sim", paste0("fold_", sprintf("%03d", ifold), ".rds"))

if(!file.exists(outfile) || force) { # if the file doesn't exist or we want to force re-do it
  if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)
  partial_dir = here::here("results", "simulations", "survey_sim", paste0("fold_", sprintf("%03d", ifold), "_partials"))
  if (!dir.exists(partial_dir)) dir.create(partial_dir)

  temp = settings_new[ifold,]

  # check to see if old fold exists
  if(temp$inf_level == 0 || temp$inf_level == 15) {
    old_fold = settings %>%
      filter(scenario == temp$scenario,
             family == temp$family,
             len == temp$len,
             In == temp$In,
             snr_b == temp$snr_b,
             snr_eps == temp$snr_eps,
             strata_sigma == temp$strata_sigma,
             strata_scale == temp$strata_scale,
             inf_level == temp$inf_level / 3) %>% # old inf level is 0 or 15
      pull(fold)
    old_file =  here::here("results", "simulations", "survey_sim3", paste0("fold_", sprintf("%03d", old_fold), ".rds"))
    old_partials = here::here("results", "simulations", "survey_sim3", paste0("fold_", sprintf("%03d", old_fold), "_partials"))
    if(file.exists(old_file)) { # if it exists, just copy it over
      file.copy(from = old_file, to = outfile)
      message("copying over old file")
      quit(save = "no")
    }
    if(dir.exists(old_partials)) { # copy over old partials
      file.copy(
        from = list.files(old_partials, full.names = TRUE),
        to   = partial_dir,
        recursive = FALSE
      )
      message("copying over old partials")
    }
  }

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
  completed_iters = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$") %>%
    str_extract("\\d+") %>%
    as.integer()

  for (iter in 1:nsim) {
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


  unlink(partial_dir, recursive = TRUE)
}
