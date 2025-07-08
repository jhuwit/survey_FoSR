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
source(here::here("R", "00_data_gen_function_twostage.R"))
source(here::here("R", "utils.R"))

ifold = get_fold()
fname = paste0("sim_fold_", sprintf("%03d", ifold), ".rds")

nsim = 250
B = 500
parallel = TRUE
force = TRUE

ncores = parallelly::availableCores() - 1
fit_types = c('weighted', 'weighted', 'weighted', 'weighted', 'unweighted', 'unweighted')
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'weighted', 'unweighted', 'weighted', 'unweighted')

options(survey.lonely.psu = "adjust")


settings1 = expand_grid(snr_b = c(0.1, 0.5, 1, 2.5, 5, 10),
                       snr_eps = c(0.1, 0.5, 1, 2.5, 5),
                       psu_re = c(.5, .25),
                       In = c(100, 250, 500)) %>%
  mutate(len = 50,
         family = "gaussian",
         scen = 1)

settings2 = expand_grid(snr_b = c(0.1, 0.5, 1, 2.5, 5, 10),
                        # snr_eps = c(0.1, 0.5, 1, 2.5, 5),
                        psu_re = 0.5,
                        family = c("poisson", "binomial"),
                        In = c(100, 250, 500)) %>%
  mutate(len = 50,
         scen = 1)

settings3 = expand_grid(snr_b = c(0.1, 0.5, 1, 2.5, 5, 10),
                        snr_eps = c(0.1, 0.5, 1, 2.5, 5),
                        psu_re = c(.5),
                        In = c(100, 250, 500)) %>%
  mutate(len = 100,
         family = "gaussian",
         scen = 1)

settings4 = expand_grid(snr_b = c(0.1, 0.5, 1, 2.5, 5, 10),
                        # snr_eps = c(0.1, 0.5, 1, 2.5, 5),
                        psu_re = 0.5,
                        family = c("poisson", "binomial"),
                        In = c(100, 250, 500))  %>%
  mutate(len = 100,
         scen = 1)

settings5 = expand_grid(snr_b = c(0.1, 0.5, 1, 2.5, 5, 10),
                        snr_eps = c(0.1, 0.5, 1, 2.5, 5),
                        psu_re = c(.5),
                        In = c(100, 250, 500)) %>%
  mutate(len = 50,
         family = "gaussian",
         scen = 2)

settings6 = expand_grid(snr_b = c(0.1, 0.5, 1, 2.5, 5, 10),
                        # snr_eps = c(0.1, 0.5, 1, 2.5, 5),
                        psu_re = 0.5,
                        family = c("poisson", "binomial"),
                        In = c(100, 250, 500)) %>%
  mutate(len = 50,
         scen = 2)

settings7 = expand_grid(snr_b = c(0.1, 0.5, 1, 2.5, 5, 10),
                        snr_eps = c(0.1, 0.5, 1, 2.5, 5),
                        psu_re = c(.5),
                        In = c(100, 250, 500)) %>%
  mutate(len = 100,
         family = "gaussian",
         scen = 2)

settings8 = expand_grid(snr_b = c(0.1, 0.5, 1, 2.5, 5, 10),
                        # snr_eps = c(0.1, 0.5, 1, 2.5, 5),
                        psu_re = 0.5,
                        family = c("poisson", "binomial"),
                        In = c(100, 250, 500))  %>%
  mutate(len = 100,
         scen = 2)

settings = bind_rows(settings1, settings2,
                     settings3, settings4,
                     settings5, settings6,
                     settings7, settings8)

# just use PSU RE = 0.5

if(!file.exists(here::here("results", "simulations", "two_stage_sim", fname)) || force) {


  temp = settings[ifold,] # set settings for this simulation

  strata_sigma = 0.05
  psu_sigma = sqrt(strata_sigma ^ 2 * temp$psu_re)

  lst = generate_superpopulation(scenario = temp$scen,
                                 family = temp$family,
                                 I = 10e6,
                                 psu_sigma = psu_sigma,
                                 snr_b = temp$snr_b,
                                 snr_eps = temp$snr_eps)



  sim_res <- list() ## store simulation results
  for (iter in 1:nsim) {
    set.seed(iter)
    x = try({
      data <- sample_from_population(
        X_des = lst$X_des,
        Y_obs = lst$Y_obs,
        ## dimension of the functional domain
        I_n = temp$In,
        num_strata = 30,
        stratum_assignments = lst$stratum_assignments,
        psu_assignments = lst$psu_assignments,
        dirichlet_probs = lst$dirichlet_probs
      )

      beta_true = lst$beta_true
      betaTilde = map(.x = fit_types, .f = get_betatilde, data = data, family = family)

      betaHat = map(.x = betaTilde, get_betahat, L = temp$len)

      if (parallel) {
        plan(multisession, workers = ncores)
        res = future_map2(
          .x = boot_types,
          .y = betaHat,
          .f = run_boots,
          data = data,
          family = family,
          num_boots = B,
          set_seed = TRUE,
          seed = 2025,
          L = temp$len,
          samp_stages = c("PPSWR", "Poisson"),
          .options = furrr_options(seed = TRUE)
        )

        cis = future_map2(
          .x = res,
          .y = betaHat,
          .f = get_cis,
          L = temp$len,
          smooth_for_ci = TRUE,
          .options = furrr_options(seed = TRUE)
        )
        plan(sequential)
      } else{
        res = map2(
          .x = boot_types,
          .y = betaHat,
          .f = run_boots,
          data = data,
          family = family,
          num_boots = B,
          set_seed = TRUE,
          samp_stages = c("PPSWR", "Poisson"),
          seed = 2025,
          L = temp$len)

        cis = map2(
          .x = res,
          .y = betaHat,
          .f = get_cis,
          L = temp$len,
          smooth_for_ci = TRUE)
      }

      stats = map2(
        .x = cis,
        .y = boot_types,
        .f = get_coverage_stats,
        beta_true = beta_true,
        L = temp$len
      ) %>%
        list_rbind() %>%
        mutate(n = nrow(data), n_boot = B)

      sim_res[[iter]] <- stats
    })
    rm(x)
    print(iter)
  }

  result =
    sim_res %>%
    keep(., is.data.frame) %>%
    list_rbind(names_to = "id")

  if (!dir.exists(here::here("results", "simulations", "two_stage_sim"))) {
    dir.create(here::here("results", "simulations", "two_stage_sim"), recursive = TRUE)
  }

  write_rds(result, here::here("results", "simulations", "two_stage_sim", fname))
}

