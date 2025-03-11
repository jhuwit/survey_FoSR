library(future)
library(furrr)
library(tidyverse)
force = TRUE
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "00_data_gen_function_corr2.R"))
source(here::here("R", "utils.R"))
sim_settings = read_rds(here::here("sim_data", "sim_settings.rds"))
ifold = get_fold()
fname = paste0("sim_fold_", sprintf("%03d", ifold), ".rds")

nsim = 500
parallel = TRUE
ncores = parallelly::availableCores() - 2
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'Preston', 'weighted', 'none')
fit_types = c('design', 'design', 'design', 'weighted', 'none')
options(survey.lonely.psu = "adjust")

# ifold = 1
if(!file.exists(here::here("results", "simulations", fname)) ||
   force) {


  tmp_settings = sim_settings %>% filter(fold == ifold) # set settings for this simulation

  # I_n scenario sigma sampled_strata num_boots family    fold
  # <dbl>    <dbl> <dbl>          <dbl>     <dbl> <chr>    <int>
  #   1    30        1   0.5             10       100 gaussian     1
  scenario = tmp_settings$scenario
  SNR_sigma = tmp_settings$sigma
  family = tmp_settings$family
  B = tmp_settings$num_boots
  I_n = tmp_settings$I_n
  sampled_strata = tmp_settings$sampled_strata

  lst = generate_superpopulation(scenario = scenario,
                                 family = family)


  sim_res <- list() ## store simulation results
  for (iter in 1:nsim) {
    set.seed(iter)
    x = try({
      data <- simulate_survey_fosr_corr(
        X_des = lst$X_des,
        Y_obs = lst$Y_obs,
        ## dimension of the functional domain
        I_n = I_n,
        # subjects in each psu-strata combination
        # number PSU
        num_strata = 30,
        sampled_strata = sampled_strata,
        stratum_assignments = lst$stratum_assignments,
        psu_assignments = lst$psu_assignments,
        dirichlet_probs = lst$dirichlet_probs
      )

      beta_true = lst$beta_true
      betaTilde = map(.x = fit_types, .f = get_betatilde, data = data, family = family)

      betaHat = map(.x = betaTilde, get_betahat)

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
          L = 50,
          .options = furrr_options(seed = TRUE)
        )

        cis = future_map2(
          .x = res,
          .y = betaHat,
          .f = get_cis,
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
          seed = 2025,
          L = 50)

        cis = map2(
          .x = res,
          .y = betaHat,
          .f = get_cis,
          smooth_for_ci = TRUE)
      }

      stats = map2(
        .x = cis,
        .y = boot_types,
        .f = get_coverage_stats,
        beta_true = beta_true
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

  if (!dir.exists(here::here("results", "simulations"))) {
    dir.create(here::here("results", "simulations"), recursive = TRUE)
  }

  write_rds(result, here::here("results", "simulations", fname))
}

