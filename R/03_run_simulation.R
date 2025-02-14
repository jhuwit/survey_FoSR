library(future)
library(furrr)
library(tidyverse)
force = TRUE
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "00_data_gen_function.R"))
source(here::here("R", "utils.R"))
sim_settings = read_rds(here::here("sim_data", "sim_settings.rds"))
ifold = get_fold()
fname = paste0("sim_fold_", sprintf("%03d", ifold), ".rds")

nsim = 500
parallel = TRUE
ncores = parallelly::availableCores() - 1
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'Preston', 'weighted', 'none')
fit_types = c('design', 'design', 'design', 'weighted', 'none')
options(survey.lonely.psu = "adjust")


if(!file.exists(here::here("results", "simulations", fname)) ||
   force) {
  tmp_settings = sim_settings %>% filter(fold == ifold)
  scenario = tmp_settings$scenario
  SNR_sigma = tmp_settings$sigma
  family = tmp_settings$family
  B = tmp_settings$num_boots
  I_n = tmp_settings$I_n
  sampled_strata = tmp_settings$sampled_strata

  L = 50
  grid = seq(0, 1, length = L)
  beta_true = matrix(NA, 2, L) # time-dependent beta-true trajectory

  if (scenario == 1) {
    beta_true[1, ] = -0.15 - 0.1 * sin(2 * grid * pi) - 0.1 * cos(2 * grid *
                                                                    pi)
    beta_true[2, ] = dnorm(grid, .6, .15) / 20
  } else if (scenario == 2) {
    beta_true[1, ] = 0.53 + 0.06 * sin(3 * grid * pi) - 0.03 * cos(6.5 * grid *
                                                                     pi)
    beta_true[2, ] = dnorm(grid, 0.2, .1) / 60 + dnorm(grid, 0.35, .1) /
      200 -
      dnorm(grid, 0.65, .06) / 250 + dnorm(grid, 1, .07) / 60
  }
  rownames(beta_true) <- c("Intercept", "x")

  I = 10e6 # superpopulation size
  set.seed(4574)
  X_des = cbind(1, rnorm(I, 0, 2)) # design matrix: Intercept + predicted variables
  true_fixef <- X_des %*% beta_true


  set.seed(4574)
  if (family == "gaussian") {
    Y_obs <- matrix(
      rnorm(
        n = nrow(true_fixef) * L,
        mean = as.vector(t(true_fixef)),
        sd = SNR_sigma
      ),
      nrow = nrow(true_fixef),
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "binomial") {
    p_true = plogis(true_fixef)
    Y_obs <- matrix(
      rbinom(
        n = nrow(p_true) * L,
        size = 1,
        prob = as.vector(t(p_true))
      ),
      nrow = nrow(p_true),
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "poisson") {
    lam_true <- exp(true_fixef)
    Y_obs <- matrix(
      rpois(n = nrow(lam_true) * L, lambda = as.vector(t(lam_true))),
      nrow = nrow(lam_true),
      ncol = L,
      byrow = TRUE
    )
  }

  sim_res <- list() ## store simulation results
  for (iter in 1:nsim) {
    set.seed(iter)
    x = try({
      data <- simulate_survey_fosr(
        X_des = X_des,
        Y_obs = Y_obs,
        beta_true = beta_true,
        true_fixef = true_fixef,
        L = 50,
        ## dimension of the functional domain
        I_n = I_n,
        # subjects in each psu-strata combination
        alpha = 2,
        # number PSU
        num_strata = 30,
        sampled_strata = sampled_strata
      )

      ## pre-process simulated data
      Y_mat <- data[, grepl('Y.', colnames(data))]
      dat.fit <- data.frame(Y = Y_mat, data[, !grepl('Y.', colnames(data))])
      model_formula = as.formula(paste0('Y~', 'X'))
      out_index <- grep(paste0("^", model_formula[2]), names(data))

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
          out_index = out_index,
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
          L = 50,
          out_index = out_index)

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

