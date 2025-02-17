rm(list = ls())
library(future)
library(furrr)
library(tidyverse)
force = FALSE
source(here::here("R", "01_sim_functions_ml.R"))
source(here::here("R", "00_data_gen_function_ml.R"))
source(here::here("R", "utils.R"))
sim_settings = read_rds(here::here("sim_data", "sim_settings.rds"))
ifold = get_fold()
fname = paste0("sim_fold_", sprintf("%03d", ifold), ".rds")

nsim = 500
parallel = FALSE
ncores = parallelly::availableCores() - 1
boot_types = c('BRR', 'bootstrap', 'JKn', 'weighted', 'none')
fit_types = c('svy2lme', 'svy2lme', 'svy2lme', 'lmer_wtd', 'lmer')
options(survey.lonely.psu = "adjust")


if(!file.exists(here::here("results", "simulations", fname)) ||
   force) {

  scenario <- 1 ## indicate the scenario to use
  L = 50
  I = 10e5
  J = 5
  SNR_B = 0.5
  SNR_sigma = 1
  family = "gaussian"
  set.seed(123)
  J_subj <- pmax(rpois(I, J), 1) # number of visits per subject
  n <- sum(J_subj) # total number of visits and subjects
  set.seed(123)
  # design matrix with random effects
  X_des = cbind(1, rnorm(n, 0, 2)) # design matrix (n x 2)

  ## specify true fixed effects functions
  grid <- seq(0, 1, length = L)
  beta_true <- matrix(NA, 2, L) # time-dependent beta-true trajectory
  if(scenario == 1){
    beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
    beta_true[2,] = dnorm(grid, .6, .15)/20
  }else if(scenario == 2){
    beta_true[1,] <- 0.53 + 0.06*sin(3*grid*pi) - 0.03*cos(6.5*grid*pi)
    beta_true[2,] <- dnorm(grid, 0.2, .1)/60 + dnorm(grid, 0.35, .1)/200 -
      dnorm(grid, 0.65, .06)/250 + dnorm(grid, 1, .07)/60
  }
  rownames(beta_true) <- c("Intercept", "x")

  # random effects
  psi_true <- matrix(NA, 2, L)
  psi_true[1,] <- (1.5 - sin(2*grid*pi) - cos(2*grid*pi) )
  psi_true[1,] <- psi_true[1,] / sqrt(sum(psi_true[1,]^2))
  psi_true[2,] <- sin(4*grid*pi)
  psi_true[2,] <- psi_true[2,] / sqrt(sum(psi_true[2,]^2))
  ## not sure if we want this!
  psi2_true <- t(bs(grid, df = 2, degree = 1, intercept = TRUE))
  psi2_true <- psi2_true / sqrt(rowSums(psi2_true^2))



  fixef <- X_des %*% beta_true # fixed effects (n x L)
  set.seed(123)
  c_true <- rmvnorm(I, mean = rep(0, 2), sigma = diag(c(3, 1.5))) ## simulate score function, I x 2
  b_true <- c_true %*% psi_true # b true for random eff (I x L)

  id_visit_vec = rep(1:I, J_subj)
  X_df = # add a dataframe with ID indices for later sampling probs
    X_des %>%
    as.data.frame() %>%
    mutate(id = id_visit_vec) %>%
    group_by(id) %>%
    mutate(visit = row_number()) %>%
    ungroup() %>%
    mutate(index = row_number())


  ranef <- matrix(NA, nrow = n, ncol = L) # initialize random effects matrix
  ranef <- b_true[id_visit_vec, ]
  visit_effect_scale <- 0.3  # For example, reduce the effect by 70%
  if(!is.null(psi2_true)){ ## by default do not add subject-visit random deviation
    c2_true <- rmvnorm(n, mean = rep(0, 2), sigma = diag(c(3, 1.5)))
    ranef <- ranef + visit_effect_scale * (c2_true %*% psi2_true)
  }

  ranef <- sd(fixef)/sd(ranef)/SNR_B*ranef ## adjust for relative importance of random effects

  eta_true <- fixef + ranef # n x L

  Y_obs <- matrix(NA, n, L) # y observations: add noise to eta
  set.seed(123)
  if (family == "gaussian") {
    sd_signal <- sd(eta_true)
    Y_obs <- matrix(
      rnorm(n = nrow(eta_true) * L,
            mean = as.vector(t(eta_true)),
            sd = sd_signal/SNR_sigma),
      nrow = nrow(eta_true),
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "binomial") {
    p_true = plogis(eta_true)
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
    lam_true <- exp(eta_true)
    Y_obs <- matrix(
      rpois(n = nrow(lam_true) * L, lambda = as.vector(t(lam_true))),
      nrow = nrow(lam_true),
      ncol = L,
      byrow = TRUE
    )
  }


  nsim =
  sim_res <- list() ## store simulation results
  for (iter in 1:nsim) {
    set.seed(iter)
    x = try({
      data =
        simulate_survey_fosr_ml(X_des, # design matrix
                                X_df, # design matrix with ID and subject indices
                                eta_true, # true fixed effects with random eff
                                Y_obs,
                                J_subj, # num visits per person
                                I_n = 30, # subjects in each psu-strata combination
                                alpha = 2, # number PSU
                                L = 50, # dimension of the functional domain
                                num_strata = 30,  # total strata
                                sampled_strata = 10)


      ## pre-process simulated data
      # Y_mat <- data[, grepl('Y.', colnames(data))]
      # dat.fit <- data.frame(Y = Y_mat, data[, !grepl('Y.', colnames(data))])
      # model_formula = as.formula(paste0('Y~', 'X'))
      # out_index <- grep(paste0("^", model_formula[2]), names(data))

      betaTilde = map(.x = fit_types,
                      .f = get_betatilde,
                      data = data,
                      formula = Y ~ X + (1 | ID), argvals = NULL)

      map(.x = betaTilde,
          .f = \(x) x %>% as_tibble() %>%
            mutate(beta = c("beta_0", "beta_1")) %>%
            pivot_longer(cols = -beta) %>%
            mutate(L = as.numeric(name))) %>%
        bind_rows(.id = "id") %>%
        ggplot(aes(x = L, y = value, color = id)) +
        geom_line() +
        facet_wrap(.~beta)
      # betaTilde = get_betatilde(data, formula = Y ~ X + (1 | ID), argvals = NULL)
      betaHat = map(.f = get_betahat,
                    .x = betaTilde)

      if (parallel) {
        plan(multisession, workers = ncores)
        res = future_map2(
          .x = boot_types,
          .y = betaHat,
          .f = run_boots,
          data = data,
          family = family,
          num_boots = 20,
          set_seed = TRUE,
          seed = 2025,
          L = 50,
          .options = furrr_options(seed = TRUE)
        )
        # haven't done yet
        cis = future_map2(
          .x = res,
          .y = betaHat,
          .f = get_cis,
          smooth_for_ci = TRUE,
          .options = furrr_options(seed = TRUE)
        )
        plan(sequential)
      } else{
        res = map(
          .x = boot_types,
          .f = run_boots,
          data = data,
          betaHat = betaHat,
          num_boots = B,
          set_seed = TRUE,
          seed = 2025,
          family = "gaussian",
          L = 50
        )

        cis = map(
          .x = res,
          .f = get_cis,
          betaHat = betaHat,
          .progress = TRUE
        )
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

