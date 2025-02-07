source(here::here("R", "00_data_gen_function.R"))
source(here::here("R", "01_sim_functions.R"))
nsim <- 1
B <- 500
parallel = TRUE
library(future)
library(furrr)
ncores = parallel::detectCores() - 1
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'Preston', 'no_design')

sim_res <- list() ## store simulation results
for(iter in 1:nsim) {
  # iter = 1
  set.seed(iter)
  data <- simulate_survey_fosr()

  ## pre-process simulated data
  Y_mat <- data[, grepl('Y.', colnames(data))]
  dat.fit <- data.frame(Y = Y_mat, data[, !grepl('Y.', colnames(data))])
  model_formula = as.formula(paste0('Y~', 'X'))
  out_index <- grep(paste0("^", model_formula[2]), names(data))


  betaTilde = get_betatilde(data)
  betaHat = get_betahat(betaTilde)

  if (parallel) {
    plan(multisession, workers = ncores)
    res = future_map(
      .x = boot_types,
      .f = run_boots,
      data = data,
      betaHat = betaHat,
      num_boots = B,
      set_seed = TRUE,
      seed = 2025,
      L = 50,
      out_index = out_index
    )

    cis = future_map(.x = res,
                     .f = get_cis,
                     betaHat = betaHat)
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
      L = 50,
      out_index = out_index
    )

    cis = map(.x = res,
              .f = get_cis,
              betaHat = betaHat,
              .progress = TRUE)
  }

  stats = map2(
    .x = cis,
    .y = boot_types,
    .f = get_coverage_stats,
    beta_true = beta_true
  ) %>%
    list_rbind() %>%
    mutate(n = nrow(data))

  sim_res[[iter]] <- stats
  rm(stats)
  rm(data)
}



#
# ## Clean results
# ################################################################################
# ## Obtain MISE, coverage, and computing time
# ################################################################################
# ## MISE
# MISE_lfosr3s <- lapply(sim_local, '[[', 1) %>% bind_rows()
# colMeans(MISE_lfosr3s)
#
# ## coverage
# ### joint
# cover_lfosr3s <- t(lapply(sim_local, '[[', 2) %>% bind_cols())
# ### point-wise
# cover_lfosr3s_pw <- array(NA, dim = c(nsim, L, p))
# for(i in 1:nsim){
#   cover_lfosr3s_pw[i,,1] <- sim_local[[i]]$cover_lfosr3s_pw[2,] # intercept
#   cover_lfosr3s_pw[i,,2] <- sim_local[[i]]$cover_lfosr3s_pw[2,] # fixed effects
# }
# ### print results
# print(paste0("lfosr3s coverage (joint CI): ", colMeans(cover_lfosr3s)[2]))
# print(paste0("lfosr3s coverage (pointwise CI): ", apply(cover_lfosr3s_pw, 3, mean)[2]))
#
# ## time
# time_lfosr3s <- lapply(sim_local, '[[', 4) %>% bind_rows()
# print(apply(time_lfosr3s, 2, median))
