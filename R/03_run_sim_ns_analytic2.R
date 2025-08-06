### generate PA data


library(here)
library(devtools)
library(tictoc)
library(fastFMM)
library(haven)
library(dplyr)
library(survey)
library(progress)
library(lme4)
library(mgcv)
library(ggplot2)
library(furrr)
library(future)
library(gridExtra)
library(refund)
library(tidyverse)
library(parallel)
source(here::here("R", "sim_functions_unwt.R"))
ncores_total = parallelly::availableCores() - 1
force = TRUE


get_var = function(data, family = "gaussian",
                   model_formula = as.formula(paste0('Y~', 'X'))){


  out_index <- grep(paste0("^", model_formula[2]), names(data)) # indices that start with the outcome name
  if(length(out_index) != 1){ # functional observations stored in multiple columns
    L <- length(out_index)
  }else{ # observations stored as a matrix in one column using the I() function
    L <- ncol(data[,out_index])
  }
  argvals <- 1:L

  unimm <- function(l) {
    data$Yl <- unclass(data[, out_index][, l])
    fit_uni <- suppressMessages(glm(
      formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
      family = family,
      control = glm.control(maxit = 5000),
      data = data
    ))
    sd = summary(fit_uni)$coef[,2]
    return(list(sd = sd))
  }

  massmm <- lapply(argvals, unimm)


  # Obtain betaTilde, fixed effects estimates
  sds <- t(do.call(rbind, lapply(massmm, '[[', 1)))
  colnames(sds) <- argvals
  return(sds)
}


get_cis_analytic = function(betaHat, betaTilde, L = 50, data, nknots = 35){
  p = nrow(betaHat) # number of beta coefs
  designmat = cbind(rep(1, nrow(data)), data$X) # design matrix X
  Y_mat = data %>% # Y matrix
    select(starts_with("Y")) %>%
    as.matrix()

  Y_pred = designmat %*% betaTilde # y = x beta
  resid_mat = Y_mat - Y_pred # residual matrix
  resid_mat_smooth <- matrix(NA, nrow = nrow(resid_mat), ncol = ncol(resid_mat))

  argvals = 1:ncol(resid_mat)
  # for (i in 1:nrow(resid_mat)) {
  #   # Fit GAM to residual curve for subject i
  #   gam_fit <- gam(resid_mat[i, ] ~ s(argvals, bs = "tp", k = nknots + 1),
  #                  method = "REML")
  #   resid_mat_smooth[i, ] <- fitted(gam_fit)
  # }

  for (i in 1:nrow(resid_mat)) {
    # Fit GAM to residual curve for subject i
    # resid_mat_smooth[i, ] <- stats::smooth.spline(x = argvals, y = resid_mat[i,], cv = TRUE)$y
    resid_mat_smooth[i, ] <- stats::smooth.spline(x = argvals, y = resid_mat[i,], spar = .75)$y
  }

  # resid_mat %>%
  #   as.data.frame()  %>%
  #   mutate(rn = row_number()) %>%
  #   slice(1:10) %>%
  #   pivot_longer(cols = starts_with("Y")) %>%
  #   mutate(n = as.numeric(sub(".*Y", "", name))) %>%
  #   ggplot(aes(x = n, y = value, color = factor(rn))) +
  #   geom_line() +
  #   theme(legend.position = "none")
  #
  resid_mat_smooth %>%
    as.data.frame()  %>%
    mutate(type = "smooth") %>%
    mutate(rn = row_number()) %>%
    slice(1:10) %>%
    pivot_longer(cols = starts_with("V")) %>%
    mutate(n = as.numeric(sub(".*V", "", name))) %>%
    ggplot(aes(x = n, y = value, color = factor(rn))) +
    geom_line() +
    theme(legend.position = "none")

  R_hat = matrix(NA, nrow = L, ncol = L) # covariance matrix of Ys
  for (l1 in 1:L) {
    for (l2 in 1:L) {
      if (l1 == l2) {
        R_hat[l1, l2] = (1 / (n - p)) * sum((resid_mat_smooth[,l1]^2)) # if l1 = l2, use residual variance as estimate
      } else{
        R_hat[l1, l2] = mean(resid_mat_smooth[,l1] * resid_mat_smooth[,l2]) # method of moments estimator
      }
    }
  }


  XTX_inv <- solve(t(designmat) %*% designmat) # precompute
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  for (l1 in 1:L) {
    for (l2 in 1:L) {
      betaHat.var[l1, l2, ] <- diag((XTX_inv %*% t(designmat) * R_hat[l1, l2]) %*% designmat %*% XTX_inv)
    }
  }

  # smooth
  betaHat.var.smooth <- array(NA, dim = dim(betaHat.var))

  for (r in 1:nrow(betaHat)) {
    betaHat.var.smooth[, , r] <- refund::fbps(betaHat.var[, , r])$Yhat
    # Ensure positive definite
    # betaHat.var[, , r] <- pmax(0, betaHat.var[, , r])
    betaHat.var.smooth[, , r] <- eigenval_trim(betaHat.var.smooth[, , r])
  }

  qn <- rep(0, length = nrow(betaHat))
  N <- 10000 ## sample size in simulation-based approach
  for (i in 1:length(qn)) {
    Sigma <- betaHat.var.smooth[, , i]
    x_sample <- mvtnorm::rmvnorm(N, mean = betaHat[i, ], sigma = Sigma)
    un <- rep(NA, N)
    for (j in 1:N) {
      un[j] <- max(abs((x_sample[j, ] - betaHat[i, ]) / sqrt(diag(Sigma))))
    }
    qn[i] <- quantile(un, 0.95)
  }

  return(list(
    betaHat = betaHat,
    betaHat.var = betaHat.var.smooth,
    qn = qn
  ))
}



source(here::here("R", "utils.R"))
sim_settings = read_rds(here::here("sim_data", "sim_settings_nonsurvey.rds"))


keep_folds =
  sim_settings %>%
  filter(num_boots == 500) %>%
  group_by(family, scenario, L, n) %>%
  slice(1) %>%
  ungroup() %>%
  pull(fold)

sim_settings =
  sim_settings %>%
  filter(num_boots == 500) %>%
  filter(family == "gaussian" | fold %in% keep_folds)

sim_settings %>%
  filter(family != "gaussian") %>%
  pull(fold) %>%
  paste(collapse = ",")

ifold = get_fold()
ifold = 1
# ifold = 111
fname = paste0("sim_fold_", sprintf("%03d", ifold), "_analytic_smooth75.rds")


  tmp_settings = sim_settings %>% filter(fold == ifold) # set settings for this simulation


  scen = tmp_settings$scenario
  sigma = tmp_settings$sigma
  family = tmp_settings$family
  B = tmp_settings$num_boots
  samp_size = tmp_settings$n
  nsims = 500
  len = tmp_settings$L
  # if(samp_size > 1000 || len > 1000){
  #   run_bayes = FALSE
  # } else{
  #   run_bayes = TRUE
  # }
  run_bayes = FALSE
  # if(family == "gaussian"){ # if gaussian, we don't need to bootstrap so don't need inner parallelization
    ncores_outer = ncores_total
  # } else{
  # ncores_outer = ceiling(ncores_total / 2)
  # ncores_inner = ncores_total - ncores_outer

  plan(multisession, workers = ncores_outer)
  sim_res <- future_map(1:nsims, function(iter) {
    x = try({
      data = generate_population(
        n = samp_size,
        L = len,
        scenario = scen,
        family = family,
        seed = iter,
        SNR_sigma = sigma
      )
      sample_data = data$sample_data

      #### fui
      # if(family != "gaussian") {
      tic()

      beta_tilde =
        sample_data %>%
        get_betatilde(., family = family)


      b_hat = beta_tilde %>%
        get_betahat(., L = len)

      ci = get_cis_analytic(
        betaHat = b_hat,
        betaTilde = beta_tilde,
        L = len,
        data = sample_data
      )
      elapsed_time = toc(quiet = TRUE)
      stats_fui = get_coverage_stats_fui(ci, beta_true = data$beta_true, L = len) %>%
        mutate(time = unname(elapsed_time$toc - elapsed_time$tic))


      elapsed_time <- NULL


      ret =
        stats_fui %>% mutate(type = "FUI")
      rm(stats_fui, elapsed_time, sample_data, b_hat, ci, data)

    # print(iter)
    return(ret)
    })
    x
    },
.options = furrr_options(seed = TRUE),
.progress = TRUE)

  # print(sim_res)

  result =
    sim_res %>%
    keep(., is.data.frame) %>%
    list_rbind(names_to = "iter") %>%
    select(type, iter, var, everything())

  # check
  # result %>%
  #   ggplot(aes(x = var, y = time, color = type)) +
  #   geom_boxplot()
  #
  # result %>%
  #   pivot_longer(cols = c(MISE, mean_pw_se, cover_pw, time)) %>%
  #   ggplot(aes(x = var, y = value, color = type)) +
  #   geom_boxplot(outlier.shape = NA)+
  #   facet_wrap(.~name, scales = "free_y")

  if (!dir.exists(here::here("results", "simulations", "non_survey"))) {
    dir.create(here::here("results", "simulations", "non_survey"), recursive = TRUE)
  }

  write_rds(result,
            here::here("results", "simulations", "non_survey", fname))

plan(sequential)

get_cis_analytic = function(betaHat, betaTilde, L = 50, data, nknots = 35){
  p = nrow(betaHat) # number of beta coefs
  designmat = cbind(rep(1, nrow(data)), data$X) # design matrix X
  Y_mat = data %>% # Y matrix
    select(starts_with("Y")) %>%
    as.matrix()

  Y_pred = designmat %*% betaTilde # y = x beta
  resid_mat = Y_mat - Y_pred # residual matrix
  resid_mat_smooth <- matrix(NA, nrow = nrow(resid_mat), ncol = ncol(resid_mat))

  argvals = 1:ncol(resid_mat)
  # for (i in 1:nrow(resid_mat)) {
  #   # Fit GAM to residual curve for subject i
  #   gam_fit <- gam(resid_mat[i, ] ~ s(argvals, bs = "tp", k = nknots + 1),
  #                  method = "REML")
  #   resid_mat_smooth[i, ] <- fitted(gam_fit)
  # }

  for (i in 1:nrow(resid_mat)) {
    # Fit GAM to residual curve for subject i
    # resid_mat_smooth[i, ] <- stats::smooth.spline(x = argvals, y = resid_mat[i,], cv = TRUE)$y
    resid_mat_smooth[i, ] <- stats::smooth.spline(x = argvals, y = resid_mat[i,], spar = .5)$y
  }



  R_hat = matrix(NA, nrow = L, ncol = L) # covariance matrix of Ys
  for (l1 in 1:L) {
    for (l2 in 1:L) {
      if (l1 == l2) {
        R_hat[l1, l2] = (1 / (n - p)) * sum((resid_mat_smooth[,l1]^2)) # if l1 = l2, use residual variance as estimate
      } else{
        R_hat[l1, l2] = mean(resid_mat_smooth[,l1] * resid_mat_smooth[,l2]) # method of moments estimator
      }
    }
  }


  XTX_inv <- solve(t(designmat) %*% designmat) # precompute
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  for (l1 in 1:L) {
    for (l2 in 1:L) {
      betaHat.var[l1, l2, ] <- diag((XTX_inv %*% t(designmat) * R_hat[l1, l2]) %*% designmat %*% XTX_inv)
    }
  }

  # smooth
  betaHat.var.smooth <- array(NA, dim = dim(betaHat.var))

  for (r in 1:nrow(betaHat)) {
    betaHat.var.smooth[, , r] <- refund::fbps(betaHat.var[, , r])$Yhat
    # Ensure positive definite
    # betaHat.var[, , r] <- pmax(0, betaHat.var[, , r])
    betaHat.var.smooth[, , r] <- eigenval_trim(betaHat.var.smooth[, , r])
  }

  qn <- rep(0, length = nrow(betaHat))
  N <- 10000 ## sample size in simulation-based approach
  for (i in 1:length(qn)) {
    Sigma <- betaHat.var.smooth[, , i]
    x_sample <- mvtnorm::rmvnorm(N, mean = betaHat[i, ], sigma = Sigma)
    un <- rep(NA, N)
    for (j in 1:N) {
      un[j] <- max(abs((x_sample[j, ] - betaHat[i, ]) / sqrt(diag(Sigma))))
    }
    qn[i] <- quantile(un, 0.95)
  }

  return(list(
    betaHat = betaHat,
    betaHat.var = betaHat.var.smooth,
    qn = qn
  ))
}



source(here::here("R", "utils.R"))
sim_settings = read_rds(here::here("sim_data", "sim_settings_nonsurvey.rds"))


keep_folds =
  sim_settings %>%
  filter(num_boots == 500) %>%
  group_by(family, scenario, L, n) %>%
  slice(1) %>%
  ungroup() %>%
  pull(fold)

sim_settings =
  sim_settings %>%
  filter(num_boots == 500) %>%
  filter(family == "gaussian" | fold %in% keep_folds)

sim_settings %>%
  filter(family != "gaussian") %>%
  pull(fold) %>%
  paste(collapse = ",")

ifold = get_fold()
ifold = 1
# ifold = 111
fname = paste0("sim_fold_", sprintf("%03d", ifold), "_analytic_smooth50.rds")


tmp_settings = sim_settings %>% filter(fold == ifold) # set settings for this simulation


scen = tmp_settings$scenario
sigma = tmp_settings$sigma
family = tmp_settings$family
B = tmp_settings$num_boots
samp_size = tmp_settings$n
nsims = 500
len = tmp_settings$L
# if(samp_size > 1000 || len > 1000){
#   run_bayes = FALSE
# } else{
#   run_bayes = TRUE
# }
run_bayes = FALSE
# if(family == "gaussian"){ # if gaussian, we don't need to bootstrap so don't need inner parallelization
ncores_outer = ncores_total
# } else{
# ncores_outer = ceiling(ncores_total / 2)
# ncores_inner = ncores_total - ncores_outer

plan(multisession, workers = ncores_outer)
sim_res <- future_map(1:nsims, function(iter) {
  x = try({
    data = generate_population(
      n = samp_size,
      L = len,
      scenario = scen,
      family = family,
      seed = iter,
      SNR_sigma = sigma
    )
    sample_data = data$sample_data

    #### fui
    # if(family != "gaussian") {
    tic()

    beta_tilde =
      sample_data %>%
      get_betatilde(., family = family)


    b_hat = beta_tilde %>%
      get_betahat(., L = len)

    ci = get_cis_analytic(
      betaHat = b_hat,
      betaTilde = beta_tilde,
      L = len,
      data = sample_data
    )
    elapsed_time = toc(quiet = TRUE)
    stats_fui = get_coverage_stats_fui(ci, beta_true = data$beta_true, L = len) %>%
      mutate(time = unname(elapsed_time$toc - elapsed_time$tic))


    elapsed_time <- NULL


    ret =
      stats_fui %>% mutate(type = "FUI")
    rm(stats_fui, elapsed_time, sample_data, b_hat, ci, data)

    # print(iter)
    return(ret)
  })
  x
},
.options = furrr_options(seed = TRUE),
.progress = TRUE)

# print(sim_res)

result =
  sim_res %>%
  keep(., is.data.frame) %>%
  list_rbind(names_to = "iter") %>%
  select(type, iter, var, everything())

# check
# result %>%
#   ggplot(aes(x = var, y = time, color = type)) +
#   geom_boxplot()
#
# result %>%
#   pivot_longer(cols = c(MISE, mean_pw_se, cover_pw, time)) %>%
#   ggplot(aes(x = var, y = value, color = type)) +
#   geom_boxplot(outlier.shape = NA)+
#   facet_wrap(.~name, scales = "free_y")

if (!dir.exists(here::here("results", "simulations", "non_survey"))) {
  dir.create(here::here("results", "simulations", "non_survey"), recursive = TRUE)
}

write_rds(result,
          here::here("results", "simulations", "non_survey", fname))

plan(sequential)
