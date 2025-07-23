## joint CI visualization
# snr_eps snr_b    In   BRR
# <dbl> <dbl> <dbl> <dbl>
#   0.1   0.1   250 0.870
# 5     0.1   250 0.930
#  0.1   0.5   250 0.75


# visualize the data
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
library(patchwork)


nsim = 200
B = 500
force = TRUE

ncores = parallelly::availableCores() - 1
fit_types = c('weighted')
boot_types = c('BRR')

options(survey.lonely.psu = "adjust")

# just use PSU RE = 0.5

temp = tibble(
  psu_re = 0.5,
  len = 50,
  snr_b = 0.5,
  snr_eps = 0.1,
  scen = 1,
  family = "gaussian",
  In = 250
)

strata_sigma = 0.05
psu_sigma = sqrt(strata_sigma ^ 2 * temp$psu_re)

lst = generate_superpopulation(
  scenario = temp$scen,
  family = temp$family,
  I = 10e6,
  L = temp$len,
  psu_sigma = psu_sigma,
  snr_b = temp$snr_b,
  snr_eps = temp$snr_eps
)

L = 50
nknots_min = NULL
nknots_min_cov = 35
mult_fac = 1.3
nknots <- min(round(L / 2), nknots_min)
nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
nknots_fpca <- min(round(L / 2), 35)
argvals = 1:L

pdf(here::here("joint_ci_13.pdf"))
for(iter in 1:nsim){
  x = try({
    set.seed(iter)
    data <- sample_from_population_wor(
      X_des = lst$X_des,
      Y_obs = lst$Y_obs,
      L = temp$len,
      ## dimension of the functional domain
      I_n = temp$In,
      num_strata = 30,
      stratum_assignments = lst$stratum_assignments,
      psu_assignments = lst$psu_assignments,
      dirichlet_probs = lst$dirichlet_probs,
      seed = iter
    )

    beta_true = lst$beta_true
    betaTilde = map(
      .x = fit_types,
      .f = get_betatilde,
      data = data,
      family = temp$family
    )

    betaHat = map(.x = betaTilde, get_betahat, L = temp$len)


    plan(multisession, workers = ncores)
    res = future_map2(
      .x = boot_types,
      .y = betaHat,
      .f = run_boots,
      data = data,
      family = temp$family,
      num_boots = B,
      set_seed = TRUE,
      seed = 2025 + iter,
      L = temp$len,
      samp_stages = c("PPSWOR", "Poisson"),
      .options = furrr_options(seed = TRUE),
      .progress = TRUE
    )
    plan(sequential)



    res_x = res[[1]][2, , ]

    btrue_df =
      beta_true %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    bhat_df =
      betaHat[[1]] %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    btilde_df =
      betaTilde[[1]] %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    # Plot the raw estimates for 10 repetitions
    p1 =
      res_x %>%
      as_tibble() %>%
      mutate(l = row_number()) %>%
      pivot_longer(cols = starts_with("V"), names_to = "rep", values_to = "value")  %>%
      mutate(id = sub(".*V", "", rep) %>% as.numeric) %>%
      filter(id <= 10) %>%
      ggplot(aes(x = l, y = value, color = factor(rep))) +
      theme(legend.position = "none") +
      geom_line() +
      # Add labeled lines for custom legend
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      geom_line(data = btilde_df,
                aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      # Manual legend
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Raw Estimate" = "blue")) +
      labs(x = "Functional Domain", y = "Coefficient", color = "Line Type",
           title = "Raw Beta Coef Estimates for 10 Repetitions") +
      theme_light()

    betaHat_boot <- array(NA, dim = c(2, nrow(res_x), ncol(res_x)))

    # smooth the bootstrapped estimates
    for (b in 1:B) {
      betaHat_boot[, , b] <- t(apply(res[[1]][, , b], 1, function(x)
        gam(x ~ s(
          argvals, bs = "tp", k = (nknots + 1)
        ), method = "GCV.Cp")$fitted.values))
    }

    boot_df =
      betaHat_boot[2,,] %>%
      as_tibble() %>%
      mutate(l = row_number()) %>%
      pivot_longer(cols = starts_with("V"), names_to = "rep", values_to = "value")  %>%
      mutate(id = sub(".*V", "", rep) %>% as.numeric) %>%
      filter(id <= 10)
    p2 =
      boot_df %>%
      ggplot(aes(x = l, y = value, color = factor(rep))) +
      theme(legend.position = "none") +
      geom_line() +
      # Add labeled lines for custom legend
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      geom_line(data = btilde_df,
                aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      # Manual legend
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Raw Estimate" = "blue")) +
      labs(x = "Functional Domain", y = "Coefficient", color = "Line Type",
           title = "Smoothed Beta Coef Estimates for 10 Repetitions") +
      theme_light()


    betaHat.var <- array(NA, dim = c(L, L, 2))

    for (r in 1:2) {
      betaHat.var[, , r] <- mult_fac * var(t(betaHat_boot[r, , ]))
    }
    var_x = betaHat.var[,,2]

    bhat_df2 = betaHat[[1]][2,] %>%
      as_tibble() %>%
      rename(beta = value) %>%
      mutate(l = row_number())

    N <- 10000
    qn <- rep(0, length = nrow(betaHat[[1]]))
    set.seed(456)
    i = 2

    est_bs <- t(betaHat_boot[i, , ])


    fit_fpca <- suppressWarnings(refund::fpca.face(est_bs, knots = nknots_fpca)) # suppress sqrt(Eigen$values) NaNs
    ## extract estimated eigenfunctions/eigenvalues
    phi <- fit_fpca$efunctions
    lambda <- fit_fpca$evalues
    K <- length(fit_fpca$evalues)

    ## simulate random coefficients
    theta <- matrix(stats::rnorm(N * K), nrow = N, ncol = K) # generate independent standard normals
    if (K == 1) {
      theta <- theta * sqrt(lambda) # scale to have appropriate variance
      X_new <- tcrossprod(theta, phi) # simulate new functions
    } else{
      theta <- theta %*% diag(sqrt(lambda)) # scale to have appropriate variance
      X_new <- tcrossprod(theta, phi) # simulate new functions
    }
    x_sample <- X_new + t(fit_fpca$mu %o% rep(1, N)) # add back in the mean function
    Sigma_sd <- Rfast::colVars(x_sample, std = TRUE, na.rm = FALSE) # standard deviation: apply(x_sample, 2, sd)
    x_mean <- colMeans(est_bs)
    un <- rep(NA, N)
    for (j in 1:N) {
      un[j] <- max(abs((x_sample[j, ] - x_mean) / Sigma_sd))
    }
    qn[i] <- stats::quantile(un, 0.95)
    p3 =
      un %>%
      as_tibble() %>%
      ggplot(aes(x = value)) +
      geom_histogram(color = "black") +
      geom_vline(aes(xintercept = stats::quantile(un, 0.95)), color = "red") +
      labs(x = "Value", y= "Count", title = "Histogram of Max Statisics") +
      theme_light()


    cma = qn[2]
    var_x_df =
      diag(var_x) %>%
      as_tibble() %>%
      mutate(l = row_number()) %>%
      left_join(bhat_df2) %>%
      mutate(min = beta - 1.96 * sqrt(value),
             max = beta + 1.96 * sqrt(value),
             cma_min  = beta - cma * sqrt(value),
             cma_max = beta + cma * sqrt(value))


    p4 = var_x_df %>%
      ggplot() +
      geom_ribbon(aes(x = l, ymin = cma_min, ymax = cma_max), fill = "darkgrey") +
      geom_ribbon(aes(x = l, ymin = min, ymax = max), fill = "lightgrey") +
      theme_light() +
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      # geom_line(data = btilde_df %>% filter(name == "X"),
      #           aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Raw Estimate" = "blue"))


    # get some coverage stats
    margin_pw = 1.96 * sqrt(diag(betaHat.var[,,2]))
    margin_joint = cma * sqrt(diag(betaHat.var[,,2]))
    upper_j = betaHat[[1]][2,] + margin_joint
    lower_j = betaHat[[1]][2,] - margin_joint
    upper = betaHat[[1]][2,] + margin_pw
    lower = betaHat[[1]][2,] - margin_pw
    p = 2
    true = beta_true[p,]
    covered <- (lower_j < true) & (upper_j > true)
    covered_pw <- (lower < true) & (upper > true)
    cover_joint =  all(covered)  # Joint coverage check
    cover_pw = sum(covered_pw)


    p4 = var_x_df %>%
      ggplot() +
      geom_ribbon(aes(x = l, ymin = cma_min, ymax = cma_max), fill = "darkgrey") +
      geom_ribbon(aes(x = l, ymin = min, ymax = max), fill = "lightgrey") +
      theme_light() +
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      # geom_line(data = btilde_df %>% filter(name == "X"),
      #           aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Raw Estimate" = "blue"))  +
      labs(title = "Confidence Intervals",color = "Line Type", x = "Functional Domain", y = "Coefficient",
           subtitle = paste0("Joint covered = ", cover_joint, " PW covered = ", cover_pw))

    print((p1 + p2 + plot_layout(guides = "collect", axes = "collect"))/ (p3 + p4))
  })
  rm(x)
}

dev.off()

L = 50
nknots_min = NULL
nknots_min_cov = 35
mult_fac = 1.2
nknots <- min(round(L / 2), nknots_min)
nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
nknots_fpca <- min(round(L / 2), 35)
argvals = 1:L

pdf(here::here("joint_ci_12.pdf"))
for(iter in 1:nsim){
  x = try({
    set.seed(iter)
    data <- sample_from_population_wor(
      X_des = lst$X_des,
      Y_obs = lst$Y_obs,
      L = temp$len,
      ## dimension of the functional domain
      I_n = temp$In,
      num_strata = 30,
      stratum_assignments = lst$stratum_assignments,
      psu_assignments = lst$psu_assignments,
      dirichlet_probs = lst$dirichlet_probs,
      seed = iter
    )

    beta_true = lst$beta_true
    betaTilde = map(
      .x = fit_types,
      .f = get_betatilde,
      data = data,
      family = temp$family
    )

    betaHat = map(.x = betaTilde, get_betahat, L = temp$len)


    plan(multisession, workers = ncores)
    res = future_map2(
      .x = boot_types,
      .y = betaHat,
      .f = run_boots,
      data = data,
      family = temp$family,
      num_boots = B,
      set_seed = TRUE,
      seed = 2025 + iter,
      L = temp$len,
      samp_stages = c("PPSWOR", "Poisson"),
      .options = furrr_options(seed = TRUE),
      .progress = TRUE
    )
    plan(sequential)



    res_x = res[[1]][2, , ]

    btrue_df =
      beta_true %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    bhat_df =
      betaHat[[1]] %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    btilde_df =
      betaTilde[[1]] %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    # Plot the raw estimates for 10 repetitions
    p1 =
      res_x %>%
      as_tibble() %>%
      mutate(l = row_number()) %>%
      pivot_longer(cols = starts_with("V"), names_to = "rep", values_to = "value")  %>%
      mutate(id = sub(".*V", "", rep) %>% as.numeric) %>%
      filter(id <= 10) %>%
      ggplot(aes(x = l, y = value, color = factor(rep))) +
      theme(legend.position = "none") +
      geom_line() +
      # Add labeled lines for custom legend
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      geom_line(data = btilde_df,
                aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      # Manual legend
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Raw Estimate" = "blue")) +
      labs(x = "Functional Domain", y = "Coefficient", color = "Line Type",
           title = "Raw Beta Coef Estimates for 10 Repetitions") +
      theme_light()

    betaHat_boot <- array(NA, dim = c(2, nrow(res_x), ncol(res_x)))

    # smooth the bootstrapped estimates
    for (b in 1:B) {
      betaHat_boot[, , b] <- t(apply(res[[1]][, , b], 1, function(x)
        gam(x ~ s(
          argvals, bs = "tp", k = (nknots + 1)
        ), method = "GCV.Cp")$fitted.values))
    }

    boot_df =
      betaHat_boot[2,,] %>%
      as_tibble() %>%
      mutate(l = row_number()) %>%
      pivot_longer(cols = starts_with("V"), names_to = "rep", values_to = "value")  %>%
      mutate(id = sub(".*V", "", rep) %>% as.numeric) %>%
      filter(id <= 10)
    p2 =
      boot_df %>%
      ggplot(aes(x = l, y = value, color = factor(rep))) +
      theme(legend.position = "none") +
      geom_line() +
      # Add labeled lines for custom legend
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      geom_line(data = btilde_df,
                aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      # Manual legend
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Raw Estimate" = "blue")) +
      labs(x = "Functional Domain", y = "Coefficient", color = "Line Type",
           title = "Smoothed Beta Coef Estimates for 10 Repetitions") +
      theme_light()


    betaHat.var <- array(NA, dim = c(L, L, 2))

    for (r in 1:2) {
      betaHat.var[, , r] <- mult_fac * var(t(betaHat_boot[r, , ]))
    }
    var_x = betaHat.var[,,2]

    bhat_df2 = betaHat[[1]][2,] %>%
      as_tibble() %>%
      rename(beta = value) %>%
      mutate(l = row_number())

    N <- 10000
    qn <- rep(0, length = nrow(betaHat[[1]]))
    set.seed(456)
    i = 2

    est_bs <- t(betaHat_boot[i, , ])


    fit_fpca <- suppressWarnings(refund::fpca.face(est_bs, knots = nknots_fpca)) # suppress sqrt(Eigen$values) NaNs
    ## extract estimated eigenfunctions/eigenvalues
    phi <- fit_fpca$efunctions
    lambda <- fit_fpca$evalues
    K <- length(fit_fpca$evalues)

    ## simulate random coefficients
    theta <- matrix(stats::rnorm(N * K), nrow = N, ncol = K) # generate independent standard normals
    if (K == 1) {
      theta <- theta * sqrt(lambda) # scale to have appropriate variance
      X_new <- tcrossprod(theta, phi) # simulate new functions
    } else{
      theta <- theta %*% diag(sqrt(lambda)) # scale to have appropriate variance
      X_new <- tcrossprod(theta, phi) # simulate new functions
    }
    x_sample <- X_new + t(fit_fpca$mu %o% rep(1, N)) # add back in the mean function
    Sigma_sd <- Rfast::colVars(x_sample, std = TRUE, na.rm = FALSE) # standard deviation: apply(x_sample, 2, sd)
    x_mean <- colMeans(est_bs)
    un <- rep(NA, N)
    for (j in 1:N) {
      un[j] <- max(abs((x_sample[j, ] - x_mean) / Sigma_sd))
    }
    qn[i] <- stats::quantile(un, 0.95)
    p3 =
      un %>%
      as_tibble() %>%
      ggplot(aes(x = value)) +
      geom_histogram(color = "black") +
      geom_vline(aes(xintercept = stats::quantile(un, 0.95)), color = "red") +
      labs(x = "Value", y= "Count", title = "Histogram of Max Statisics") +
      theme_light()


    cma = qn[2]
    var_x_df =
      diag(var_x) %>%
      as_tibble() %>%
      mutate(l = row_number()) %>%
      left_join(bhat_df2) %>%
      mutate(min = beta - 1.96 * sqrt(value),
             max = beta + 1.96 * sqrt(value),
             cma_min  = beta - cma * sqrt(value),
             cma_max = beta + cma * sqrt(value))


    p4 = var_x_df %>%
      ggplot() +
      geom_ribbon(aes(x = l, ymin = cma_min, ymax = cma_max), fill = "darkgrey") +
      geom_ribbon(aes(x = l, ymin = min, ymax = max), fill = "lightgrey") +
      theme_light() +
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      # geom_line(data = btilde_df %>% filter(name == "X"),
      #           aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Raw Estimate" = "blue"))


    # get some coverage stats
    margin_pw = 1.96 * sqrt(diag(betaHat.var[,,2]))
    margin_joint = cma * sqrt(diag(betaHat.var[,,2]))
    upper_j = betaHat[[1]][2,] + margin_joint
    lower_j = betaHat[[1]][2,] - margin_joint
    upper = betaHat[[1]][2,] + margin_pw
    lower = betaHat[[1]][2,] - margin_pw
    p = 2
    true = beta_true[p,]
    covered <- (lower_j < true) & (upper_j > true)
    covered_pw <- (lower < true) & (upper > true)
    cover_joint =  all(covered)  # Joint coverage check
    cover_pw = sum(covered_pw)


    p4 = var_x_df %>%
      ggplot() +
      geom_ribbon(aes(x = l, ymin = cma_min, ymax = cma_max), fill = "darkgrey") +
      geom_ribbon(aes(x = l, ymin = min, ymax = max), fill = "lightgrey") +
      theme_light() +
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      # geom_line(data = btilde_df %>% filter(name == "X"),
      #           aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Raw Estimate" = "blue"))  +
      labs(title = "Confidence Intervals",color = "Line Type", x = "Functional Domain", y = "Coefficient",
           subtitle = paste0("Joint covered = ", cover_joint, " PW covered = ", cover_pw))

    print((p1 + p2 + plot_layout(guides = "collect", axes = "collect"))/ (p3 + p4))
  })
  rm(x)
}

dev.off()


pdf(here::here("joint_ci_13.pdf"))
for(iter in 1:nsim){
  x = try({
    set.seed(iter)
    data <- sample_from_population_wor(
      X_des = lst$X_des,
      Y_obs = lst$Y_obs,
      L = temp$len,
      ## dimension of the functional domain
      I_n = temp$In,
      num_strata = 30,
      stratum_assignments = lst$stratum_assignments,
      psu_assignments = lst$psu_assignments,
      dirichlet_probs = lst$dirichlet_probs,
      seed = iter
    )

    beta_true = lst$beta_true
    betaTilde = map(
      .x = fit_types,
      .f = get_betatilde,
      data = data,
      family = temp$family
    )

    betaHat = map(.x = betaTilde, get_betahat, L = temp$len)


    plan(multisession, workers = ncores)
    res = future_map2(
      .x = boot_types,
      .y = betaHat,
      .f = run_boots,
      data = data,
      family = temp$family,
      num_boots = B,
      set_seed = TRUE,
      seed = 2025,
      L = temp$len,
      samp_stages = c("PPSWOR", "Poisson"),
      .options = furrr_options(seed = TRUE),
      .progress = TRUE
    )
    plan(sequential)



    res_x = res[[1]][2, , ]

    res_est = map(.x = res,
                  .f = function(bt) { apply(bt, 1, rowMeans) %>% t()})

    betaHat2 = res_est[[1]]

    btrue_df =
      beta_true %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    bhat_df =
      betaHat[[1]] %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    bhat2_df =
      betaHat2 %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    btilde_df =
      betaTilde[[1]] %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    # Plot the raw estimates for 10 repetitions
    p1 =
      res_x %>%
      as_tibble() %>%
      mutate(l = row_number()) %>%
      pivot_longer(cols = starts_with("V"), names_to = "rep", values_to = "value")  %>%
      mutate(id = sub(".*V", "", rep) %>% as.numeric) %>%
      filter(id <= 10) %>%
      ggplot(aes(x = l, y = value, color = factor(rep))) +
      theme(legend.position = "none") +
      geom_line() +
      # Add labeled lines for custom legend
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      geom_line(data = bhat2_df,
                aes(x = l, y = value, color = "Smoothed 2 Estimate"), linewidth = 1) +
      geom_line(data = btilde_df,
                aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      # Manual legend
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Smoothed 2 Estimate" = "grey",
                                    "Raw Estimate" = "blue")) +
      labs(x = "Functional Domain", y = "Coefficient", color = "Line Type",
           title = "Raw Beta Coef Estimates for 10 Repetitions") +
      theme_light()

    betaHat_boot <- array(NA, dim = c(2, nrow(res_x), ncol(res_x)))

    # smooth the bootstrapped estimates
    for (b in 1:B) {
      betaHat_boot[, , b] <- t(apply(res[[1]][, , b], 1, function(x)
        gam(x ~ s(
          argvals, bs = "tp", k = (nknots + 1)
        ), method = "GCV.Cp")$fitted.values))
    }

    boot_df =
      betaHat_boot[2,,] %>%
      as_tibble() %>%
      mutate(l = row_number()) %>%
      pivot_longer(cols = starts_with("V"), names_to = "rep", values_to = "value")  %>%
      mutate(id = sub(".*V", "", rep) %>% as.numeric) %>%
      filter(id <= 10)
    p2 =
      boot_df %>%
      ggplot(aes(x = l, y = value, color = factor(rep))) +
      theme(legend.position = "none") +
      geom_line() +
      # Add labeled lines for custom legend
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      geom_line(data = btilde_df,
                aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      # Manual legend
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Raw Estimate" = "blue")) +
      labs(x = "Functional Domain", y = "Coefficient", color = "Line Type",
           title = "Smoothed Beta Coef Estimates for 10 Repetitions") +
      theme_light()


    betaHat.var <- array(NA, dim = c(L, L, 2))

    for (r in 1:2) {
      betaHat.var[, , r] <- mult_fac * var(t(betaHat_boot[r, , ]))
    }
    var_x = betaHat.var[,,2]

    bhat_df2 = betaHat[[1]][2,] %>%
      as_tibble() %>%
      rename(beta = value) %>%
      mutate(l = row_number())

    N <- 10000
    qn <- rep(0, length = nrow(betaHat[[1]]))
    set.seed(456)
    i = 2

    est_bs <- t(betaHat_boot[i, , ])


    fit_fpca <- suppressWarnings(refund::fpca.face(est_bs, knots = nknots_fpca)) # suppress sqrt(Eigen$values) NaNs
    ## extract estimated eigenfunctions/eigenvalues
    phi <- fit_fpca$efunctions
    lambda <- fit_fpca$evalues
    K <- length(fit_fpca$evalues)

    ## simulate random coefficients
    theta <- matrix(stats::rnorm(N * K), nrow = N, ncol = K) # generate independent standard normals
    if (K == 1) {
      theta <- theta * sqrt(lambda) # scale to have appropriate variance
      X_new <- tcrossprod(theta, phi) # simulate new functions
    } else{
      theta <- theta %*% diag(sqrt(lambda)) # scale to have appropriate variance
      X_new <- tcrossprod(theta, phi) # simulate new functions
    }
    x_sample <- X_new + t(fit_fpca$mu %o% rep(1, N)) # add back in the mean function
    Sigma_sd <- Rfast::colVars(x_sample, std = TRUE, na.rm = FALSE) # standard deviation: apply(x_sample, 2, sd)
    x_mean <- colMeans(est_bs)
    un <- rep(NA, N)
    for (j in 1:N) {
      un[j] <- max(abs((x_sample[j, ] - x_mean) / Sigma_sd))
    }
    qn[i] <- stats::quantile(un, 0.95)
    p3 =
      un %>%
      as_tibble() %>%
      ggplot(aes(x = value)) +
      geom_histogram(color = "black") +
      geom_vline(aes(xintercept = stats::quantile(un, 0.95)), color = "red") +
      labs(x = "Value", y= "Count", title = "Histogram of Max Statisics") +
      theme_light()


    cma = qn[2]
    var_x_df =
      diag(var_x) %>%
      as_tibble() %>%
      mutate(l = row_number()) %>%
      left_join(bhat_df2) %>%
      mutate(min = beta - 1.96 * sqrt(value),
             max = beta + 1.96 * sqrt(value),
             cma_min  = beta - cma * sqrt(value),
             cma_max = beta + cma * sqrt(value))


    p4 = var_x_df %>%
      ggplot() +
      geom_ribbon(aes(x = l, ymin = cma_min, ymax = cma_max), fill = "darkgrey") +
      geom_ribbon(aes(x = l, ymin = min, ymax = max), fill = "lightgrey") +
      theme_light() +
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      # geom_line(data = btilde_df %>% filter(name == "X"),
      #           aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Raw Estimate" = "blue"))


    # get some coverage stats
    margin_pw = 1.96 * sqrt(diag(betaHat.var[,,2]))
    margin_joint = cma * sqrt(diag(betaHat.var[,,2]))
    upper_j = betaHat[[1]][2,] + margin_joint
    lower_j = betaHat[[1]][2,] - margin_joint
    upper = betaHat[[1]][2,] + margin_pw
    lower = betaHat[[1]][2,] - margin_pw
    p = 2
    true = beta_true[p,]
    covered <- (lower_j < true) & (upper_j > true)
    covered_pw <- (lower < true) & (upper > true)
    cover_joint =  all(covered)  # Joint coverage check
    cover_pw = sum(covered_pw)


    p4 = var_x_df %>%
      ggplot() +
      geom_ribbon(aes(x = l, ymin = cma_min, ymax = cma_max), fill = "darkgrey") +
      geom_ribbon(aes(x = l, ymin = min, ymax = max), fill = "lightgrey") +
      theme_light() +
      geom_line(data = btrue_df,
                aes(x = l, y = value, color = "Truth"), linewidth = 1) +
      geom_line(data = bhat_df,
                aes(x = l, y = value, color = "Smoothed Estimate"), linewidth = 1) +
      # geom_line(data = btilde_df %>% filter(name == "X"),
      #           aes(x = l, y = value, color = "Raw Estimate"), linewidth = 1) +
      scale_color_manual(values = c("Truth" = "red",
                                    "Smoothed Estimate" = "black",
                                    "Raw Estimate" = "blue"))  +
      labs(title = "Confidence Intervals",color = "Line Type", x = "Functional Domain", y = "Coefficient",
           subtitle = paste0("Joint covered = ", cover_joint, " PW covered = ", cover_pw))

    print((p1 + p2 + plot_layout(guides = "collect", axes = "collect"))/ (p3 + p4))
  })
  rm(x)
}

dev.off()

