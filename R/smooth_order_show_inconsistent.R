library(ggplot2)
library(dplyr)
library(mgcv)
library(tidyr)
library(purrr)
library(future)
library(furrr)
library(readr)
force = FALSE
theme_set(theme_light(base_size = 14))

# data generating funciton, same as our simulation
generate_data = function(I = 1000, seed = 2213, L = 50, sd = 2, snr_eps = 1){
  set.seed(seed)
  X_des  = cbind(1, rnorm(I, 0, 2))

  ## simulate true beta based on scenarios
  grid  = seq(0, 1, length = L)
  beta_fixed  = matrix(NA, 2, L)
  beta_fixed[1, ]  = -0.15 - 0.1 * sin(2 * pi * grid) - 0.1 * cos(2 * pi * grid)
  beta_fixed[2, ]  = dnorm(grid, 0.6, 0.15) / 20

  lp = X_des %*% beta_fixed
  sd_lp = sd(as.vector(lp))
  sigma = sd_lp / snr_eps
  set.seed(seed)
  # Y = lp + matrix(rnorm(I * L, sd = sd), I, L)
  Y_obs = matrix(
    rnorm(n = I * L,
          mean = as.vector(t(lp)),
          sd = sigma), # need to use t to put in correct order
    nrow = I,
    ncol = L,
    byrow = TRUE
  )
  return(list(Y = Y_obs,
              X = X_des,
              bt = beta_fixed))
}

# function to calculate bias from smooth first and smooth second processes
calculate_bias = function(I, L = 50, nknots = min(50/2, 35), seed) {
  data = generate_data(I = I, L = L, seed = seed)
  argvals = 1:L
  X = data$X
  Y = data$Y
  bt = data$bt
  # smooth second
  beta_tilde = solve(t(X) %*%X) %*% t(X) %*% Y
  beta_hat_ss = t(apply(beta_tilde, 1, function(x) {
    gam(x ~ s(argvals, bs = "tp", k = nknots + 1),
        method = "GCV.Cp")$fitted.values
  }))

  # smooth first
  smooth_Y = t(apply(Y, 1, function(x) {
    gam(x ~ s(argvals, bs = "tp", k = nknots + 1),
        method = "GCV.Cp")$fitted.values
  }))
  beta_hat_sf = solve(t(X) %*% X) %*% t(X) %*% smooth_Y

  return(list(beta_hat_ss = beta_hat_ss,
              beta_hat_sf = beta_hat_sf,
              bt = bt))
}

# function to calculate bias when we use the right EDF
calculate_bias_fixed_edf = function(I, L = 50, nknots = min(50/2, 35), seed, target_edf = NULL) {
  data = generate_data(I = I, L = L, seed = seed)
  argvals = 1:L
  X = data$X
  Y = data$Y
  bt = data$bt


  # smooth second (unchanged - use GCV)
  beta_tilde = solve(t(X) %*% X) %*% t(X) %*% Y
  beta_hat_ss = t(apply(beta_tilde, 1, function(x) {
    gam(x ~ s(argvals, bs = "tp", k = nknots + 1),
        method = "GCV.Cp")$fitted.values
  }))

  # Get the EDF that smooth-second selects (for reference)
  ss_edf = apply(beta_tilde, 1, function(x) {
    sum(gam(x ~ s(argvals, bs = "tp", k = nknots + 1),
            method = "GCV.Cp")$edf)
  })

  # smooth first with GCV (original - wrong)
  smooth_Y_gcv = t(apply(Y, 1, function(x) {
    gam(x ~ s(argvals, bs = "tp", k = nknots + 1),
        method = "GCV.Cp")$fitted.values
  }))
  beta_hat_sf_gcv = solve(t(X) %*% X) %*% t(X) %*% smooth_Y_gcv

  # Get per-subject EDF from GCV smooth-first
  sf_edf_gcv = apply(Y, 1, function(x) {
    sum(gam(x ~ s(argvals, bs = "tp", k = nknots + 1),
            method = "GCV.Cp")$edf)
  })

  # smooth first with FIXED EDF (correct - use target_edf or ss_edf)
  if(is.null(target_edf)) {
    target_edf = mean(ss_edf)  # use the EDF from smooth-second
  }

  # To fix EDF, we need to find the smoothing parameter that gives target_edf
  # One approach: fit with GCV, then refit with fixed sp scaled to achieve target edf
  smooth_Y_fixed = t(apply(Y, 1, function(x) {
    # First fit to get a starting point
    fit_gcv = gam(x ~ s(argvals, bs = "tp", k = nknots + 1), method = "GCV.Cp")

    # Binary search for sp that gives target_edf
    find_sp_for_edf = function(target) {
      sp_low = 1e-10
      sp_high = 1e6

      for(i in 1:50) {
        sp_mid = sqrt(sp_low * sp_high)
        fit_tmp = gam(x ~ s(argvals, bs = "tp", k = nknots + 1),
                      sp = sp_mid, method = "GCV.Cp")
        current_edf = sum(fit_tmp$edf)

        if(abs(current_edf - target) < 0.1) break
        if(current_edf > target) {
          sp_low = sp_mid
        } else {
          sp_high = sp_mid
        }
      }
      return(sp_mid)
    }

    sp_fixed = find_sp_for_edf(target_edf)
    fit_fixed = gam(x ~ s(argvals, bs = "tp", k = nknots + 1), sp = sp_fixed)
    return(fit_fixed$fitted.values)
  }))

  beta_hat_sf_fixed = solve(t(X) %*% X) %*% t(X) %*% smooth_Y_fixed

  return(list(
    beta_hat_ss = beta_hat_ss,
    beta_hat_sf_gcv = beta_hat_sf_gcv,
    beta_hat_sf_fixed = beta_hat_sf_fixed,
    bt = bt,
    ss_edf = mean(ss_edf),
    sf_edf_gcv = mean(sf_edf_gcv)
  ))
}

run_sim_fixed_edf = function(n = 500, n_iter = 50, parallel = TRUE, L = 50){
  if(parallel){
    plan(multisession, workers = parallelly::availableCores() - 1)
    res = future_map(1:n_iter,
                     ~ calculate_bias_fixed_edf(I = n, L = L, nknots = min(L/2, 35), seed = .x),
                     .progress = TRUE)
    plan(sequential)
  } else {
    res = map(1:n_iter,
              ~ calculate_bias_fixed_edf(I = n, L = L, nknots = min(L/2, 35), seed = .x))
  }

  bt = res[[1]]$bt

  # Stack arrays
  ss_array = array(sapply(res, function(x) x$beta_hat_ss), dim = c(2, L, n_iter))
  sf_array = array(sapply(res, function(x) x$beta_hat_sf_gcv), dim = c(2, L, n_iter))
  sff_array = array(sapply(res, function(x) x$beta_hat_sf_fixed), dim = c(2, L, n_iter))

  # Means
  ss_mean = apply(ss_array, c(1, 2), mean)
  sf_gcv_mean = apply(sf_array, c(1, 2), mean)
  sf_fixed_mean = apply(sff_array, c(1, 2), mean)

  # Bias
  ss_bias = ss_mean - bt
  sf_gcv_bias = sf_gcv_mean - bt
  sf_fixed_bias = sf_fixed_mean - bt

  ss_coverage = matrix(NA, 2, L)
  sff_coverage = matrix(NA, 2, L)
  sf_coverage = matrix(NA, 2, L)
  for(j in 1:2){
    for(l in 1:L){
      ss_coverage[j, l] = mean(abs(ss_array[j, l, ] - bt[j, l]) < 1.96 * ss_sd[j, l])
      sf_coverage[j, l] = mean(abs(sf_array[j, l, ] - bt[j, l]) < 1.96 * sf_sd[j, l])
      sff_coverage[j, l] = mean(abs(sff_array[j, l, ] - bt[j, l]) < 1.96 * sf_sd[j, l])

    }
  }


  # EDF diagnostics
  mean_ss_edf = mean(sapply(res, function(x) x$ss_edf))
  mean_sf_edf_gcv = mean(sapply(res, function(x) x$sf_edf_gcv))

  return(tibble(
    s = 1:L,
    ss_bias_b0 = ss_bias[1,],
    ss_bias_b1 = ss_bias[2,],
    sf_bias_b0 = sf_gcv_bias[1,],
    sf_bias_b1 = sf_gcv_bias[2,],
    sff_bias_b0 = sf_fixed_bias[1,],
    sff_bias_b1 = sf_fixed_bias[2,],
    ss_edf = mean_ss_edf,
    sff_edf = mean_sf_edf_gcv,
    ss_cover_b0 = ss_coverage[1,],
    ss_cover_b1 = ss_coverage[2,],
    sf_cover_b0 = sf_coverage[1,],
    sf_cover_b1 = sf_coverage[2,],
    sff_cover_b0 = sff_coverage[1,],
    sff_cover_b1 = sff_coverage[2,]
  ))
}
# function to run the simulation for a given sample size
run_sim = function(n = 500, n_iter = 50, parallel = TRUE, L = 50){
  if(parallel){
    plan(multisession, workers = parallelly::availableCores() - 1)
    res = future_map(.x = 1:n_iter, .f = calculate_bias, I = n, L = 50, nknots = min(L/2, 35),
                     .progress = TRUE)
    plan(sequential)
  } else{
    res = map(.x = 1:n_iter, .f = calculate_bias, I = n, L = 50, nknots = min(L/2, 35))
  }

  bt = res[[1]]$bt
  ss_array = array(sapply(res, function(x) x$beta_hat_ss),
                   dim = c(2, L, n_iter))

  sf_array = array(sapply(res, function(x) x$beta_hat_sf),
                   dim = c(2, L, n_iter))

  ss_mean = apply(ss_array, c(1, 2), mean)
  sf_mean = apply(sf_array, c(1, 2), mean)
  ss_sd = apply(ss_array, c(1, 2), sd)
  sf_sd = apply(sf_array, c(1, 2), sd)

  ss_bias = ss_mean - bt
  sf_bias = sf_mean - bt

  # get coverage
  ss_coverage = matrix(NA, 2, L)
  sf_coverage = matrix(NA, 2, L)

  for(j in 1:2){
    for(l in 1:L){
      ss_coverage[j, l] = mean(abs(ss_array[j, l, ] - bt[j, l]) < 1.96 * ss_sd[j, l])
      sf_coverage[j, l] = mean(abs(sf_array[j, l, ] - bt[j, l]) < 1.96 * sf_sd[j, l])
    }
  }



  return(tibble(s = 1:L,
                ss_bias_b0 = ss_bias[1,],
                ss_bias_b1 = ss_bias[2,],
                sf_bias_b0 = sf_bias[1,],
                sf_bias_b1 = sf_bias[2,],
                ss_mean_b0 = ss_mean[1,],
                ss_mean_b1 = ss_mean[2,],
                sf_mean_b0 = sf_mean[1,],
                sf_mean_b1 = sf_mean[2,],
                ss_cover_b0 = ss_coverage[1,],
                ss_cover_b1 = ss_coverage[2,],
                sf_cover_b0 = sf_coverage[1,],
                sf_cover_b1 = sf_coverage[2,]
                ))
}

if (!file.exists(here::here("results", "simulations", "smooth_order_consistency.rds")) || force){
  n_result = map(.x = c(500, 1000, 5000, 10000),
                 .f = run_sim,
                 n_iter = 200,
                 parallel = TRUE)

  n_result =
    n_result |> bind_rows(.id = "n")

  n_result =
    n_result |>
    mutate(n = case_when(n == "1" ~ 500,
                         n == "2" ~ 1000,
                         n == "3" ~ 5000,
                         .default = 10000))

  write_rds(n_result, here::here("results", "simulations", "smooth_order_consistency.rds"))

  n_result = map(.x = c(500, 1000, 5000, 10000),
                 .f = run_sim_fixed_edf,
                 n_iter = 200,
                 parallel = TRUE)

} else n_result = read_rds(here::here("results", "simulations", "smooth_order_consistency.rds"))


######
# make plots
######
L = 50
beta_fixed  = matrix(NA, 2, L)
beta_fixed[1, ]  = -0.15 - 0.1 * sin(2 * pi * grid) - 0.1 * cos(2 * pi * grid)
beta_fixed[2, ]  = dnorm(grid, 0.6, 0.15) / 20

bt_df =
  beta_fixed |>
  as_tibble() |>
  mutate(var = c("b0", "b1")) |>
  pivot_longer(cols = -var, names_to = "s", names_transform = ~as.numeric(sub(".*V", "", .x)))

n_result |>
  pivot_longer(cols = contains("b")) |>
  filter(grepl("bias", name)) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1") |>
  mutate(smooth_order = factor(smooth_order, labels = c("Smooth first", "Smooth second"))) |>
  ggplot(aes(x = s, y = value, color = factor(n))) +
  geom_line() +
  facet_wrap(.~smooth_order) +
  labs(x = "Functional Domain", y = "Bias", title = bquote("Bias for"~beta[1])) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "darkgrey") +
  # scale_color_viridis_d(name = "Sample size", option = "B") +
  scale_color_manual(values = c("#FF7F00FF", "#19B2FFFF", "#654CFFFF", "#E51932FF"), name = "Sample Size") +
  theme(legend.position = c(0.9, 0.2))

n_result |>
  pivot_longer(cols = contains("b")) |>
  filter(grepl("bias", name)) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1") |>
  mutate(smooth_order = factor(smooth_order, labels = c("Smooth first", "Smooth second"))) |>
  ggplot(aes(x = s, y = value, color = factor(n))) +
  geom_line() +
  facet_wrap(.~smooth_order, scales = "free_y") +
  labs(x = "Functional Domain", y = "Bias", title = bquote("Bias for"~beta[1])) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "darkgrey") +
  # scale_color_viridis_d(name = "Sample size", option = "B") +
  scale_color_manual(values = c("#FF7F00FF", "#19B2FFFF", "#654CFFFF", "#E51932FF"), name = "Sample Size") +
  theme(legend.position = c(0.9, 0.2))


n_result |>
  pivot_longer(cols = contains("b")) |>
  filter(grepl("bias", name)) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0") |>
  mutate(smooth_order = factor(smooth_order, labels = c("Smooth first", "Smooth second"))) |>
  ggplot(aes(x = s, y = value, color = factor(n))) +
  geom_line() +
  facet_wrap(.~smooth_order) +
  labs(x = "Functional Domain", y = "Bias", title = bquote("Bias for"~beta[0])) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "darkgrey") +
  # scale_color_viridis_d(name = "Sample size", option = "B") +
  scale_color_manual(values = c("#FF7F00FF", "#19B2FFFF", "#654CFFFF", "#E51932FF"), name = "Sample Size") +
  theme(legend.position = c(0.9, 0.2))

n_result |>
  pivot_longer(cols = contains("b")) |>
  filter(grepl("bias", name)) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0") |>
  mutate(smooth_order = factor(smooth_order, labels = c("Smooth first", "Smooth second"))) |>
  ggplot(aes(x = s, y = value, color = factor(n))) +
  geom_line() +
  facet_wrap(.~smooth_order, scales = "free_y") +
  labs(x = "Functional Domain", y = "Bias", title = bquote("Bias for"~beta[0])) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "darkgrey") +
  # scale_color_viridis_d(name = "Sample size", option = "B") +
  scale_color_manual(values = c("#FF7F00FF", "#19B2FFFF", "#654CFFFF", "#E51932FF"), name = "Sample Size") +
  theme(legend.position = c(0.9, 0.2))

n_result |>
  filter(n == 1000) |>
  select(s, contains("mean")) |>
  pivot_longer(cols = -s)  |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0") |>
  ggplot(aes(x = s, y = value, color = smooth_order)) +
  geom_line(linewidth = .8) +
  geom_line(data = bt_df |> filter(var == "b0") |> mutate(smooth_order = "truth"), linewidth = .8) +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF", "darkgrey"), name = "",
                     labels = c("Smooth First", "Smooth Second", "Truth")) +
  labs(x = "Functional Domain", y = "Value", title = bquote("True and estimated"~beta[0])) +
  theme(legend.position = c(0.8, 0.2),
        legend.title = element_blank())

n_result |>
  filter(n == 1000) |>
  select(s, contains("mean")) |>
  pivot_longer(cols = -s)  |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1") |>
  ggplot(aes(x = s, y = value, color = smooth_order)) +
  geom_line(linewidth = .8) +
  geom_line(data = bt_df |> filter(var == "b1") |> mutate(smooth_order = "truth"), linewidth = .8) +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF", "darkgrey"), name = "",
                     labels = c("Smooth First", "Smooth Second", "Truth")) +
  labs(x = "Functional Domain", y = "Value", title = bquote("True and estimated"~beta[1])) +
  theme(legend.position = c(0.6, 0.2),
        legend.title = element_blank())


n_result |>
  filter(n == 1000) |>
  select(s, contains("cover")) |>
  pivot_longer(cols = -s)  |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0") |>
  ggplot(aes(x = s, y = value, color = smooth_order)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF"), name = "",
                     labels = c("Smooth First", "Smooth Second")) +
  labs(x = "Functional Domain", y = "Pointwise coverage", title = bquote("Empirical coverage for"~beta[0]~"at n = 1000")) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  geom_hline(aes(yintercept = 0.95), color = "darkgrey", linetype = "dashed") +
  annotate(geom = "text", label = "Nominal coverage rate = 0.95", x = 5, y = .95, color = "darkgrey", size = 3, vjust = -1)

n_result |>
  filter(n == 1000) |>
  select(s, contains("cover")) |>
  pivot_longer(cols = -s)  |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1") |>
  ggplot(aes(x = s, y = value, color = smooth_order)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF"), name = "",
                     labels = c("Smooth First", "Smooth Second")) +
  labs(x = "Functional Domain", y = "Pointwise coverage", title = bquote("Empirical coverage for"~beta[1]~"at n = 1000")) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  geom_hline(aes(yintercept = 0.95), color = "darkgrey", linetype = "dashed") +
  annotate(geom = "text", label = "Nominal coverage rate = 0.95", x = 5, y = .95, color = "darkgrey", size = 3, vjust = -1)

# First, get the bias data and determine scaling
bias_data = n_result |>
  filter(n == 1000) |>
  select(s, contains("bias")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1")

# Scaling parameters: map bias range to ~[0, 0.5] so it doesn't overlap too much with coverage
bias_range = range(bias_data$value)
scale_factor = 0.4 / max(abs(bias_range))  # scale to fit in lower part of plot
offset = 0.5  # center the bias around 0.5 on the primary axis

# Coverage data
cover_data = n_result |>
  filter(n == 1000) |>
  select(s, contains("cover")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1")

# Plot
ggplot() +
  # Coverage lines (primary axis)
  geom_line(data = cover_data,
            aes(x = s, y = value, color = smooth_order),
            linewidth = 1) +
  # Bias lines (scaled to primary axis)
  geom_line(data = bias_data,
            aes(x = s, y = value * scale_factor + offset, linetype = smooth_order),
            color = "grey40", linewidth = 0.8) +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF"), name = "Coverage",
                     labels = c("Smooth First", "Smooth Second")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Bias",
                        labels = c("Smooth First", "Smooth Second")) +
  scale_y_continuous(
    name = "Pointwise coverage",
    limits = c(0, 1),
    sec.axis = sec_axis(~ (. - offset) / scale_factor, name = "Bias")
  ) +
  labs(x = "Functional Domain",
       title = bquote("Empirical coverage for" ~ beta[1] ~ "at n = 1000")) +
  geom_hline(yintercept = 0.95, color = "darkgrey", linetype = "dashed") +
  geom_hline(yintercept = offset, color = "grey40", linetype = "dotted", alpha = 0.5) +  # zero line for bias
  annotate(geom = "text", label = "Nominal = 0.95", x = 5, y = 0.95,
           color = "darkgrey", size = 3, vjust = -1) +
  theme(legend.position = "bottom")

bias_data = n_result |>
  filter(n == 1000) |>
  select(s, contains("bias")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0")

# Scaling parameters: map bias range to ~[0, 0.5] so it doesn't overlap too much with coverage
bias_range = range(bias_data$value)
scale_factor = 0.4 / max(abs(bias_range))  # scale to fit in lower part of plot
offset = 0.5  # center the bias around 0.5 on the primary axis

# Coverage data
cover_data = n_result |>
  filter(n == 1000) |>
  select(s, contains("cover")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0")

# Plot
ggplot() +
  # Coverage lines (primary axis)
  geom_line(data = cover_data,
            aes(x = s, y = value, color = smooth_order),
            linewidth = 1) +
  # Bias lines (scaled to primary axis)
  geom_line(data = bias_data,
            aes(x = s, y = value * scale_factor + offset, linetype = smooth_order),
            color = "grey40", linewidth = 0.8) +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF"), name = "Coverage",
                     labels = c("Smooth First", "Smooth Second")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Bias",
                        labels = c("Smooth First", "Smooth Second")) +
  scale_y_continuous(
    name = "Pointwise coverage",
    limits = c(0, 1),
    sec.axis = sec_axis(~ (. - offset) / scale_factor, name = "Bias")
  ) +
  labs(x = "Functional Domain",
       title = bquote("Empirical coverage for" ~ beta[0] ~ "at n = 1000")) +
  geom_hline(yintercept = 0.95, color = "darkgrey", linetype = "dashed") +
  geom_hline(yintercept = offset, color = "grey40", linetype = "dotted", alpha = 0.5) +  # zero line for bias
  annotate(geom = "text", label = "Nominal = 0.95", x = 5, y = 0.95,
           color = "darkgrey", size = 3, vjust = -1) +
  theme(legend.position = "bottom")
