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

theme_set(theme_light())
nsim = 15
B = 500
force = TRUE

get_plot_df = function(sim_res, tname){
  num_var = nrow(sim_res$betaHat)
  map_dfr(.x = 1:num_var,
          .f = function(r){
            data.frame(s = 1:length(sim_res$betaHat[r,]),
                       beta = sim_res$betaHat[r,],
                       lower = sim_res$betaHat[r, ] - 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
                       upper = sim_res$betaHat[r, ] + 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
                       lower.joint = sim_res$betaHat[r, ] - sim_res$qn[r] *
                         sqrt(diag(sim_res$betaHat.var[, , r])),
                       upper.joint = sim_res$betaHat[r, ] + sim_res$qn[r] * sqrt(diag(sim_res$betaHat.var[, , r]))) %>%
              mutate(name = sim_res$betaHat %>% rownames() %>% .[r])
          }) %>%
    mutate(varname = tname)
}

fit_types = c('weighted', 'unweighted')
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'weighted', 'unweighted')

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
mult_fac = 1.2
nknots <- min(round(L / 2), nknots_min)
nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
nknots_fpca <- min(round(L / 2), 35)
argvals = 1:L
iter = 1
pdf(here::here("wtd_debug.pdf"))
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

    data$weight = data$weight / mean(data$weight)

    beta_true = lst$beta_true
    betaTilde = map(
      .x = fit_types,
      .f = get_betatilde,
      data = data,
      family = temp$family
    )

    betaHat = map(.x = betaTilde, get_betahat, L = temp$len)


    plan(multisession, workers = 8)
    res = future_map(
      .x = boot_types,
      .f = run_boots_fast_internal,
      betaHat = betaHat[[1]],
      data = data,
      family = temp$family,
      num_boots = B,
      seed = 2025,
      L = temp$len,
      samp_stages = c("PPSWOR", "Poisson"),
      .options = furrr_options(seed = TRUE),
      ncores = 2
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

    plan(sequential)



    # res_x = res[[1]][2, , ]

    btrue_df =
      beta_true %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    bhat_df_wt =
      betaHat[[1]] %>%
      as_tibble() %>%
      mutate(beta = c("beta_0", "beta_1")) %>%
      pivot_longer(cols = -beta) %>%
      mutate(x = as.numeric(sub(".*V", "", name))) %>%
      filter(beta == "beta_1") %>%
      select(-beta) %>%
      rename(l = x)

    bhat_df_unwt =
      betaHat[[2]] %>%
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

    p1 = btrue_df %>%
      ggplot(aes(x = l, y = value, color = "Truth")) +
      geom_line() +
      geom_line(data = bhat_df_unwt, aes(x = l, y = value, color = "Unweighted")) +
      geom_line(data = bhat_df_wt, aes(x = l, y = value, color = "Weighted")) +
      scale_color_manual(values = c("Truth" = "black",
                                    "Unweighted" = "red",
                                    "Weighted" = "blue")) +
                           labs(y= "beta") +
      scale_y_continuous(limits = c(-.1, .2))

    plot_dfs = map2(.x = cis, .f = get_plot_df,
                   .y = boot_types)

    p2 = plot_dfs[[1]] %>%
      filter(name == "X") %>%
      ggplot() +
      geom_ribbon(aes(x  = s, ymin = lower.joint, ymax = upper.joint), fill = "darkgrey") +
      geom_ribbon(aes(x  = s, ymin = lower, ymax = upper), fill = "lightgrey") +
      geom_line(aes(x = s, y = beta, color = "Estimate")) +
      geom_line(data = btrue_df, aes(x=l, y = value, color = "Truth")) +
      labs(title = "BRR", x = "l") +
      scale_color_manual(values = c("Truth" = "black",
                                    "Estimate" = "blue")) +

      scale_y_continuous(limits = c(-.1, .2))

    p3 = plot_dfs[[2]] %>%
      filter(name == "X") %>%
      ggplot() +
      geom_ribbon(aes(x  = s, ymin = lower.joint, ymax = upper.joint), fill = "darkgrey") +
      geom_ribbon(aes(x  = s, ymin = lower, ymax = upper), fill = "lightgrey") +
      geom_line(aes(x = s, y = beta, color = "Estimate")) +
      geom_line(data = btrue_df, aes(x=l, y = value, color = "Truth")) +
      labs(title = "RWYB",
           x = "l") +
      scale_color_manual(values = c("Truth" = "black",
                                    "Estimate" = "blue")) +

      scale_y_continuous(limits = c(-.1, .2))

    p4 = plot_dfs[[3]] %>%
      filter(name == "X") %>%
      ggplot() +
      geom_ribbon(aes(x  = s, ymin = lower.joint, ymax = upper.joint), fill = "darkgrey") +
      geom_ribbon(aes(x  = s, ymin = lower, ymax = upper), fill = "lightgrey") +
      geom_line(aes(x = s, y = beta, color = "Estimate")) +
      geom_line(data = btrue_df, aes(x=l, y = value, color = "Truth")) +
      labs(title = "Weighted",
           x = "l") +
      scale_color_manual(values = c("Truth" = "black",
                                    "Estimate" = "blue")) +

      scale_y_continuous(limits = c(-.1, .2))

    p5 = plot_dfs[[4]] %>%
      filter(name == "X") %>%
      ggplot() +
      geom_ribbon(aes(x  = s, ymin = lower.joint, ymax = upper.joint), fill = "darkgrey") +
      geom_ribbon(aes(x  = s, ymin = lower, ymax = upper), fill = "lightgrey") +
      geom_line(aes(x = s, y = beta, color = "Estimate")) +
      geom_line(data = btrue_df, aes(x=l, y = value, color = "Truth")) +
      labs(title = "Unweighted",
           x = "l") +
      scale_color_manual(values = c("Truth" = "black",
                                    "Estimate" = "red")) +
      scale_y_continuous(limits = c(-.1, .2))


    print((p1 | ((p2 + p3) / (p4 + p5)) + plot_layout(guides = "collect", axes = "collect")))
  })
  rm(x)
}
dev.off()
