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
fname = paste0("sim_fold_", sprintf("%03d", ifold), "_boot.rds")



if(!file.exists(here::here("results", "simulations", "non_survey", fname)) ||
   force) {
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
  #   ncores_outer = ncores_total
  # } else{
    ncores_outer = ceiling(ncores_total / 2)
    ncores_inner = ncores_total - ncores_outer

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
        b_hat =
          sample_data %>%
          get_betatilde(., family = family) %>%
          get_betahat(., L = len)

        boots = run_boots(
          b_hat,
          set_seed = TRUE,
          data = sample_data,
          seed = iter,
          family = family,
          num_boots = B,
          L = len,
          parallel = TRUE,
          ncores_inner = ncores_inner
        )

        ci = get_cis_boot(boots, b_hat, smooth_for_ci = TRUE, L = len)
        elapsed_time = toc(quiet = TRUE)
        stats_fui = get_coverage_stats_fui(ci, beta_true = data$beta_true, L = len) %>%
          mutate(time = unname(elapsed_time$toc - elapsed_time$tic))
        rm(boots)
      # } else {
      #   tic()
      #   beta_tilde =
      #     sample_data %>%
      #     get_betatilde(., family = family)
      #
      #
      #   b_hat = beta_tilde %>%
      #     get_betahat(., L = len)
      #
      #   ci = get_cis_analytic(betaHat = b_hat, betaTilde = beta_tilde, L = len, data = sample_data)
      #   elapsed_time = toc(quiet = TRUE)
      #   stats_fui = get_coverage_stats_fui(ci, beta_true = data$beta_true, L = len) %>%
      #     mutate(time = unname(elapsed_time$toc - elapsed_time$tic))
      # }

      elapsed_time <- NULL

      #### famm
      Y_mat =
        sample_data %>%
        select(starts_with("Y")) %>%
        as.matrix()

      tic()
      fit_pffr =
        sample_data %>%
        select(-starts_with("Y")) %>%
        mutate(Y_mat = I(Y_mat)) %>%
        pffr(
          Y_mat ~ X,
          data = .,
          algorithm = "bam",
          family = family,
          method = "fREML",
          bs.yindex = list(
            bs = "ps",
            k = 15,
            m = c(2, 1)
          )
        )
      elapsed_time = toc(quiet = TRUE)
      stats_famm = get_coverage_stats_famm(fit_pffr, beta_true = data$beta_true, L = len) %>%
        mutate(time = unname(elapsed_time$toc - elapsed_time$tic))

      ## bayes
      if(family == "gaussian" && run_bayes){ # can only run bayes for gaussian models
        Y_mat =
          sample_data %>%
          select(starts_with("Y")) %>%
          as.matrix()
        df = sample_data %>%
          select(-starts_with("Y")) %>%
          mutate(Ym = Y_mat)
        tic()
        fit_bayes = refund::bayes_fosr(Ym ~ X, est.method = "Gibbs", data = df)
        elapsed_time = toc(quiet = TRUE)
        rm(df)

        stats_bayes = get_coverage_stats_bayes(fit_bayes, beta_true = data$beta_true, L = len) %>%
          mutate(time = unname(elapsed_time$toc - elapsed_time$tic))

        ret =
          stats_famm %>%
          mutate(type = "FAMM") %>%
          bind_rows(stats_fui %>% mutate(type = "FUI")) %>%
          bind_rows(stats_bayes %>% mutate(type = "Bayes"))
        rm(stats_famm, stats_fui, elapsed_time, Y_mat, sample_data, b_hat, ci, data, fit_pffr, fit_bayes, stats_bayes)
      } else {
        ret =
          stats_famm %>%
          mutate(type = "FAMM") %>%
          bind_rows(stats_fui %>% mutate(type = "FUI"))
        rm(stats_famm, stats_fui, elapsed_time, Y_mat, sample_data, b_hat, ci, data, fit_pffr)
      }
      # print(iter)
      return(ret)
    })
    x
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE)

  print(sim_res)

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
}
plan(sequential)
