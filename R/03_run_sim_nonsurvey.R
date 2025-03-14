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
source(here::here("R", "sim_functions_unwt.R"))
ncores = parallelly::availableCores() - 1

force = FALSE

source(here::here("R", "utils.R"))
sim_settings = read_rds(here::here("sim_data", "sim_settings_nonsurvey.rds"))
ifold = get_fold()
ifold = 2
fname = paste0("sim_fold_", sprintf("%03d", ifold), ".rds")


# B = 500
# scen = 1
# fam = "gaussian"
# sigma = 1
# samp_size = 500

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
  plan(multisession, workers = ncores)
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
      tic()
      b_hat =
        sample_data %>%
        get_betatilde() %>%
        get_betahat(., L = len)

      boots = run_boots(
        b_hat,
        set_seed = TRUE,
        data = sample_data,
        seed = iter,
        family = family,
        num_boots = B,
        L = len
      )

      ci = get_cis(boots, b_hat, smooth_for_ci = TRUE, L = len)
      elapsed_time = toc(quiet = TRUE)
      stats_fui = get_coverage_stats_fui(ci, beta_true = data$beta_true, L = len) %>%
        mutate(time = unname(elapsed_time$toc - elapsed_time$tic))

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
          Y_mat ~ X + s(ID, bs = "re"),
          data = .,
          algorithm = "bam",
          discrete = TRUE,
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

      ret =
        stats_famm %>%
        mutate(type = "FAMM") %>%
        bind_rows(stats_fui %>% mutate(type = "FUI"))
      rm(stats_famm, stats_fui, elapsed_time, Y_mat, boots, sample_data, b_hat, ci, data, fit_pffr)
      return(ret)
    })
    x
  },
  .options = furrr_options(seed = TRUE))

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
  #   pivot_longer(cols = c(MISE, mean_pw_se, cover_pw)) %>%
  #   ggplot(aes(x = var, y = value, color = type)) +
  #   geom_boxplot(outlier.shape = NA)+
  #   facet_wrap(.~name, scales = "free_y")

  if (!dir.exists(here::here("results", "simulations", "non_survey"))) {
    dir.create(
      here::here("results", "simulations", "non_survey", recursive = TRUE))
    }

  write_rds(result,
            here::here("results", "simulations", "non_survey", fname))
}
