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
source(here::here("R_cp", "sim_functions_unwt.R"))
source(here::here("R_cp", "01_sim_functions.R"))
source(here::here("R_cp", "00_data_gen_function_ff.R"))
source(here::here("R_cp", "create_nonsurvey_settings.R"))
source(here::here("R_cp", "utils.R"))


force = TRUE
force_iter = FALSE
nsim = 200
B = 500

ifold = get_fold()
# ifold = 1
outfile = here::here("results", "simulations", "non_survey", paste0("fold_", sprintf("%03d", ifold), ".rds"))

if(!file.exists(outfile) || force) { # if the file doesn't exist or we want to force re-do it
  if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)
  partial_dir = here::here("results", "simulations", "non_survey", paste0("fold_", sprintf("%03d", ifold), "_partials"))
  if (!dir.exists(partial_dir)) dir.create(partial_dir, recursive = TRUE)

  temp = sim_settings %>% filter(fold == ifold) # set settings for this simulation

  # check for completed iters
  completed_iters = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$") %>%
    str_extract("\\d+") %>%
    as.integer()

  for (iter in 1:nsim) {
    if (iter %in% completed_iters && !force_iter) {
      message("Skipping completed iteration: ", iter)
      next
    }

    set.seed(iter)
    x = try({
      data = generate_population_nonsurvey(
        n = temp$n,
        L = temp$L,
        scenario = temp$scenario,
        family = temp$family,
        seed = iter,
        snr_eps = temp$sigma
      )
      sample_data = data$sample_data

      tic()
      b_hat =
        sample_data %>%
        get_betatilde(., family = temp$family, type = "unweighted") %>%
        get_betahat(., L = temp$L)

      boots = run_boots_xfast(data = sample_data,
                            betaHat = b_hat,
                            family = temp$family,
                            seed = 2025,
                            L = temp$L,
                            num_boots = temp$num_boots,
                            boot_type = "unweighted")

      ci = get_cis(boots, betaHat = b_hat, smooth_for_ci = TRUE, L = temp$L,
                   mult_fac = 1.2)
      elapsed_time = toc(quiet = TRUE)
      stats_fui = get_coverage_stats_fui(ci, beta_true = data$beta_true, L = temp$L) %>%
        mutate(time = unname(elapsed_time$toc - elapsed_time$tic))
      rm(boots)

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
          family = temp$family,
          method = "fREML",
          bs.yindex = list(
            bs = "ps",
            k = 15,
            m = c(2, 1)
          )
        )
      elapsed_time = toc(quiet = TRUE)
      stats_famm = get_coverage_stats_famm(fit_pffr, beta_true = data$beta_true, L = temp$L) %>%
        mutate(time = unname(elapsed_time$toc - elapsed_time$tic))

      ## bayes

      stats =
        stats_famm %>%
        mutate(type = "FAMM") %>%
        bind_rows(stats_fui %>% mutate(type = "FUI"))

      write_rds(stats, file = file.path(partial_dir, sprintf("iter_%03d.rds", iter)))

    })
    rm(x)
    print(iter)
  }


  partial_files = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$", full.names = TRUE)

  sim_res = map(partial_files, read_rds) %>%
    keep(., is.data.frame) %>%
    list_rbind(names_to = "id")

  write_rds(sim_res, file = outfile)


  unlink(partial_dir, recursive = TRUE)
}
