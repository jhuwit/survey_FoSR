## smoothing parameter sensitivity analysis

library(future)
library(furrr)
library(purrr)
library(here)
library(readr)
library(stringr)
library(devtools)
# library(fastFMM)
# library(haven)
library(survey)
library(gtools)
library(progress)
library(lme4)
library(mgcv)
# library(tidyfun)
library(mvtnorm)
library(refund)
library(svrep)

source(here::here("R", "create_survey_settings.R"))
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "01_sim_functions_smoothsens.R"))
source(here::here("R", "00_data_gen_function_ff.R"))
source(here::here("R", "00_data_gen_function_fast.R"))
source(here::here("R", "utils.R"))


force = FALSE
force_iter = FALSE
B = 500
# ncores = parallelly::availableCores()
ncores = future::availableCores()

fit_types = c('weighted')
boot_types = c('BRR')
options(survey.lonely.psu = "adjust")


ifold = get_fold()

temp = smaller_grid[ifold,]

nsim = switch(as.character(temp$len),
              "50" = 200,
              "100" = 200,
              "1440" = 100)

outfile = here::here("results", "simulations", "survey_sim_smoothsens", paste0("fold_", sprintf("%03d", ifold), ".rds"))
if (file.exists(outfile)) {
  x = read_rds(outfile)
  if(nrow(x) >= nsim * length(boot_types) * 2) {
    message("Output file already has results for all iterations. Exiting.")
    quit(save = "no", status = 0)
  } else if (force) {
    message("Output file exists but force = TRUE, so re-running all iterations.")
  } else {
    message("Output file exists but does not have results for all iterations. Continuing to run missing iterations.")
  }}

if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)
partial_dir = here::here("results", "simulations", "survey_sim_smoothsens", paste0("fold_", sprintf("%03d", ifold), "_partials"))
if (!dir.exists(partial_dir)) dir.create(partial_dir)
partial_files = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$", full.names = TRUE)

if(length(partial_files) >= nsim) {
  sim_res = map(partial_files, \(x) read_rds(x) %>% mutate(id = sub(".*iter_(.+)\\.rds.*", "\\1", basename(x)))) %>%
    keep(., is.data.frame) %>%
    list_rbind()

  write_rds(sim_res, file = outfile)

  unlink(partial_dir, recursive = TRUE)
  message("All iterations complete, writing and exiting.")
  quit(save = "no", status = 0)
}

if(temp$len == 1440) {
  options(future.globals.maxSize = 100 * 1024^3)  # Equivalent to 1.13 GiB

  structure = generate_population_structure(
    I = 10e6,
    L = temp$len,
    family = temp$family,
    seed = 4574,
    num_strata = 30,
    strata_sigma = temp$strata_sigma,
    psu_factor = 0.5,
    strata_scale = temp$strata_scale,
    snr_b = temp$snr_b,
    snr_eps = temp$snr_eps
  )
  structure$global_sd_lp = compute_global_sd_lp_fast(structure)
  structure$snr_eps = temp$snr_eps  # match input

  for (iter in seq_len(nsim)) {
    if(file.exists(outfile)) {
      x = read_rds(outfile)
      c1 = unique(x$id) %>% as.integer()
      c2 = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$") %>%
        str_extract("\\d+") %>%
        as.integer()
      completed_iters = c(c1, c2) %>% unique()
    } else {
      completed_iters = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$") %>%
        str_extract("\\d+") %>%
        as.integer()
    }
    if (iter %in% completed_iters && !force_iter) {
      message("Skipping completed iteration: ", iter)
      next
    }

    set.seed(iter)
    x = try({
      data = sample_by_strata_stream_fast(
        structure,
        I_n = temp$In,
        num_selected_psu = 2,
        seed = iter,
        inf_level = temp$inf_level,
        compression = 2,
        family = temp$family
      )


      beta_true = structure$beta_fixed

       betaTilde = get_betatilde(data = data,
                                family = temp$family,
                                type = "weighted")

      betaHat = get_betahat_sensitivity(betaTilde,
                                        L = temp$len,
                                        nknots_spec = temp$nknots_spec,
                                        splines = temp$splines,
                                        method = temp$method)

      res = run_boots_xfast(data = data,
                            boot_type = "BRR",
                            betaHat = betaHat,
                            family = temp$family,
                            num_boots = B,
                            L = temp$len,
                            samp_stages = c("PPSWOR", "Poisson"))

      cis = get_cis_sensitivity(betaTilde_boot = res,
                                betaHat = betaHat,
                                L = temp$len,
                                smooth_for_ci = TRUE,
                                nknots_spec = temp$nknots_spec,
                                splines = temp$splines,
                                method = temp$method)

      stats = get_coverage_stats(mod_output = cis,
                                 beta_true = bh_data,
                                 name = iter,
                                 L = temp$len) |>
        mutate(n = nrow(data), n_boot = B)


      write_rds(stats, file = file.path(partial_dir, sprintf("iter_%03d.rds", iter)))

    })
    rm(x)
    print(iter)
  }
} else {

  lst = generate_superpopulation(
    family = temp$family,
    I = 10e4,
    L = temp$len,
    snr_b = temp$snr_b,
    snr_eps = temp$snr_eps,
    strata_sigma = temp$strata_sigma,
    strata_scale = temp$strata_scale,
    seed = 111
  )

  if(file.exists(here::here(partial_dir, "bh_data.rds"))) {
    bh_data = read_rds(here::here(partial_dir, "bh_data.rds"))
  } else {
    df = data.frame(X = lst$X_des[,2],
                    Y = lst$Y_obs)

    bt_data = get_betatilde(family = temp$family,
                            type = "unweighted",
                            data = df)
    rm(df); gc()
    bh_data  = get_betahat(bt_data, L = temp$len)

    write_rds(bh_data, here::here(partial_dir, "bh_data.rds"))
  }


  for (iter in seq_len(nsim)) {

    if(file.exists(outfile)) {
      x = read_rds(outfile)
      c1 = unique(x$id) %>% as.integer()
      c2 = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$") %>%
        str_extract("\\d+") %>%
        as.integer()
      completed_iters = c(c1, c2) %>% unique()
    } else {
      completed_iters = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$") %>%
        str_extract("\\d+") %>%
        as.integer()
    }


    if (iter %in% completed_iters && !force_iter) {
      message("Skipping completed iteration: ", iter)
      next
    }

    set.seed(iter)
    x = try({
      data = sample_from_population_wor(
        X_des = lst$X_des,
        Y_obs = lst$Y_obs,
        L = temp$len,
        I_n = temp$In,
        num_strata = 30,
        stratum_assignments = lst$stratum_assignments,
        psu_assignments = lst$psu_assignments,
        dirichlet_probs = lst$dirichlet_probs,
        seed = iter,
        inf_level = temp$inf_level,
        compression = 2,
        family = temp$family
      )
      print(summary(data$weight))

      betaTilde = get_betatilde(data = data,
                                family = temp$family,
                                type = "weighted")

      betaHat = get_betahat_sensitivity(betaTilde,
                                        L = temp$len,
                                        nknots_spec = temp$nknots_spec,
                                        splines = temp$splines,
                                        method = temp$method)

      res = run_boots_xfast(data = data,
                            boot_type = "BRR",
                            betaHat = betaHat,
                            family = temp$family,
                            num_boots = B,
                            L = temp$len,
                            samp_stages = c("PPSWOR", "Poisson"))

    cis = get_cis_sensitivity(betaTilde_boot = res,
                              betaHat = betaHat,
                              L = temp$len,
                              smooth_for_ci = TRUE,
                              nknots_spec = temp$nknots_spec,
                              splines = temp$splines,
                              method = temp$method)

      stats = get_coverage_stats(mod_output = cis,
                                 beta_true = bh_data,
                                 name = iter,
                                 L = temp$len) |>
        mutate(n = nrow(data), n_boot = B)



      write_rds(stats, file = file.path(partial_dir, sprintf("iter_%03d.rds", iter)))

    })
    rm(x)
    print(iter)
  }
}


partial_files = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$", full.names = TRUE)

sim_res = map(partial_files, \(x) read_rds(x) %>% mutate(id = sub(".*iter_(.+)\\.rds.*", "\\1", basename(x)))) %>%
  keep(., is.data.frame) %>%
  list_rbind()

if (file.exists(outfile)) {
  x = read_rds(outfile)
  sim_res = sim_res %>% bind_rows(x)
}
write_rds(sim_res, file = outfile)


unlink(partial_dir, recursive = TRUE)

