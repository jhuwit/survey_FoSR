library(future)
library(furrr)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(here)
library(gtools)
library(devtools)
library(survey)
library(progress)
library(lme4)
library(mgcv)
library(mvtnorm)
library(refund)
library(svrep)

source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "00_data_gen_function_ml.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "create_survey_settings_ml.R"))


force = FALSE
force_iter = FALSE
B = 500
# ncores = parallelly::availableCores()
ncores = future::availableCores()

fit_types = c('weighted', 'unweighted')
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'weighted', 'unweighted')
options(survey.lonely.psu = "adjust")


ifold = get_fold()


temp = settings[ifold,]


nsim = 200
outfile = here::here("results", "simulations", "survey_sim_ml", paste0("fold_", sprintf("%03d", ifold), ".rds"))
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
partial_dir = here::here("results", "simulations", "survey_sim_ml", paste0("fold_", sprintf("%03d", ifold), "_partials"))
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



lst = generate_superpopulation_ml(
  family = temp$family,
  I = 10e6,
  L = temp$len,
  snr_b = temp$snr_b,
  snr_eps = temp$snr_eps,
  strata_sigma = temp$strata_sigma,
  strata_scale = temp$strata_scale,
  seed = 111,
  snr_u = temp$snr_u
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

plan(multisession, workers = ncores)


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
    data = sample_from_population_wor_ml(
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
      family = temp$family,
      subj_vec = lst$subj_vec,
      visits_per_subj =  lst$visits_per_subj
    )
    print(summary(data$weight))


    ##################
    # do some basic data vis
    ##################

    # data |>
    #   filter(ID %in% unique(data$ID)[1:5]) |>
    #   pivot_longer(cols = starts_with("Y")) |>
    #   mutate(s = sub(".*Y", "", name) |> as.numeric()) |>
    #   ggplot(aes(x = s, y = value, color = factor(ID), group = visit)) +
    #   # geom_point() +
    #   geom_smooth(se = FALSE) +
    #   facet_wrap(~ID)



    betaTilde = map(
      .x = fit_types,
      .f = get_betatilde,
      data = data,
      family = temp$family
    )

    betaHat = map(.x = betaTilde, get_betahat, L = temp$len)


    res = future_map(
      .x = boot_types,
      .f = run_boots_xfast,
      betaHat = betaHat[[1]], # don't actually use betahat in calculation so we can just use the weighted one
      data = data,
      family = temp$family,
      num_boots = B,
      seed = 2025,
      L = temp$len,
      samp_stages = c("PPSWOR", "Poisson"),
      .options = furrr_options(seed = TRUE)
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

    stats = future_map2(
      .x = cis,
      .y = boot_types,
      .f = get_coverage_stats,
      beta_true = bh_data,
      L = temp$len
    ) %>%
      list_rbind() %>%
      mutate(n = nrow(data), n_boot = B)


    write_rds(stats, file = file.path(partial_dir, sprintf("iter_%03d.rds", iter)))

  })
  rm(x)
  print(iter)
}
plan(sequential)


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

