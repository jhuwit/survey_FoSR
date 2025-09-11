## empirical sim
library(tidyverse)
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "utils.R"))
library(future)
library(furrr)
force = FALSE
force_iter = FALSE
pa_df = read_rds(here::here("sim_data", "pa_df_persub.rds"))
L = 96
n_sim = 1000
n_boot = 500
compression = 2
nsim = 200

ifold = get_fold()

run_one_simulation = function(iter){
  pik = sampling::inclusionprobabilities(inclusion_probs, n_sim)
  set.seed(iter)
  ind = as.logical(sampling::UPpoisson(pik))
  which(ind) %>% length

  means = means %>%
    mutate(p_i = pik)

  sample_df =
    means %>% slice(which(ind)) %>%
    select(SEQN, psu, strata, gender, age, p_i) %>%
    mutate(weight = 1/p_i)

  dat_sample =
    tibble(
      ID = sample_df$SEQN,
      X = sample_df$gender,
      weight = sample_df$weight,
      psu = sample_df$psu,
      strata = sample_df$strata
    )

  Y_sample =
    pa_means_hour %>%
    filter(SEQN %in% dat_sample$ID) %>%
    select(starts_with("hour"))

  colnames(Y_sample) = paste0("Y", 1:L)

  ss_data = cbind(dat_sample, Y_sample)

  ## fit model on full data to get "true" effect
  bt = get_betatilde(data = full_data, family = "gaussian", type = "unweighted")
  beta_true = get_betahat(bt, L = L) ## these are the true fns
  #
  # bh %>%
  #   as_tibble() %>%
  #   mutate(type = c("Intercept", "X")) %>%
  #   janitor::clean_names() %>%
  #   pivot_longer(cols = starts_with("x")) %>%
  #   mutate(xval = as.numeric(sub(".*x", "", name))) %>%
  #   ggplot(aes(x = xval, y = value)) +
  #   geom_smooth() +
  #   facet_wrap(.~type, scales = "free_y")
  ## now fit on sample
  betaTilde = map(
    .x = fit_types,
    .f = get_betatilde,
    data = ss_data,
    family = "gaussian"
  )

  betaHat = map(.x = betaTilde, get_betahat, L = L)


  res = map(
    .x = boot_types,
    .f = run_boots_fast,
    betaHat = betaHat[[1]],
    # don't actually use betahat in calculation so we can just use the weighted one
    data = ss_data,
    family = "gaussian",
    num_boots = n_boot,
    seed = 2025,
    L = L,
    samp_stages = c("Poisson")
  )

  fit_list = list(betaHat[[1]], betaHat[[1]], betaHat[[2]])

  cis = map2(
    .x = res,
    .y = fit_list,
    .f = get_cis,
    L = L,
    smooth_for_ci = TRUE
  )

  stats = map2(
    .x = cis,
    .y = boot_types,
    .f = get_coverage_stats,
    beta_true = beta_true,
    L = L
  ) %>%
    list_rbind() %>%
    mutate(n = nrow(data))

  stats
}

weight_types = c("uniform", "nh_weights", "mims_weights", "combo_weights")
fit_types = c('weighted', 'unweighted')
boot_types = c('BRR', 'weighted', 'unweighted')
options(survey.lonely.psu = "adjust")

weight_curr = weight_types[ifold]
outfile = here::here("results", "simulations", "empirical_sim", paste0("fold_", sprintf("%03d", ifold), ".rds"))

if(!file.exists(outfile) || force) { # if the file doesn't exist or we want to force re-do it
  if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)

  partial_dir = here::here("results", "simulations", "empirical_sim", paste0("fold_", sprintf("%03d", ifold), "_partials"))
  if (!dir.exists(partial_dir)) dir.create(partial_dir)


  means =
    pa_df %>%
    pivot_longer(cols = starts_with("min")) %>%
    group_by(SEQN, data_release_cycle, psu, strata, gender,
             age = age_in_years_at_screening,
             weight = full_sample_2_year_mec_exam_weight) %>%
    summarize(mean_mims = mean(value, na.rm = TRUE),
              .groups = "drop")

  means =
    means %>%
    mutate(weight = weight / 2,
           weight_sc = scale(weight, center = TRUE, scale = TRUE))

  dat_full = tibble(
    ID = means$SEQN,
    X = means$gender)


  pa_means_hour =
    pa_df %>%
    select(SEQN, starts_with("min")) %>%
    pivot_longer(cols = starts_with("min")) %>%
    mutate(min = sub(".*min\\_", "", name) %>% as.numeric) %>%
    mutate(hour = ceiling(min / 15)) %>%
    group_by(hour, SEQN) %>%
    summarize(mims = sum(value, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = "hour", values_from = "mims", names_prefix = "hour_")
  # need to arrange first
  Y = as.matrix(pa_means_hour %>% select(starts_with("hour")))
  colnames(Y) = paste0("Y", 1:L)

  full_data = cbind(dat_full, Y)


  if(weight_curr == "uniform"){
    inclusion_probs = rep(1 / nrow(means), nrow(means)) # uniform sampling
  } else if(weight_curr == "nh_weights"){
    incl_score = means$weight_sc
    score_compressed  = pmax(pmin(incl_score, compression), -1 * compression)
    inclusion_probs = plogis(score_compressed)
  } else if (weight_curr == "mims_weights") {
    incl_score = scale(means$mean_mims, center = TRUE, scale = TRUE)
    score_compressed  = pmax(pmin(incl_score, compression), -1 * compression)
    inclusion_probs = plogis(score_compressed)
  } else if (weight_curr == "combo_weights") {
    incl_score = (0.5 * means$weight_sc) +
      (0.5 * scale(means$mean_mims, center = TRUE, scale = TRUE))
    score_compressed  = pmax(pmin(incl_score, compression), -1 * compression)
    inclusion_probs = plogis(score_compressed)
  }

  plan(multisession, workers = 9)

  sim_res_list = future_map(
    .x = 1:nsim,
    .f = run_one_simulation,
    .options = furrr_options(seed = TRUE)
  )


  sim_res = sim_res_list %>%
    keep(., is.data.frame) %>%
    list_rbind(names_to = "id")

  write_rds(sim_res, file = outfile)


}

