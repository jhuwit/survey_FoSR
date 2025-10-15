## empirical sim
library(tidyverse)
source(here::here("R_cp", "01_sim_functions.R"))
source(here::here("R", "utils.R"))
library(future)
library(furrr)
force = FALSE
pa_df = read_rds(here::here("data", "mims_covariates.rds"))
L = 1440
n_sim = 1000
n_boot = 500
compression = 2
nsim = 200

settings = expand_grid(iter = 1:nsim,
                       weight_type = c("uniform", "nh_weights", "mims_weights", "combo_weights"))
ifold = get_fold()
weight_curr = settings[ifold,]$weight_type
iter = settings[ifold,]$iter

outfile = here::here("results", "simulations", "empirical_sim_1440", paste0("fold_", sprintf("%03d", ifold), ".rds"))

if(!file.exists(outfile) || force) { # if the file doesn't exist or we want to force re-do it
  if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)

  fit_types = c('weighted', 'unweighted')
  boot_types = c('BRR', 'weighted', 'unweighted')
  options(survey.lonely.psu = "adjust")


  means =
    pa_df %>%
    pivot_longer(cols = starts_with("min")) %>%
    group_by(SEQN,
             data_release_cycle,
             psu,
             strata,
             gender,
             age = age_in_years_at_screening,
             weight = full_sample_2_year_mec_exam_weight) %>%
    summarize(mean_mims = mean(value, na.rm = TRUE), .groups = "drop")

  means =
    means %>%
    mutate(weight = weight / 2,
           weight_sc = scale(weight, center = TRUE, scale = TRUE))

  dat_full = tibble(ID = means$SEQN, X = means$gender)


  pa_means =
    pa_df %>%
    select(SEQN, starts_with("min"))

  Y = as.matrix(pa_means %>% select(starts_with("min")))
  colnames(Y) = paste0("Y", 1:L)

  full_data = cbind(dat_full, Y)


  if (weight_curr == "uniform") {
    inclusion_probs = rep(1 / nrow(means), nrow(means)) # uniform sampling
  } else if (weight_curr == "nh_weights") {
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


  pik = sampling::inclusionprobabilities(inclusion_probs, n_sim)
  set.seed(iter)
  ind = as.logical(sampling::UPpoisson(pik))
  # which(ind) %>% length

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
    pa_means %>%
    filter(SEQN %in% dat_sample$ID) %>%
    select(starts_with("min"))

  colnames(Y_sample) = paste0("Y", 1:L)

  ss_data = cbind(dat_sample, Y_sample)

  ## fit model on full data to get "true" effect
  bt = get_betatilde(data = full_data, family = "gaussian", type = "unweighted")
  beta_true = get_betahat(bt, L = L, nknots_min = 15) ## these are the true fns
  # beta_true %>%
  #   as_tibble() %>%
  #   mutate(type = c("Intercept", "X")) %>%
  #   janitor::clean_names() %>%
  #   pivot_longer(cols = starts_with("x")) %>%
  #   mutate(xval = as.numeric(sub(".*x", "", name))) %>%
  #   ggplot(aes(x = xval, y = value)) +
  #   geom_smooth() +
  #   facet_wrap(.~type, scales = "free_y")


  betaTilde = map(
    .x = fit_types,
    .f = get_betatilde,
    data = ss_data,
    family = "gaussian"
  )

  betaHat = map(.x = betaTilde, get_betahat, L = L, nknots_min = 15)


  res = map(
    .x = boot_types,
    .f = run_boots,
    data = ss_data,
    weights = ss_data$weight,
    family = "gaussian",
    num_boots = n_boot,
    seed = 2025,
    parallel = TRUE,
    n_cores = parallelly::availableCores() -1
  )


  fit_list = list(betaHat[[1]], betaHat[[1]], betaHat[[2]])

  cis = map2(
    .x = res,
    .y = fit_list,
    .f = get_cis_par,
    L = L,
    smooth_for_ci = TRUE,
    nknots_min = 15,
    nknots_min_fpca = 15,
    parallel = TRUE,
    n_cores = parallelly::availableCores() -1
  )

  stats = map2(
    .x = cis,
    .y = boot_types,
    .f = get_coverage_stats,
    beta_true = beta_true,
    L = L
  ) %>%
    list_rbind() %>%
    mutate(n = nrow(data),
           iter = iter)



    write_rds(stats, file = outfile)

}
