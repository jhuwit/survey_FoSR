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
source(here::here("R", "utils.R"))

settings = expand_grid(slope = "random",
                       snr_b = c(2.5, 5),
                       snr_eps = c(.125, .25, .5),
                       psu_re = 0.25,
                       I_n =c(100, 250, 500))

ifold = get_fold()
ifold = 10
## issue is some strata only have 1 PSU - need to fix!
# ifold = 10
temp = settings[ifold, ]
I = 10e6
L = 50
family = "gaussian"
seed = 4576
num_strata = sampled_strata = 30
num_selected_psu = 2
strata_sigma = 0.05
psu_sigma = sqrt(strata_sigma ^ 2 * temp$psu_re)
sd_beta = 0.125
sampled_strata = 30
I_n = temp$I_n
snr_b = temp$snr_b
snr_eps = temp$snr_eps
nsim = 50

set.seed(seed)
X_des <- cbind(1, rnorm(I, 0, 2))

grid <- seq(0, 1, length = L)
beta_fixed <- matrix(NA, 2, L)

beta_fixed[1, ] <- -0.15 - 0.1 * sin(2 * pi * grid) - 0.1 * cos(2 * pi * grid)
beta_fixed[2, ] <- dnorm(grid, 0.6, 0.15) / 20


rownames(beta_fixed) <- c("Intercept", "x")

nbasis <- 5
basis <- fda::create.bspline.basis(c(0, 1), nbasis)
Phi <- fda::eval.basis(grid, basis)
## Step two: assign PSU and strata
set.seed(seed)
dirichlet_probs <- gtools::rdirichlet(1, rep(1, num_strata))
set.seed(seed)
stratum_assignments <- sample(1:num_strata, I, replace = TRUE, prob = dirichlet_probs)
psu_assignments <- rep(NA, I)

for (s in 1:num_strata) {
  # s = 1
  set.seed(seed + s)
  num_in_strata = sum(stratum_assignments == s)
  num_psu = round(runif(1, 100, 200), 0)
  set.seed(seed + s)
  dps = gtools::rdirichlet(1, rep(1, num_psu))
  set.seed(seed + s)
  psu_in_stratum = sample(1:num_psu,
                          num_in_strata,
                          replace = TRUE,
                          prob = dps)
  psu_assignments[stratum_assignments == s] <- paste0(s, "_", psu_in_stratum)
}


set.seed(seed)
strata_scores <- matrix(rnorm(num_strata * nbasis, 0, strata_sigma), num_strata, nbasis)
strata_random_effects <- strata_scores %*% t(Phi)

total_psu = length(unique(psu_assignments))

set.seed(seed)
psu_scores <- matrix(
  rnorm(total_psu * nbasis, 0, psu_sigma),
  total_psu,
  nbasis
)
psu_random_effects <- psu_scores %*% t(Phi)



strata_effects_indiv <- strata_random_effects[stratum_assignments, ]
psu_effects_indiv <- psu_random_effects[as.numeric(factor(psu_assignments)), ]
random_effects <- strata_effects_indiv + psu_effects_indiv

rm(strata_effects_indiv, psu_effects_indiv)
## Add stratum-specific slope modifications
set.seed(seed)
stratum_scaling <- rnorm(num_strata, mean = 1, sd = sd_beta)

beta1_by_stratum <- matrix(rep(stratum_scaling, each = L), nrow = num_strata) *
  matrix(rep(beta_fixed[2, ], times = num_strata),
         nrow = num_strata,
         byrow = TRUE)

# Assign to individuals
beta1_by_indiv <- beta1_by_stratum[stratum_assignments, ]


if (temp$slope == "fixed") {
  fixef_signal <- matrix(rep(beta_fixed[1, ], I), nrow = I, byrow = TRUE) +
    X_des[, 2] * beta1_by_indiv
  ranef <- sd(fixef_signal) / sd(random_effects) / snr_b * random_effects
  rm(random_effects)
} else {
  fixef_signal <- matrix(rep(beta_fixed[1, ], I), nrow = I, byrow = TRUE) +
    X_des[, 2] * matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)

  # New: include stratum-specific slope variation in the random effects
  slope_re <- (stratum_scaling[stratum_assignments] - 1) *
    matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)
  ranef <- slope_re + random_effects
  ranef <- sd(fixef_signal) / sd(ranef) / snr_b * ranef
  rm(random_effects)
}

sd_lp = sd(as.vector(t(fixef_signal)) + as.vector(t(ranef)))
sigma = sd_lp / snr_eps

# generate outcomes
set.seed(seed)
Y_obs =  matrix(
  rnorm(
    n = I * L,
    mean = as.vector(t(fixef_signal)) + as.vector(t(ranef)),
    sd = sigma
  ),
  nrow = I,
  ncol = L,
  byrow = TRUE
)

X1 = X_des[, 2]
sim_res = list()
iter = 1
for (iter in 1:nsim) {
  set.seed(iter)
  selected_strata = sample(1:num_strata,
                           sampled_strata,
                           replace = FALSE,
                           prob = dirichlet_probs)


  p_strata_design <- dirichlet_probs
  p_strata_design[!(1:num_strata %in% selected_strata)] <- 0
  p_strata_design <- p_strata_design / sum(p_strata_design)

  final_sample <- c()
  p_i_given_strata <- rep(NA, I)
  psus <- c()

  for (strata in 1:num_strata) {
    inds_in_stratum <- which(stratum_assignments == strata)

    # Get PSU sizes in this stratum
    psu_sizes <- table(psu_assignments[inds_in_stratum])
    psu_ids <- names(psu_sizes)

    # Sample PSUs WITH replacement using PPS
    selected_psus <- sample(psu_ids,
                            size = num_selected_psu,
                            replace = TRUE,
                            prob = psu_sizes)

    psu_probs <- psu_sizes / sum(psu_sizes)  # PPS

    psu_prob_selected <- 1 - (1 - psu_probs)^num_selected_psu
    names(psu_prob_selected) <- psu_ids
    for (psu in selected_psus) {
      inds_in_psu <- which(psu_assignments == psu &
                             stratum_assignments == strata)

      # Informative sampling based on X1 or PC scores
      incl_score <- X1[inds_in_psu]  # Or your combination like: -0.6 * score1 + ...
      inclusion_probs <- plogis(incl_score)
      inclusion_probs <- inclusion_probs / sum(inclusion_probs)
      inclusion_probs <- inclusion_probs * I_n * 2
      inclusion_probs[inclusion_probs > 1] <- 1

      sampled_units <- inds_in_psu[rbinom(length(inds_in_psu), 1, inclusion_probs) == 1]

      final_sample <- c(final_sample, sampled_units)
      psus <- c(psus, rep(psu, length(sampled_units)))
      p_psu <- psu_prob_selected[psu]

      p_i_given_strata[inds_in_psu] <- p_psu * inclusion_probs
    }
  }
  if(num_strata != sampled_strata){
    p_i <- p_strata_design[stratum_assignments] * p_i_given_strata
  } else{
    p_i = p_i_given_strata
  }
  survey_weights <- 1 / p_i

  dat.sim <- data.frame(
    ID = final_sample,
    X = X1[final_sample],
    strata = stratum_assignments[final_sample],
    psu = sub(".*\\_", "", psus),
    weight = survey_weights[final_sample]
  ) %>%
    filter(!is.na(weight)) %>%
    mutate(weight = weight / mean(weight))

  Y_sample <- data.frame(Y_obs[final_sample, ])
  colnames(Y_sample) <- paste0("Y", 1:L)

  data =  cbind(dat.sim, Y_sample)
  rm(dat.sim, Y_sample)
  parallel = TRUE
  ncores = parallelly::availableCores() - 2
  boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'weighted', 'none')
  fit_types = c('design', 'design', 'weighted', 'none')
  options(survey.lonely.psu = "adjust")

  betaTilde = map(
    .x = fit_types,
    .f = get_betatilde,
    data = data,
    family = family
  )

  betaHat = map(.x = betaTilde, get_betahat)

  if(parallel){
    plan(multisession, workers = ncores)
    res = future_map2(
      .x = boot_types,
      .y = betaHat,
      .f = run_boots,
      data = data,
      samp_stages = c("PPSWR", "Poisson"),
      family = family,
      # num_boots = 100,
      num_boots = 500,
      set_seed = TRUE,
      seed = 2025,
      L = 50,
      .options = furrr_options(seed = TRUE)
    )


    # run_boots(data = data, boot_type = "Rao-Wu-Yue-Beaumont", betaHat = betaHat[[1]],
    #           num_boots = 10, set_seed = TRUE, samp_stages = c("PPSWR", "Poisson"))
    #
    #
    # run_boots(data = data, boot_type = "BRR", betaHat = betaHat[[1]],
    #           num_boots = 10, set_seed = TRUE, samp_stages = c("PPSWR", "Poisson"))
    cis = future_map2(
      .x = res,
      .y = betaHat,
      .f = get_cis,
      smooth_for_ci = TRUE,
      .options = furrr_options(seed = TRUE)
    )
    plan(sequential)
  } else{
    res = map2(
      .x = boot_types,
      .y = betaHat,
      .f = run_boots,
      data = data,
      family = family,
      # num_boots = 100,
      num_boots = 500,
      set_seed = TRUE,
      seed = 2025,
      L = 50
    )

    cis = map2(
      .x = res,
      .y = betaHat,
      .f = get_cis,
      smooth_for_ci = TRUE
    )
  }



  stats = map2(
    .x = cis,
    .y = boot_types,
    .f = get_coverage_stats,
    beta_true = beta_fixed
  ) %>%
    list_rbind()

  rm(cis, res, betaHat, betaTilde)
  sim_res[[iter]] <- stats
  print(iter)
}
sim_res_df = bind_rows(sim_res, .id = "iter")

fold = sprintf("%02d", ifold)
if(!dir.exists(here::here("results", "snr"))){
  dir.create(here::here("results", "snr"), recursive = TRUE)
}
write_rds(sim_res_df, here::here("results", "snr", paste0("data2_fold_", fold, ".rds")))


#####
# get_plot_df = function(sim_res, tname){
#   num_var = nrow(sim_res$betaHat)
#   map_dfr(.x = 1:num_var,
#           .f = function(r){
#             data.frame(s = 1:length(sim_res$betaHat[r,]),
#                        beta = sim_res$betaHat[r,],
#                        lower = sim_res$betaHat[r, ] - 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
#                        upper = sim_res$betaHat[r, ] + 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
#                        lower.joint = sim_res$betaHat[r, ] - sim_res$qn[r] *
#                          sqrt(diag(sim_res$betaHat.var[, , r])),
#                        upper.joint = sim_res$betaHat[r, ] + sim_res$qn[r] * sqrt(diag(sim_res$betaHat.var[, , r]))) %>%
#               mutate(name = sim_res$betaHat %>% rownames() %>% .[r])
#           }) %>%
#     mutate(varname = tname)
# }
# plot_objs = map2(.x = cis, .y = boot_types, .f = get_plot_df) %>% bind_rows()
#
# btrue =
#   beta_fixed %>%
#   as_tibble() %>%
#   mutate(beta = c("beta_0", "beta_1")) %>%
#   pivot_longer(cols = -beta) %>%
#   mutate(x = as.numeric(sub(".*V", "", name)),
#          type = "betaTrue")
#
# btrue_df =
#   btrue %>%
#   mutate(name = if_else(beta == "beta_0", "(Intercept)", "X")) %>%
#   select(-beta) %>%
#   rename(s = x, beta = value)
# plot_objs %>%
#   ggplot(aes(x = s, y = beta)) +
#   geom_line()  +
#   facet_grid(varname~name) +
#   geom_ribbon(aes(x = s, ymin = lower, ymax = upper), alpha = .3,
#               fill = "gray10") +
#   geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint), alpha = .3,
#               fill = "gray20") +
#   labs(x = "L", y = "Estimate") +
#   theme_classic() +
#   geom_line(data = btrue_df, aes(x = s, y = beta), color = "red")
#
