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


temp = settings[ifold, ]
I = 10e6
L = 50
family = "gaussian"
seed = 4576
num_strata = 30
psu_per_strata = 2
strata_sigma = 0.05
psu_sigma = sqrt(strata_sigma ^ 2 * temp$psu_re)
sd_beta = 0.125
sampled_strata = 30
I_n = temp$I_n
snr_b = temp$snr_b
snr_eps = temp$snr_eps
nsim = 200

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
X_des <- cbind(1, rnorm(I, 0, 2))

set.seed(seed)
dirichlet_probs <- gtools::rdirichlet(1, rep(1, num_strata))
set.seed(seed)
stratum_assignments <- sample(1:num_strata, I, replace = TRUE, prob = dirichlet_probs)

psu_assignments <- rep(NA, I)
for (s in 1:num_strata) {
  set.seed(seed + s)
  psus_in_stratum <- sample(1:psu_per_strata, sum(stratum_assignments == s), replace = TRUE)
  psu_assignments[stratum_assignments == s] <- paste0(s, "_", psus_in_stratum)
}

set.seed(seed)
strata_scores <- matrix(rnorm(num_strata * nbasis, 0, strata_sigma), num_strata, nbasis)
strata_random_effects <- strata_scores %*% t(Phi)

set.seed(seed)
psu_scores <- matrix(
  rnorm(num_strata * psu_per_strata * nbasis, 0, psu_sigma),
  num_strata * psu_per_strata,
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
I = nrow(X_des)
sim_res = list()
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
  for (strata in selected_strata) {
    inds <- which(stratum_assignments == strata)

    # Induce informative sampling based on X and PC scores
    # incl_score <- -0.6 * score1[inds] + 0.3 * score2[inds] + 0.4 * score3[inds] + 0.6 * X1[inds]
    incl_score <-  X1[inds]

    inclusion_probs <- plogis(incl_score)
    inclusion_probs <- inclusion_probs / sum(inclusion_probs)
    inclusion_probs <- inclusion_probs * I_n * 2
    inclusion_probs[inclusion_probs > 1] <- 1

    sampled_in_strata <- inds[rbinom(length(inds), 1, inclusion_probs) == 1]
    final_sample <- c(final_sample, sampled_in_strata)
    psus <- c(psus, psu_assignments[sampled_in_strata])
    p_i_given_strata[inds] <- inclusion_probs
  }

  p_i <- p_strata_design[stratum_assignments] * p_i_given_strata
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
  ncores = parallelly::availableCores() - 1
  boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'Preston', 'weighted', 'none')
  fit_types = c('design', 'design', 'design', 'weighted', 'none')
  options(survey.lonely.psu = "adjust")

  betaTilde = map(
    .x = fit_types,
    .f = get_betatilde,
    data = data,
    family = family
  )

  betaHat = map(.x = betaTilde, get_betahat)

  plan(multisession, workers = ncores)
  res = future_map2(
    .x = boot_types,
    .y = betaHat,
    .f = run_boots,
    data = data,
    family = family,
    # num_boots = 100,
    num_boots = 500,
    set_seed = TRUE,
    seed = 2025,
    L = 50,
    .options = furrr_options(seed = TRUE)
  )

  cis = future_map2(
    .x = res,
    .y = betaHat,
    .f = get_cis,
    smooth_for_ci = TRUE,
    .options = furrr_options(seed = TRUE)
  )
  plan(sequential)


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
write_rds(sim_res_df, here::here("results", "snr", paste0("fold_", fold, ".rds")))
