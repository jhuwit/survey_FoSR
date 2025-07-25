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

generate_superpopulation = function(I = 10e6, # size of superpopulation
                                    L = 50, # length of functional domain
                                    scenario = 1, # shape of function
                                    family = "gaussian",
                                    seed = 4574,
                                    num_strata = 30, # num total strata
                                    strata_sigma = 0.05,
                                    psu_sigma = 0.025,
                                    sd_beta = 0.125, # strata scaling factor
                                    snr_b = 1, # signal noise ratio for random to fixed effects
                                    snr_eps = 1, # signal to noise for gaussian
                                    strata_d = 1,
                                    psu_d = 10
){
  stopifnot("scenario must be either 1 or 2" = scenario %in% c(1, 2))
  stopifnot("family must be either 'gaussian', 'poisson', or 'binomial'" = family %in% c("gaussian", "poisson", "binomial"))


  X_des <- cbind(1, rnorm(I, 0, 2))

  ## simulate true beta
  grid <- seq(0, 1, length = L)
  beta_fixed <- matrix(NA, 2, L)
  beta_fixed[1, ] <- -0.15 - 0.1 * sin(2 * pi * grid) - 0.1 * cos(2 * pi * grid)
  beta_fixed[2, ] <- dnorm(grid, 0.6, 0.15) / 20
  rownames(beta_fixed) <- c("Intercept", "x")

  ## assign individuals to PSU and strata
  set.seed(seed)
  dirichlet_probs <- gtools::rdirichlet(1, rep(strata_d, num_strata))
  set.seed(seed)
  stratum_assignments <- sample(1:num_strata, I, replace = TRUE, prob = dirichlet_probs)
  psu_assignments <- rep(NA, I)

  for (s in 1:num_strata) {
    set.seed(seed + s)
    num_in_strata = sum(stratum_assignments == s)
    num_psu = round(runif(1, 75, 125), 0)
    set.seed(seed + s)
    dps = gtools::rdirichlet(1, rep(psu_d, num_psu))
    set.seed(seed + s)
    psu_in_stratum = sample(1:num_psu,
                            num_in_strata,
                            replace = TRUE,
                            prob = dps)
    psu_assignments[stratum_assignments == s] <- paste0(s, "_", psu_in_stratum)
  }

  ## create psu and strata-specific random effects
  nbasis <- 5
  basis <- fda::create.bspline.basis(c(0, 1), nbasis)
  Phi <- fda::eval.basis(grid, basis)

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

  ## add stratum-specific slope modifications
  set.seed(seed)
  stratum_scaling <- rnorm(num_strata, mean = 1, sd = sd_beta)

  beta1_by_stratum <- matrix(rep(stratum_scaling, each = L), nrow = num_strata) *
    matrix(rep(beta_fixed[2, ], times = num_strata),
           nrow = num_strata,
           byrow = TRUE)

  # assign to individuals
  beta1_by_indiv <- beta1_by_stratum[stratum_assignments, ]

  beta1_mean = colMeans(beta1_by_indiv)
  argvals = 1:L
  kk = L / 2 + 1
  beta1_mean_sm = gam(beta1_mean ~ s(argvals, bs = "tp", k = kk),
                      method = "GCV.Cp")$fitted.values




  # adjust random effect based on signal to noise parameters

  fixef_signal <- matrix(rep(beta_fixed[1, ], I), nrow = I, byrow = TRUE) +
    X_des[, 2] * matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)

  # include stratum-specific slope variation in the random effects
  slope_re <- (stratum_scaling[stratum_assignments] - 1) *
    matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)
  ranef <- slope_re + random_effects
  ranef <- sd(fixef_signal) / sd(ranef) / snr_b * ranef
  rm(random_effects)


  # generate outcomes
  set.seed(seed)
  if (family == "gaussian") {
    sd_lp = sd(as.vector(t(fixef_signal)) + as.vector(t(ranef)))
    sigma = sd_lp / snr_eps

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
  } else if(family == "binomial") {
    p_true = plogis(as.vector(fixef_signal) + as.vector(t(ranef)))
    Y_obs <- matrix(
      rbinom(
        n = I * L,
        size = 1,
        prob = as.vector(t(p_true))
      ),
      nrow = I,
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "poisson"){
    lam_true = exp(as.vector(fixef_signal) + as.vector(t(ranef)))
    Y_obs <- matrix(
      rpois(n = I * L,
            lambda = as.vector(t(lam_true))),
      nrow = I,
      ncol = L,
      byrow = TRUE
    )
  }

  return(list(Y_obs = Y_obs,
              X_des = X_des,
              stratum_assignments = stratum_assignments,
              psu_assignments = psu_assignments,
              dirichlet_probs = dirichlet_probs,
              beta_true = beta_fixed,
              b1_mean = beta1_mean_sm))
}
theme_set(theme_light())
nsim = 5

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





L = 50
nknots_min = NULL
nknots_min_cov = 35
mult_fac = 1.2
nknots <- min(round(L / 2), nknots_min)
nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
nknots_fpca <- min(round(L / 2), 35)
argvals = 1:L
iter = 1

settings = expand_grid(strata_d = c(1, 2,3, 4),
                       psu_d = c(1, 2, 3, 4, 5, 7, 10))
# pdf(here::here("wtd_debug_test.pdf"))
get_bias = function(iter, lst, fit_types = c("weighted", "unweighted")){
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
  data$weight = pmin(data$weight, quantile(data$weight, .99))

  beta_true = lst$beta_true
  b1m = lst$b1_mean
  betaTilde = map(
    .x = fit_types,
    .f = get_betatilde,
    data = data,
    family = temp$family
  )

  betaHat = map(.x = betaTilde, get_betahat, L = temp$len)
  # get bias
  bias = map_dbl(.x = betaHat,
                 .f = \(x) sum(abs(x[2,] - b1m)^2))

  bias
}

get_bias_all = function(row) {
  tmp = settings[row,]
  s_d = tmp$strata_d
  p_d = tmp$psu_d
  lst = generate_superpopulation(
    scenario = temp$scen,
    family = temp$family,
    I = 10e5,
    L = temp$len,
    psu_sigma = psu_sigma,
    snr_b = temp$snr_b,
    snr_eps = temp$snr_eps,
    strata_d = s_d,
    psu_d = p_d
  )
  beta_true = lst$beta_true
  res = map(.x = 1:nsim, .f = get_bias, lst = lst) %>%
    map(., .f = \(x) t(x) %>% as_tibble()) %>%
    bind_rows() %>%
    magrittr::set_colnames(c("weighted", "unweighted")) %>%
    mutate(iter = 1:nsim,
           setting = row)
  res
}

final_res = map(.x = 1:nrow(settings),
                .f = get_bias_all)


write_rds(final_res,
          here::here("debug_wtd_b1m_trunc.rds"))


