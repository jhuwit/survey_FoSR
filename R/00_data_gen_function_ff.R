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


generate_superpopulation = function(I = 10e6, # size of superpopulation
                                    L = 50, # length of functional domain
                                    family = "gaussian",
                                    seed = 4574,
                                    num_strata = 30, # num total strata
                                    strata_sigma = 0.05,
                                    psu_factor = 0.5,
                                    strata_scale = 0.125, # strata scaling factor
                                    snr_b = 1, # signal noise ratio for random to fixed effects
                                    snr_eps = 1 # signal to noise for gaussian
){
  stopifnot("family must be either 'gaussian', 'poisson', or 'binomial'" = family %in% c("gaussian", "poisson", "binomial"),
            "I must be greater than 1000" = I > 1000,
            "I must be greater than or equal to 25" = L >= 25,
            "num_strata must be greater than 0" = num_strata > 0,
            "strata_sigma, psu_factor and strata_scale must be greater than or equal to 0" = all(c(strata_sigma, psu_factor, strata_scale) >= 0),
            "signal to noise ratios must be greater than 0" = snr_b > 0 | is.na(snr_b),
            "signal to noise ratios must be greater than 0" = snr_eps > 0 | is.na(snr_eps))


  set.seed(seed)
  X_des  = cbind(1, rnorm(I, 0, 2))

  ## simulate true beta based on scenarios
  grid  = seq(0, 1, length = L)
  beta_fixed  = matrix(NA, 2, L)
  beta_fixed[1, ]  = -0.15 - 0.1 * sin(2 * pi * grid) - 0.1 * cos(2 * pi * grid)
  beta_fixed[2, ]  = dnorm(grid, 0.6, 0.15) / 20

  rownames(beta_fixed)  = c("Intercept", "x")

  ## assign individuals to strata
  set.seed(seed)
  dirichlet_probs  = gtools::rdirichlet(1, rep(4, num_strata)) # generate dirichlet probabilities
  set.seed(seed)
  stratum_assignments  = sample(1:num_strata, I, replace = TRUE, prob = dirichlet_probs) # generate stratum assignments
  psu_assignments  = rep(NA, I)

  # assign individuals to PSUs - between 75 and 125 psus per stratum
  for (s in 1:num_strata) {
    set.seed(seed + s)
    num_in_strata = sum(stratum_assignments == s)
    num_psu = round(runif(1, 75, 125), 0)
    set.seed(seed + s)
    dps = gtools::rdirichlet(1, rep(10, num_psu))
    set.seed(seed + s)
    psu_in_stratum = sample(1:num_psu,
                            num_in_strata,
                            replace = TRUE,
                            prob = dps)
    psu_assignments[stratum_assignments == s]  = paste0(s, "_", psu_in_stratum)
  }

  # case where there's no strata-specific noise
  if (strata_sigma == 0 & strata_scale == 0){
    lin_pred  = matrix(rep(beta_fixed[1, ], I), nrow = I, byrow = TRUE) + X_des[, 2] * matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)
  } else if (strata_sigma == 0 & strata_scale > 0) {
    set.seed(seed)
    stratum_scaling  = rnorm(num_strata, mean = 1, sd = strata_scale)

    beta1_by_stratum  = matrix(rep(stratum_scaling, each = L), nrow = num_strata) *
      matrix(rep(beta_fixed[2, ], times = num_strata),
             nrow = num_strata,
             byrow = TRUE)

    # assign to individuals
    beta1_by_indiv  = beta1_by_stratum[stratum_assignments, ]
    fixef_signal  = matrix(rep(beta_fixed[1, ], I), nrow = I, byrow = TRUE) +
      X_des[, 2] * matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)

    slope_re  = (stratum_scaling[stratum_assignments] - 1) *
      matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)
    ranef  = slope_re
    ranef  = sd(as.vector(fixef_signal)) / sd(as.vector(ranef)) / snr_b * ranef
    rm(slope_re)
    lin_pred = fixef_signal + ranef

  } else if (strata_sigma > 0 & strata_scale == 0) {
    psu_sigma = sqrt(strata_sigma ^ 2 * psu_factor)

    ## create psu and strata-specific random effects
    nbasis  = 5
    basis  = fda::create.bspline.basis(c(0, 1), nbasis)
    Phi  = fda::eval.basis(grid, basis)

    set.seed(seed)
    strata_scores  = matrix(rnorm(num_strata * nbasis, 0, strata_sigma), num_strata, nbasis)
    strata_random_effects  = strata_scores %*% t(Phi)

    total_psu = length(unique(psu_assignments))

    set.seed(seed)
    psu_scores  = matrix(
      rnorm(total_psu * nbasis, 0, psu_sigma),
      total_psu,
      nbasis
    )
    psu_random_effects  = psu_scores %*% t(Phi)


    strata_effects_indiv  = strata_random_effects[stratum_assignments, ]
    psu_effects_indiv  = psu_random_effects[as.numeric(factor(psu_assignments)), ]
    random_effects  = strata_effects_indiv + psu_effects_indiv

    rm(strata_effects_indiv, psu_effects_indiv)

    fixef_signal  = matrix(rep(beta_fixed[1, ], I), nrow = I, byrow = TRUE) +
      X_des[, 2] * matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)
    ranef  = sd(as.vector(fixef_signal)) / sd(as.vector(random_effects)) / snr_b * random_effects
    rm(random_effects)
    lin_pred = fixef_signal + ranef

  } else { # random effects and slope modification
    psu_sigma = sqrt(strata_sigma ^ 2 * psu_factor)

    ## create psu and strata-specific random effects
    nbasis  = 5
    basis  = fda::create.bspline.basis(c(0, 1), nbasis)
    Phi  = fda::eval.basis(grid, basis)

    set.seed(seed)
    strata_scores  = matrix(rnorm(num_strata * nbasis, 0, strata_sigma), num_strata, nbasis)
    strata_random_effects  = strata_scores %*% t(Phi)

    total_psu = length(unique(psu_assignments))

    set.seed(seed)
    psu_scores  = matrix(
      rnorm(total_psu * nbasis, 0, psu_sigma),
      total_psu,
      nbasis
    )
    psu_random_effects  = psu_scores %*% t(Phi)


    strata_effects_indiv  = strata_random_effects[stratum_assignments, ]
    psu_effects_indiv  = psu_random_effects[as.numeric(factor(psu_assignments)), ]
    random_effects  = strata_effects_indiv + psu_effects_indiv

    rm(strata_effects_indiv, psu_effects_indiv)

    ## add stratum-specific slope modifications
    set.seed(seed)
    stratum_scaling  = rnorm(num_strata, mean = 1, sd = strata_scale)

    beta1_by_stratum  = matrix(rep(stratum_scaling, each = L), nrow = num_strata, byrow = TRUE) *
      matrix(rep(beta_fixed[2, ], times = num_strata),
             nrow = num_strata,
             byrow = TRUE)

    # assign to individuals
    beta1_by_indiv  = beta1_by_stratum[stratum_assignments, ]

    # adjust random effect based on signal to noise parameters
    fixef_signal  = matrix(rep(beta_fixed[1, ], I), nrow = I, byrow = TRUE) +
      X_des[, 2] * matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)

    # include stratum-specific slope variation in the random effects
    slope_re  = (stratum_scaling[stratum_assignments] - 1) *
      matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)
    ranef  = slope_re + random_effects
    ranef  = sd(as.vector(fixef_signal)) / sd(as.vector(ranef)) / snr_b * ranef
    rm(random_effects)
    lin_pred = fixef_signal + ranef
  }


  # temp = matrix(c(1, 2, 3, 4, 5, 6), ncol = 3, byrow = TRUE) # making sure I am doing dimensions right!
  # matrix(rnorm(n = 6, mean = as.vector(t(temp)), sd = .001), nrow = 2, ncol = 3, byrow = TRUE)

  # temp = matrix(c(10, -5, 10, -5, -5, 10), ncol = 3, byrow = TRUE) # making sure I am doing dimensions right!
  # matrix(rbinom(n = 6, size = 1, prob = plogis(as.vector(t(temp)))), nrow = 2, ncol = 3, byrow = TRUE)
  # rbinom(n = 6, size = 1, prob = plogis(as.vector(t(temp))))

  # lin pred is n x L
  # generate outcomes
  if (family == "gaussian") {
    sd_lp = sd(as.vector(lin_pred))
    sigma = sd_lp / snr_eps
    set.seed(seed)
    Y_obs = matrix(
      rnorm(n = I * L,
            mean = as.vector(t(lin_pred)),
            sd = sigma), # need to use t to put in correct order
      nrow = I,
      ncol = L,
      byrow = TRUE
    )
  } else if(family == "binomial") {
    p_true = plogis(as.vector(t(lin_pred)))
    set.seed(seed)
    Y_obs  = matrix(
      rbinom(
        n = I * L,
        size = 1,
        prob = p_true
      ),
      nrow = I,
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "poisson"){
    lam_true = exp(as.vector(t(lin_pred)))
    set.seed(seed)
    Y_obs  = matrix(
      rpois(n = I * L,
            lambda = lam_true),
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
              beta_true = beta_fixed))
}


get_p_i = function(i, probs) probs[i] * (1 + sum((probs[-i]) / (1-probs[-i])))


sample_from_population_wor = function(X_des, # design matrix
                                      Y_obs, # y matrix
                                      I_n = 500, # subjects in each psu-strata combination
                                      num_strata = 30,  # total strata
                                      stratum_assignments, # assignment to ea strata
                                      num_selected_psu = 2,
                                      dirichlet_probs, # stratum probabilities
                                      psu_assignments, # assignment to psu (w/in strata)
                                      L = 50, # length of fnl domain
                                      seed = 1,
                                      inf_level = 1,
                                      compression = 3,
                                      family = "gaussian"
){
  I = nrow(X_des)
  X1 = X_des[, 2]
  # select strata (use all for now)
  selected_strata = 1:num_strata
  p_strata_design  = dirichlet_probs[selected_strata]

  final_sample  = c() # to store final sample ids
  p1  = rep(NA, I) # stage 1 selection probability
  p2  = rep(NA, I) # stage 2 selection probability
  p_overall  = rep(NA, I) # overall selection probability
  psus  = c() # to store PSUs

  # within each strata, select PSU with replacement using PPS
  for (strata in 1:num_strata) {
    # strata = 1
    # individuals in the strata
    inds_in_stratum  = which(stratum_assignments == strata)

    # Get PSU sizes in this stratum
    psu_sizes = table(psu_assignments[inds_in_stratum])
    psu_ids = names(psu_sizes)

    # Sample PSUs WITH replacement using PPS
    set.seed(strata + seed)
    selected_psus = sample(psu_ids,
                           size = num_selected_psu,
                           replace = FALSE,
                           prob = psu_sizes)

    psu_probs  = psu_sizes / sum(psu_sizes)  # PPS
    # probability of selection is 1 - (p(not selected both times))
    # psu_prob_selected  = 1 - (1 - psu_probs) ^ num_selected_psu
    # names(psu_prob_selected)  = psu_ids

    psu_prob_selected = map_dbl(.x = match(selected_psus, psu_ids),
                                .f = get_p_i,
                                psu_probs)

    names(psu_prob_selected)  = selected_psus
    # within each selected PSU select individuals based on X1
    for (psu in selected_psus) {
      inds_in_psu  = which(psu_assignments == psu &
                             stratum_assignments == strata)
      if (inf_level == 0) {
        # Uniform sampling
        n = length(inds_in_psu)
        inclusion_probs = rep(1 / n, n)
      } else {
        # Compute mean outcome in PSU
        y_mean = rowMeans(Y_obs[inds_in_psu, ])

        # Compute inclusion score depending on family
        incl_score = switch(family,
                             "gaussian" = y_mean * inf_level,
                             "poisson"  = log(y_mean) * inf_level,
                             "binomial" = qlogis(pmin(pmax(y_mean, 1e-6), 1 - 1e-6)) * inf_level,
                             stop("Unknown family"))

        # Apply compression and map to probabilities
        score_compressed = pmax(pmin(incl_score, compression), -compression)
        inclusion_probs = plogis(score_compressed)


      }

      inclusion_probs_adj  = inclusion_probs / sum(inclusion_probs) * I_n
      inclusion_probs_adj[inclusion_probs_adj > 1]  = 1

      set.seed(strata + seed + which(selected_psus == psu)) # ensure reproducibility
      sampled_units  = inds_in_psu[rbinom(length(inds_in_psu), 1, inclusion_probs_adj) == 1]

      final_sample  = c(final_sample, sampled_units)
      psus  = c(psus, rep(psu, length(sampled_units)))
      p_psu  = psu_prob_selected[which(names(psu_prob_selected) == psu)]

      p1[inds_in_psu]  = p_psu
      p2[inds_in_psu]  = inclusion_probs_adj
      p_overall[inds_in_psu]  = p_psu * inclusion_probs_adj
    }
  }


  survey_weights  = 1 / p_overall


  dat.sim  = data.frame(
    ID = final_sample,
    X = X1[final_sample],
    strata = stratum_assignments[final_sample],
    psu = sub(".*\\_", "", psus),
    weight = survey_weights[final_sample],
    p_stage1 = p1[final_sample],
    p_stage2 = p2[final_sample]
  )


  Y_sample  = data.frame(Y_obs[final_sample, ])
  colnames(Y_sample)  = paste0("Y", 1:L)

  data =  cbind(dat.sim, Y_sample)
  return(data)
}

generate_population_nonsurvey = function(n = 500, # size of superpopulation
                               L = 50, # length of functional domain
                               scenario = 1, # shape of function
                               family = "gaussian",
                               seed = 4574,
                               snr_eps = 1
){
  stopifnot("scenario must be either 1 or 2" = scenario %in% c(1, 2))
  stopifnot("family must be either 'gaussian', 'poisson', or 'binomial'" = family %in% c("gaussian", "poisson", "binomial"))

  grid = seq(0, 1, length = L)
  beta_true = matrix(NA, 2, L) # time-dependent beta-true trajectory

  if (scenario == 1) {
    beta_true[1, ] = -0.15 - 0.1 * sin(2 * grid * pi) - 0.1 * cos(2 * grid *
                                                                    pi)
    beta_true[2, ] = dnorm(grid, .6, .15) / 20
  } else if (scenario == 2) {
    beta_true[1, ] = 0.53 + 0.06 * sin(3 * grid * pi) - 0.03 * cos(6.5 * grid *
                                                                     pi)
    beta_true[2, ] = dnorm(grid, 0.2, .1) / 60 + dnorm(grid, 0.35, .1) /
      200 -
      dnorm(grid, 0.65, .06) / 250 + dnorm(grid, 1, .07) / 60
  }
  rownames(beta_true) <- c("Intercept", "x")
  set.seed(seed)
  X_des = cbind(1, rnorm(n, 0, 2)) # design matrix: Intercept + predicted variables
  true_fixef <- X_des %*% beta_true

  set.seed(seed + 1)
  if (family == "gaussian") {
    sd_lp = sd(as.vector(true_fixef))
    sigma = sd_lp / snr_eps
    Y_obs <- matrix(
      rnorm(
        n = nrow(true_fixef) * L,
        mean = as.vector(t(true_fixef)),
        sd = sigma
      ),
      nrow = nrow(true_fixef),
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "binomial") {
    p_true = plogis(true_fixef)
    Y_obs <- matrix(
      rbinom(
        n = nrow(p_true) * L,
        size = 1,
        prob = as.vector(t(p_true))
      ),
      nrow = nrow(p_true),
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "poisson") {
    lam_true <- exp(true_fixef)
    Y_obs <- matrix(
      rpois(n = nrow(lam_true) * L, lambda = as.vector(t(lam_true))),
      nrow = nrow(lam_true),
      ncol = L,
      byrow = TRUE
    )
  }

  dat.sim <- data.frame(
    ID = 1:n,
    X = X_des[,2]
  )

  # Generate fixed effects and true Y.
  Y_df <- data.frame(Y_obs)
  colnames(Y_df) <- paste0("Y", 1:L)

  sample_data <- cbind(dat.sim, Y_df)

  return(list(sample_data = sample_data, beta_true = beta_true))
}

