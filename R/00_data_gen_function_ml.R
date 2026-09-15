library(future)
library(furrr)
# library(tidyverse)
library(here)
library(devtools)
# library(fastFMM)
library(dplyr)
library(survey)
library(progress)
library(lme4)
# library(paletteer)
library(mgcv)
# library(ggplot2)
# library(gridExtra)
# library(tidyverse)
# library(tidyfun)
library(mvtnorm)
library(refund)
library(svrep)


generate_superpopulation_ml = function(I = 10e6, # size of superpopulation
                                    L = 50, # length of functional domain
                                    family = "gaussian",
                                    seed = 4574,
                                    num_strata = 30, # num total strata
                                    strata_sigma = 0.05,
                                    psu_factor = 0.5,
                                    strata_scale = 0.125, # strata scaling factor
                                    snr_b = 1, # signal noise ratio for random to fixed effects
                                    snr_eps = 1, # signal to noise for gaussian
                                    snr_u = 1, # signal to noise ratio for individual REs to strata/PSU REs
                                    J = 5 # mean # visits per person
){
  stopifnot("family must be either 'gaussian', 'poisson', or 'binomial'" = family %in% c("gaussian", "poisson", "binomial"),
            "I must be greater than 1000" = I > 1000,
            "L must be greater than or equal to 25" = L >= 25,
            "num_strata must be greater than 0" = num_strata > 0,
            "strata_sigma, psu_factor and strata_scale must be greater than or equal to 0" = all(c(strata_sigma, psu_factor, strata_scale) >= 0),
            "signal to noise ratios must be greater than 0" = snr_b > 0 | is.na(snr_b),
            "signal to noise ratios must be greater than 0" = snr_eps > 0 | is.na(snr_eps))
  set.seed(seed)
  visits_per_subj = pmax(runif(I, min = 0, max = 7), 1)
  subj_vec = rep(1:I, visits_per_subj)
  n = sum(visits_per_subj)

  set.seed(seed)
  X_des  = cbind(1, rnorm(n, 0, 2))

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

  stratum_assignments_visit = rep(stratum_assignments, times = visits_per_subj)
  psu_assignments_visit = rep(psu_assignments, times = visits_per_subj)

  # case where there's no strata-specific noise
  if (strata_sigma == 0 & strata_scale == 0){
    # lin_pred  = matrix(rep(beta_fixed[1, ], I), nrow = I, byrow = TRUE) + X_des[, 2] * matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)
    fixef_signal  = matrix(rep(beta_fixed[1, ], n), nrow = n, byrow = TRUE) + X_des[, 2] * matrix(rep(beta_fixed[2, ], n), nrow = n, byrow = TRUE)
    # dim(lin_pred)
  } else if (strata_sigma == 0 & strata_scale > 0) {
    set.seed(seed)
    stratum_scaling  = rnorm(num_strata, mean = 1, sd = strata_scale)

    beta1_by_stratum  = matrix(rep(stratum_scaling, each = L), nrow = num_strata) *
      matrix(rep(beta_fixed[2, ], times = num_strata),
             nrow = num_strata,
             byrow = TRUE)

    # assign to individuals
    beta1_by_indiv  = beta1_by_stratum[stratum_assignments_visit, ]
    fixef_signal  = matrix(rep(beta_fixed[1, ], n), nrow = n, byrow = TRUE) +
      X_des[, 2] * matrix(rep(beta_fixed[2, ], n), nrow = n, byrow = TRUE)

    ranef_strata_psu  = (stratum_scaling[stratum_assignments_visit] - 1) *
      matrix(rep(beta_fixed[2, ], n), nrow = n, byrow = TRUE)

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


    strata_effects_indiv  = strata_random_effects[stratum_assignments_visit, ]
    psu_effects_indiv  = psu_random_effects[as.numeric(factor(psu_assignments_visit)), ]
    ranef_strata_psu  = strata_effects_indiv + psu_effects_indiv

    rm(strata_effects_indiv, psu_effects_indiv)

    fixef_signal  = matrix(rep(beta_fixed[1, ], n), nrow = n, byrow = TRUE) +
      X_des[, 2] * matrix(rep(beta_fixed[2, ], n), nrow = n, byrow = TRUE)

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


    strata_effects_indiv  = strata_random_effects[stratum_assignments_visit, ]
    psu_effects_indiv  = psu_random_effects[as.numeric(factor(psu_assignments_visit)), ]
    ranef  = strata_effects_indiv + psu_effects_indiv

    rm(strata_effects_indiv, psu_effects_indiv)

    ## add stratum-specific slope modifications
    set.seed(seed)
    stratum_scaling  = rnorm(num_strata, mean = 1, sd = strata_scale)

    beta1_by_stratum  = matrix(rep(stratum_scaling, each = L), nrow = num_strata, byrow = TRUE) *
      matrix(rep(beta_fixed[2, ], times = num_strata),
             nrow = num_strata,
             byrow = TRUE)

    # assign to individuals
    beta1_by_indiv  = beta1_by_stratum[stratum_assignments_visit, ]

    # adjust random effect based on signal to noise parameters
    fixef_signal  = matrix(rep(beta_fixed[1, ], n), nrow = n, byrow = TRUE) +
      X_des[, 2] * matrix(rep(beta_fixed[2, ], n), nrow = n, byrow = TRUE)

    # include stratum-specific slope variation in the random effects
    slope_re  = (stratum_scaling[stratum_assignments_visit] - 1) *
      matrix(rep(beta_fixed[2, ], n), nrow = n, byrow = TRUE)
    ranef_strata_psu  = slope_re + ranef
    rm(slope_re, ranef)
  }


  psi_true = matrix(NA, 2, L)
  psi_true[1,] = (1.5 - sin(2*grid*pi) - cos(2*grid*pi) )
  psi_true[1,] = psi_true[1,] / sqrt(sum(psi_true[1,]^2))
  psi_true[2,] = sin(4*grid*pi)
  psi_true[2,] = psi_true[2,] / sqrt(sum(psi_true[2,]^2))

  set.seed(seed)
  c_true = mvtnorm::rmvnorm(I, mean = rep(0, 2), sigma = diag(c(3, 1.5))) ## simulate score function
  b_true = c_true %*% psi_true
  ranef_subj = b_true[subj_vec,]



  if (strata_sigma == 0 & strata_scale == 0){
    ranef_subj = sd(fixef_signal) / sd(ranef_subj) / snr_u * ranef_subj
    lin_pred = fixef_signal + ranef_subj
  } else {
    # first adjust ranef subj
    ranef_subj = sd(fixef_signal) / sd(ranef_subj) / snr_u * ranef_subj
    # then adjust ranef psu
    ranef_strata_psu = sd(fixef_signal) / sd(ranef_strata_psu) / snr_b * ranef_strata_psu
    lin_pred = fixef_signal + ranef_subj + ranef_strata_psu
  }

  # lin pred is n x L
  # generate outcomes
  if (family == "gaussian") {
    sd_lp = sd(as.vector(lin_pred))
    sigma = sd_lp / snr_eps
    set.seed(seed)
    Y_obs = matrix(
      rnorm(n = n,
            mean = as.vector(t(lin_pred)),
            sd = sigma), # need to use t to put in correct order
      nrow = n,
      ncol = L,
      byrow = TRUE
    )
  } else if(family == "binomial") {
    p_true = plogis(as.vector(t(lin_pred)))
    set.seed(seed)
    Y_obs  = matrix(
      rbinom(
        n = n,
        size = 1,
        prob = p_true
      ),
      nrow = n,
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "poisson"){
    lam_true = exp(as.vector(t(lin_pred)))
    set.seed(seed)
    Y_obs  = matrix(
      rpois(n = n,
            lambda = lam_true),
      nrow = n,
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
              subj_vec = subj_vec,
              visits_per_subj = visits_per_subj))
}


get_p_i = function(i, probs) probs[i] * (1 + sum((probs[-i]) / (1-probs[-i])))


sample_from_population_wor_ml = function(X_des, # design matrix
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
                                         family = "gaussian",
                                         subj_vec,
                                         visits_per_subj
){
  I = length(visits_per_subj)

  # get per subject means
  rs = rowsum(Y_obs, subj_vec)
  Y_means = rs / visits_per_subj

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
        y_mean = rowMeans(Y_means[inds_in_psu, ])

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

  # index to get the correct entries from original weight, psu, strata, etc.
  keep_rows = which(subj_vec %in% final_sample)
  Y_obs_selected = Y_obs[keep_rows,]
  subj_selected = subj_vec[keep_rows]

  psu_lookup = setNames(psus, final_sample)

  # get visit numbers
  visit_num = ave(rep(1L, length(subj_selected)), subj_selected, FUN = seq_along)


  dat.sim  = data.frame(
    ID = subj_selected,
    X = X_des[keep_rows, 2],
    strata = stratum_assignments[subj_selected],
    psu = sub(".*\\_", "", psu_lookup[as.character(subj_selected)]),
    weight = survey_weights[subj_selected],
    p_stage1 = p1[subj_selected],
    p_stage2 = p2[subj_selected],
    visit = visit_num
  )
  Y_sample = data.frame(Y_obs_selected)

  colnames(Y_sample)  = paste0("Y", 1:L)

  data =  cbind(dat.sim, Y_sample)
  return(data)
}
