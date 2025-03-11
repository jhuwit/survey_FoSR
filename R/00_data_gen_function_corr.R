required_packages <- c("svrep", "purrr", "progress", "svylme", 'dplyr', 'readr', 'refund', 'gtools')

install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

library(here)
library(devtools)
library(fastFMM)
library(haven)
library(dplyr)
library(survey)
library(progress)
library(lme4)
library(mgcv)
library(ggplot2)
library(gridExtra)

invisible(lapply(required_packages, install_and_load))



generate_superpopulation = function(I = 10e6, # size of superpopulation
                                    L = 50, # length of functional domain
                                    scenario = 1, # shape of function
                                    family = "gaussian",
                                    seed = 4574,
                                    num_strata = 30, # num total strata
                                    psu_per_strata = 2, # num PSU per strata
                                    SNR_sigma = 0.5, # signal to noise ratio (only applicable in gaussian case)
                                    sigma_strata = 0.2, # strata noise
                                    sigma_psu = 0.1 # psu noise
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
  nbasis <- 5  # Number of basis functions
  basis <- fda::create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis)
  Phi <- fda::eval.basis(grid, basis)


  set.seed(seed)
  strata_scores <- matrix(rnorm(num_strata * nbasis, mean = 0, sd = sigma_strata), num_strata, nbasis)
  set.seed(seed)
  psu_scores <- matrix(rnorm(num_strata * psu_per_strata * nbasis, mean = 0, sd = sigma_psu),
                       num_strata * psu_per_strata, nbasis)


  set.seed(seed)
  X_des = cbind(1, rnorm(I, 0, 2)) # design matrix: Intercept + predicted variables
  true_fixef <- X_des %*% beta_true

  # Assign strata
  set.seed(seed)
  dirichlet_probs <- gtools::rdirichlet(1, rep(1, num_strata))  # strata for inclusion in ea strata
  set.seed(seed)
  stratum_assignments <- sample(1:num_strata, size = I, replace = TRUE, prob = dirichlet_probs) # assign individuals to strata based on probs


  # Assign PSUs within strata
  psu_assignments <- rep(NA, I)
  for (s in 1:num_strata) {
    psus_in_stratum <- sample(1:psu_per_strata, sum(stratum_assignments == s), replace = TRUE)
    psu_assignments[stratum_assignments == s] <- paste0(s, "_", psus_in_stratum)
  }
  # psu_assignments = "strata_psu" character format
  strata_random_effects <- strata_scores %*% t(Phi)  # Shape: num_strata x L
  psu_random_effects <- psu_scores %*% t(Phi)  # Shape: num_PSU x L

  # Map strata and PSU effects to individuals
  strata_effects_indiv <- strata_random_effects[stratum_assignments, ]
  psu_effects_indiv <- psu_random_effects[as.numeric(factor(psu_assignments)), ]
  random_effects <- strata_effects_indiv + psu_effects_indiv

  set.seed(seed)
  if (family == "gaussian") {
    Y_obs <- matrix(
      rnorm(
        n = nrow(true_fixef) * L,
        mean = as.vector(t(true_fixef)) + random_effects,
        sd = SNR_sigma
      ),
      nrow = nrow(true_fixef),
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "binomial") {
    p_true = plogis(true_fixef + random_effect)
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
    lam_true <- exp(true_fixef + random_effects)
    Y_obs <- matrix(
      rpois(n = nrow(lam_true) * L, lambda = as.vector(t(lam_true))),
      nrow = nrow(lam_true),
      ncol = L,
      byrow = TRUE
    )
  }
  return(list(Y_obs = Y_obs,
              X_des = X_des,
              stratum_assignments = stratum_assignments,
              psu_assignments = psu_assignments,
              dirichlet_probs = dirichlet_probs,
              beta_true = beta_true))
}


simulate_survey_fosr_corr <- function(X_des, # design matrix
                                      Y_obs, # y matrix
                                      I_n = 100, # subjects in each psu-strata combination
                                      num_strata = 30,  # total strata
                                      sampled_strata = 30, # strata to be included in sample
                                      stratum_assignments, # assignment to ea strata
                                      psu_assignments, # assignment to psu (w/in strata)
                                      dirichlet_probs, # strata probabilities
                                      L = 50 # length of fnl domain
){

  I = nrow(X_des)
  if(sampled_strata == num_strata){ # if we are taking all of the strata, then selected strata are all
    selected_strata = 1:num_strata
  } else { # otherwise sample the selected strata proportional to size
    selected_strata = sample(1:num_strata, size = sampled_strata, replace = FALSE, prob = dirichlet_probs) %>% sort
  }

  p_strata_design <- as.vector(dirichlet_probs)

  # For strata that were not selected, set their probability to zero.
  p_strata_design[!(1:num_strata %in% selected_strata)] <- 0

  # Renormalize so that the selected strata probabilities sum to 1.
  p_strata_design <- p_strata_design / sum(p_strata_design)

  # Second-stage: Individual Sampling within each strata
  final_sample <- c() # vector for indices of final sample
  p_i_given_strata <- rep(NA, I) # probabilities of selecting each individual within strata
  psus <- c()

  for (strata in selected_strata) {
    inds <- which(stratum_assignments == strata)
    inclusion_probs <- exp(X_des[inds, 2] / 2)  # Compute individual inclusion probabilities based on covariate, dividing by 2 to reduce skewness
    inclusion_probs <- inclusion_probs / sum(inclusion_probs) # normalize
    inclusion_probs <- inclusion_probs * I_n * 2 # Scale to get the desired number of individuals in each strata / PSU combo
    # summary(inclusion_probs) # print inclusion probabilities
    inclusion_probs[inclusion_probs > 1] <- 1
    sampled_in_strata <- inds[rbinom(length(inds), 1, inclusion_probs) == 1] # sample individuals based on inclusion probabilities
    final_sample <- c(final_sample, sampled_in_strata) # add indices to final sample
    psus <- c(psus, psu_assignments[sampled_in_strata])
    p_i_given_strata[inds] <- inclusion_probs # store inclusion probabilities for each individual
  }

  # Compute final individual selection probabilities
  p_i <- p_strata_design[stratum_assignments] * p_i_given_strata
  survey_weights <- 1 / p_i  # final survey weights



  dat.sim <- data.frame(
    ID = final_sample,
    X = X_des[final_sample,2],
    strata = stratum_assignments[final_sample],
    psu = sub(".*\\_", "", psus), # we just want psu, not psu/strata combo
    weight = survey_weights[final_sample]
  ) %>%
    filter(!is.na(weight)) %>%
    mutate(weight = weight / mean(weight))

  Y_sample <- data.frame(Y_obs[final_sample,])
  colnames(Y_sample) <- paste0("Y", 1:L)

  sample_data <- cbind(dat.sim, Y_sample)

  return(sample_data)
}
