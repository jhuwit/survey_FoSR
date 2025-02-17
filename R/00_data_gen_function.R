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


## specify true fixed effects functions
# grid <- seq(0, 1, length = L)
# beta_true <- matrix(NA, 2, L) # time-dependent beta-true trajectory
# scenario <- 1 ## indicate the scenario to use
# if(scenario == 1){
#   beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
#   beta_true[2,] = dnorm(grid, .6, .15)/20
# }else if(scenario == 2){
#   beta_true[1,] <- 0.53 + 0.06*sin(3*grid*pi) - 0.03*cos(6.5*grid*pi)
#   beta_true[2,] <- dnorm(grid, 0.2, .1)/60 + dnorm(grid, 0.35, .1)/200 -
#     dnorm(grid, 0.65, .06)/250 + dnorm(grid, 1, .07)/60
# }
# rownames(beta_true) <- c("Intercept", "x")
#
# I = 10e6 # superpopulation size
# set.seed(4574)
# X_des = cbind(1, rnorm(I, 0, 2)) # design matrix: Intercept + predicted variables
# true_fixef <- X_des %*% beta_true
# L = 50



# Y_obs <- matrix(NA, nrow = I, ncol = L)
# for(i in 1:nrow(true_fixef)){
#   for(j in 1:L){
#     Y_obs[i, j] <- rnorm(1, mean = true_fixef[i, j], sd = SNR_sigma)
#   }
# }


# set.seed(4574)
# Y_obs <- matrix(
#   rnorm(n = nrow(true_fixef) * L, mean = as.vector(t(true_fixef)), sd = SNR_sigma),
#   nrow = nrow(true_fixef),
#   ncol = L,
#   byrow = TRUE
# )
#


simulate_survey_fosr <- function(X_des, # design matrix
                                 beta_true, # true beta
                                 true_fixef, # true fixed effects (Y)
                                 Y_obs, # y matrix
                                 I_n = 100, # subjects in each psu-strata combination
                                 L = 50, # dimension of the functional domain
                                 num_strata = 30,  # total strata
                                 sampled_strata = 30
){

  I = nrow(X_des)
  # Generate strata sizes using a Dirichlet distribution
  dirichlet_probs = gtools::rdirichlet(1, rep(1, num_strata)) # vector that sums to 1, of length num_strata
  if(sampled_strata == num_strata){ # if we are taking all of the strata, then selected strata are all
    selected_strata = 1:num_strata
  } else { # otherwise sample the selected strata proportional to size
    selected_strata = sample(1:num_strata, size = sampled_strata, replace = FALSE, prob = dirichlet_probs) %>% sort
  }
  stratum_assignments = sample(1:num_strata, size = I, replace = TRUE, prob = dirichlet_probs) # assign each individual in entire population to one of the strata

  p_strata_design <- as.vector(dirichlet_probs)

  # For strata that were not selected, set their probability to zero.
  p_strata_design[!(1:num_strata %in% selected_strata)] <- 0

  # Renormalize so that the selected strata probabilities sum to 1.
  p_strata_design <- p_strata_design / sum(p_strata_design)


  # Second-stage: Individual Sampling within each strata
  final_sample <- c() # vector for indices of final sample
  p_i_given_strata <- rep(NA, I) # probabilities of selecting each individual within strata
  psu = c() # vector for PSUs

  for (strata in selected_strata) {
    inds <- which(stratum_assignments == strata)
    inclusion_probs <- exp(X_des[inds, 2] / 2)  # Compute individual inclusion probabilities based on covariate, dividing by 2 to reduce skewness
    inclusion_probs <- inclusion_probs / sum(inclusion_probs) # normalize
    inclusion_probs <- inclusion_probs * I_n * 2 # Scale to get the desired number of individuals in each strata / PSU combo
    # summary(inclusion_probs) # print inclusion probabilities
    inclusion_probs[inclusion_probs > 1] <- 1
    sampled_in_strata <- inds[rbinom(length(inds), 1, inclusion_probs) == 1] # sample individuals based on inclusion probabilities
    final_sample <- c(final_sample, sampled_in_strata) # add indices to final sample
    p_i_given_strata[inds] <- inclusion_probs # store inclusion probabilities for each individual
    psu_strata <- rbinom(length(sampled_in_strata), 1, 0.5) + 1 # randomly assign PSUs in strata
    psu <- c(psu, psu_strata)
  }

  # Compute final individual selection probabilities
  p_i <- p_strata_design[stratum_assignments] * p_i_given_strata
  survey_weights <- 1 / p_i  # final survey weights

  dat.sim <- data.frame(
    ID = final_sample,
    X = X_des[final_sample,2],
    strata = stratum_assignments[final_sample], # Assume same strata here ??
    psu = psu,
    weight = survey_weights[final_sample]
  ) %>%
    filter(!is.na(weight)) %>%
    mutate(weight = weight / mean(weight))

  Y_sample <- data.frame(Y_obs[final_sample,])
  colnames(Y_sample) <- paste0("Y", 1:L)

  sample_data <- cbind(dat.sim, Y_sample)

  return(sample_data)
}
