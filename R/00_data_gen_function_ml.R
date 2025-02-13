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
library(mvtnorm)
library(splines)
invisible(lapply(required_packages, install_and_load))



simulate_survey_fosr_ml <- function(X_des, # design matrix
                                 X_df, # design matrix with ID and subject indices
                                 eta_true, # true fixed effects with random eff
                                 Y_obs,
                                 J_subj, # num visits per person
                                 I_n = 100, # subjects in each psu-strata combination
                                 alpha = 2, # number PSU
                                 L = 50, # dimension of the functional domain
                                 num_strata = 30,  # total strata
                                 sampled_strata = 30
){
  # do the selection on the ID level
  # summarize X variable by ID
  X_summ =
    X_df %>%
    group_by(id) %>%
    summarize(X2 = mean(V2))

  # Generate strata sizes using a Dirichlet distribution
  dirichlet_probs = gtools::rdirichlet(1, rep(1, num_strata)) # vector that sums to 1, of length num_strata
  if(sampled_strata == num_strata){ # if we are taking all of the strata, then selected strata are all
    selected_strata = 1:num_strata
  } else { # otherwise sample the selected strata proportional to size
    selected_strata = sample(1:num_strata, size = sampled_strata, replace = FALSE, prob = dirichlet_probs) %>% sort
  }
  stratum_assignments = sample(1:num_strata, size = I, replace = TRUE, prob = dirichlet_probs) # assign each individual in entire population to one of the strata

  p_strata = table(stratum_assignments) / sum(table(stratum_assignments))
  p_strata = c(p_strata) %>% unname()
  p_strata[which(!(1:num_strata %in% selected_strata))] <- 0 # if not selected, replace with 0
  w_strata = 1 / p_strata


  # Second-stage: Individual Sampling within each strata
  final_sample <- c() # vector for indices of final sample
  p_i_given_strata <- rep(NA, I) # probabilities of selecting each individual within strata
  psu = c() # vector for PSUs

  for (strata in selected_strata) {
    inds <- which(stratum_assignments == strata)
    inclusion_probs <- exp(X_summ$X2[inds])
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
  p_i <- p_strata[stratum_assignments] * p_i_given_strata
  survey_weights <- 1 / p_i  # Final survey weights (Inf: Not selected)

  # get the indices associated with the final sample ids
  indices =
    X_df %>%
    filter(id %in% final_sample) %>%
    pull(index)


  X_des_sample = X_des[indices,]
  eta_sample = eta_true[indices,]

  visits = X_df %>%
    filter(id %in% final_sample) %>%
    pull(visit)

  ids_long = X_df %>%
    filter(id %in% final_sample) %>%
    pull(id)

  Y_sample = Y_obs[indices,]

  strata = rep(stratum_assignments[final_sample], J_subj[final_sample]) # repeat strata assignments based on visit nums
  weights = rep(survey_weights[final_sample], J_subj[final_sample]) # repeat weights based on visit nums
  psu = rep(psu, J_subj[final_sample])
  ## here
  dat.sim <- data.frame(ID = ids_long,
                        visit = visits,
                        X = X_des_sample[,2],
                        strata = strata,
                        weight = weights,
                        psu = psu,
                        eta = I(eta_sample),
                        Y = I(Y_sample))


  sample_data <- dat.sim %>%
    filter(!is.na(weight)) # shouldn't happen but just in case

  return(sample_data)
}
