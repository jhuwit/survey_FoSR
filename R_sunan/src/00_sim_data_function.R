### Ignore original survey setting
simulate_nhanes_survey_fosr <- function(X_des, # design matrix
                                        Y_obs, # y matrix
                                        I_n = 100, # subjects in each psu-strata combination
                                        L = 1440, # dimension of the functional domain
                                        num_strata = 30,  # total strata
                                        sampled_strata = 30){
  I = nrow(X_des)
  input_X_des = X_des
  # Generate strata sizes using a Dirichlet distribution
  dirichlet_probs = gtools::rdirichlet(1, rep(1, num_strata)) # vector that sums to 1, of length num_strata
  if(sampled_strata == num_strata){ # if we are taking all of the strata, then selected strata are all
    selected_strata = 1:num_strata
  } else { # otherwise sample the selected strata proportional to size
    selected_strata = sample(1:num_strata, size = sampled_strata, replace = FALSE, prob = dirichlet_probs) %>% sort
  }
  
  # Assign each individual i entire population to one of the strata
  stratum_assignments = sample(1:num_strata, size = I, replace = TRUE, prob = dirichlet_probs) 
  
  p_strata_design <- as.vector(dirichlet_probs)
  # For strata that were not selected, set their probability to zero.
  p_strata_design[!(1:num_strata %in% selected_strata)] <- 0
  # Renormalize so that the selected strata probabilities sum to 1.
  p_strata_design <- p_strata_design / sum(p_strata_design)
  
  # Second-stage: Individual Sampling within each strata
  final_sample <- c() # vector for indices of final sample
  p_i_given_strata <- rep(NA, I) # probabilities of selecting each individual within strata
  psu = c() # vector for PSUs
  X_des[,2] = scale(X_des[,2])
  
  for (strata in selected_strata) {
    inds <- which(stratum_assignments == strata)
    inclusion_probs <- exp(X_des[inds, 2] / 2)
    inclusion_probs <- inclusion_probs / sum(inclusion_probs) # normalize
    inclusion_probs <- inclusion_probs * I_n *2# Scale to get the desired number of individuals
    
    inclusion_probs[inclusion_probs > 1] <- 1
    sampled_in_strata <- inds[rbinom(length(inds), 1, inclusion_probs) == 1] 
    # sample individuals based on inclusion probabilities
    
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
    X = input_X_des[final_sample,2],
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



### Keep orginal survey weights
simulate_nhanes_survey_fosr_w <- function(X_des, # design matrix
                                          Y_obs, # y matrix
                                          I_n = 100, # subjects in each psu-strata combination
                                          L = 1440, # dimension of the functional domain
                                          num_strata = 30,  # total strata
                                          sampled_strata = 30){
  I = nrow(X_des)
  input_X_des = X_des
  # Generate strata sizes using a Dirichlet distribution
  dirichlet_probs = gtools::rdirichlet(1, rep(1, num_strata)) # vector that sums to 1, of length num_strata
  if(sampled_strata == num_strata){ # if we are taking all of the strata, then selected strata are all
    selected_strata = 1:num_strata
  } else { # otherwise sample the selected strata proportional to size
    selected_strata = sample(1:num_strata, size = sampled_strata, replace = FALSE, prob = dirichlet_probs) %>% sort
  }
  
  # Assign each individual i entire population to one of the strata
  stratum_assignments = sample(1:num_strata, size = I, replace = TRUE, prob = dirichlet_probs) 
  
  p_strata_design <- as.vector(dirichlet_probs)
  # For strata that were not selected, set their probability to zero.
  p_strata_design[!(1:num_strata %in% selected_strata)] <- 0
  # Renormalize so that the selected strata probabilities sum to 1.
  p_strata_design <- p_strata_design / sum(p_strata_design)
  
  # Second-stage: Individual Sampling within each strata
  final_sample <- c() # vector for indices of final sample
  p_i_given_strata <- rep(NA, I) # probabilities of selecting each individual within strata
  psu = c() # vector for PSUs
  X_des[,2] = scale(X_des[,2])
  
  for (strata in selected_strata) {
    inds <- which(stratum_assignments == strata)
    inclusion_probs <- exp( (X_des[inds, 3] + X_des[inds, 2]) / 4)
    inclusion_probs <- inclusion_probs / sum(inclusion_probs) # normalize
    inclusion_probs <- inclusion_probs * I_n*2 # Scale to get the desired number of individuals
    
    inclusion_probs[inclusion_probs > 1] <- 1
    sampled_in_strata <- inds[rbinom(length(inds), 1, inclusion_probs) == 1] 
    # sample individuals based on inclusion probabilities
    
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
    X = input_X_des[final_sample,2],
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
