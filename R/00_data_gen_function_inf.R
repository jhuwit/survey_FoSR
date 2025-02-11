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

simulate_survey_fosr_informative <- function(X_des, # design matrix
                                 Y_fpca, # fpca object
                                 beta_true, # true beta
                                 true_fixef, # true fixed effects (Y)
                                 Y_obs, # y matrix
                                 I_n = 100, # subjects in each psu-strata combination
                                 alpha = 2, # number PSU
                                 L = 50, # dimension of the functional domain
                                 num_strata = 30,  # total strata
                                 sampled_strata = 30
){
  score1 = fpca_Y$scores[,1]
  score2 = fpca_Y$scores[,2]
  score3 = fpca_Y$scores[,3]
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
    # inclusion_probs <- exp(X_des[inds, 2])  # Compute individual inclusion probabilities based on covariate
    # inclusion_probs <- exp(-10 * score1 + 10 * score2 + 10 * score3)
    inclusion_probs <- exp(10 * score1[inds] + 10 * score2[inds] + 10 * score3[inds])
    # high = which(inclusion_probs >= quantile(inclusion_probs, .999))
    # low = which(inclusion_probs <= quantile(inclusion_probs, .001))
    #
    #
    # sample = Y_obs %>% as.data.frame() %>%
    #   ungroup() %>%
    #   mutate(rn = row_number(),
    #          type = if_else(rn %in% high, "high", "low")) %>%
    #   filter(rn %in% c(high, low))
    #
    # Y_mat =
    #   sample %>% select(starts_with("V")) %>% as.matrix()
    #
    # sample %>%
    #   mutate(Y = matrix(Y_mat, ncol = 50),
    #        Y_tf = tfd(Y, length = 50),
    #        Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
    #   group_by(type) %>%
    #   mutate(mean_curve = mean(Y_smooth)) %>%
    #   ungroup() %>%
    #   ggplot(aes(y = mean_curve, color = type)) +
    #   geom_spaghetti(lwd = 1.1) +
    #   # geom_spaghetti(lwd = .7) +
    #   labs(x = "L", y = "Y", title = "Mean curves and subjects by PSU")
    #

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
  # length(final_sample) # 2939, close to 30 * 50 * 2

  # Compute final individual selection probabilities
  p_i <- p_strata[stratum_assignments] * p_i_given_strata
  survey_weights <- 1 / p_i  # Final survey weights (Inf: Not selected)

  dat.sim <- data.frame(
    ID = final_sample,
    X = X_des[final_sample,2],
    strata = stratum_assignments[final_sample], # Assume same strata here ??
    psu = psu,
    weight = survey_weights[final_sample]
  )

  # Generate fixed effects and true Y.
  true_fixef_sample <- true_fixef[final_sample,]


  Y_sample <- data.frame(Y_obs[final_sample,])
  colnames(Y_sample) <- paste0("Y", 1:L)

  sample_data <- cbind(dat.sim, Y_sample)
  sample_data <- sample_data %>%
    filter(!is.na(weight)) # shouldn't happen but just in case

  return(sample_data)
}

