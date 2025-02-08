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

# Generate Population Data
I <- 10000 ## number of subjects
L <- 50 ## dimension of the functional domain
SNR_sigma = 1
family <- "gaussian" ## family of longitudinal functional data

## specify true fixed effects functions
grid <- seq(0, 1, length = L)
beta_true <- matrix(NA, 2, L) # time-dependent beta-true trajectory
scenario <- 1 ## indicate the scenario to use
if(scenario == 1){
  beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
  beta_true[2,] = dnorm(grid, .6, .15)/20
}else if(scenario == 2){
  beta_true[1,] <- 0.53 + 0.06*sin(3*grid*pi) - 0.03*cos(6.5*grid*pi)
  beta_true[2,] <- dnorm(grid, 0.2, .1)/60 + dnorm(grid, 0.35, .1)/200 -
    dnorm(grid, 0.65, .06)/250 + dnorm(grid, 1, .07)/60
}
rownames(beta_true) <- c("Intercept", "x")

# do we need to set seed anywhere?

simulate_survey_fosr <- function(I = 10e6, # total population,
                                 I_n = 50, # subjects in each psu-strata combination
                                 alpha = 2, # number PSU
                                 L = 50, # dimension of the functional domain
                                 SNR_sigma = 1,
                                 num_strata = 30 # total strata
){
  # Generate X variable
  X_des = cbind(1, rnorm(I, 0, 2))  # design matrix: Intercept + predicted variables

  # Generate strata sizes using a Dirichlet distribution
  dirichlet_probs = gtools::rdirichlet(1, rep(1, num_strata))
  stratum_assignments = sample(1:num_strata, size = I, replace = TRUE, prob = dirichlet_probs)
  p_strata = table(stratum_assignments) / sum(table(stratum_assignments))
  p_strata = c(p_strata) %>% unname()
  w_strata = 1 / p_strata


  # Second-stag: Individual Sampling within each strata
  final_sample <- c()
  p_i_given_g <- rep(0, I)
  psu = c()

  for (strata in 1:num_strata) {
    inds <- which(stratum_assignments == strata)
    inclusion_probs <- exp(X_des[inds, 2])  # Compute individual inclusion probabilities based on covariate
    inclusion_probs <- inclusion_probs / sum(inclusion_probs) # normalize
    inclusion_probs <- inclusion_probs * I_n * 2 # Scale to get the desired number of individuals in each strata / PSU combo
    inclusion_probs[inclusion_probs > 1] <- 1
    sampled_in_strata <- inds[rbinom(length(inds), 1, inclusion_probs) == 1]
    final_sample <- c(final_sample, sampled_in_strata)
    p_i_given_g[inds] <- inclusion_probs
    psu_strata <- rbinom(length(sampled_in_strata), 1, 0.5) + 1
    psu <- c(psu, psu_strata)
  }
  length(final_sample) # 2939, close to 30 * 50 * 2

  # Compute final individual selection probabilities
  p_i <- p_strata[stratum_assignments] * p_i_given_g
  survey_weights <- 1 / p_i  # Final survey weights (Inf: Not selected)

  dat.sim <- data.frame(
    ID = final_sample,
    X = X_des[final_sample,2],
    strata = stratum_assignments[final_sample], # Assume same strata here ??
    psu = psu,
    weight = survey_weights[final_sample]
  )

  # Generate fixed effects and true Y.
  true_fixef <- X_des %*% beta_true
  true_fixef_sample <- true_fixef[final_sample,]

  Y_obs <- matrix(NA, nrow(true_fixef_sample), L)
  for(i in 1:nrow(true_fixef_sample)){
    for(j in 1:L){
      Y_obs[i, j] <- rnorm(1, mean = true_fixef_sample[i, j], sd = SNR_sigma / 5)
    }
  }
  # Y_obs <- Y_obs + with(dat.sim, rnorm(alpha, sd = SNR_sigma / 10)[psu]) # not sure if we want to do this part?

  # Y_obs2 <- matrix(rnorm(nrow(true_fixef_sample) * L,
  #                       mean = as.vector(true_fixef_sample),
  #                       sd = SNR_sigma / 5),
  #                 nrow = nrow(true_fixef_sample),
  #                 ncol = L)

  # Y_obs <- true_fixef_sample + with(dat.sim, rnorm(alpha, sd = SNR_sigma)[psu]) ## add some PSU-level variability in Y

  Y_obs <- data.frame(Y_obs)
  colnames(Y_obs) <- paste0("Y", 1:L)

  sample_data <- cbind(dat.sim, Y_obs)

  return(sample_data)
}
# dat = replicate(100, simulate_survey_fosr(), simplify = FALSE)
# write_rds(dat, here::here("sim_data", "dat_f1.rds"), compress = "xz")

# data = simulate_survey_fosr()
# write_rds(dat, here::here("sim_data", "dat_f2.rds"), compress = "xz")
