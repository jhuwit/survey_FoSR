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
L <- 30 ## dimension of the functional domain
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


set.seed(209)
options(survey.lonely.psu = "adjust")

simulate_survey_fosr <- function(I = 10000, # total subjects
                                 L = 30, # dimension of the functional domain
                                 G = 10, # total cluster number
                                 Gs = 5 # total selected cluster number
                                 ){
  # Generate X variable
  X_des = cbind(1, rnorm(I, 0, 2))  # design matrix: Intercept + predicted variables

  # Generate cluster sizes using a Dirichlet distribution

  cluster_sizes <- rdirichlet(1, rep(1, G)) * I  # Normalize to sum to I
  cluster_sizes <- round(cluster_sizes)  # Round to integers
  cluster_sizes[length(cluster_sizes)] <- I - sum(cluster_sizes[-length(cluster_sizes)])
  # Adjust to ensure total population size matches I
  cluster_assignment <- rep(1:G, cluster_sizes) # Assign individuals to clusters

  # First-stage: Sample Gs cluster using PPS
  sampled_clusters <- sample(1:G, Gs, replace = FALSE, prob = cluster_sizes)
  sampled_individuals <- which(cluster_assignment %in% sampled_clusters)
  ## calculate p1 and w1
  p_g <- cluster_sizes / sum(cluster_sizes)
  w_g <- 1/p_g

  # Second-stag: Individual Sampling within selected cluster
  final_sample <- c()
  p_i_given_g <- rep(0, I)

  for (g in sampled_clusters) {
    cluster_inds <- which(cluster_assignment == g)
    scaling_factor <- length(cluster_inds)/10 # TBD

    inclusion_probs <- exp(X_des[cluster_inds, 2])  # Compute individual inclusion probabilities
    inclusion_probs <- inclusion_probs / sum(inclusion_probs) * scaling_factor# Normalize probabilities
    inclusion_probs[inclusion_probs > 1] <- 1

    sampled_in_cluster <- cluster_inds[rbinom(length(cluster_inds), 1, inclusion_probs) == 1]
    final_sample <- c(final_sample, sampled_in_cluster)
    p_i_given_g[cluster_inds] <- inclusion_probs
  }
  length(final_sample) # N = 482

  # Compute final individual selection probabilities
  p_i <- p_g[cluster_assignment] * p_i_given_g
  survey_weights <- 1 / p_i  # Final survey weights (Inf: Not selected)

  dat.sim <- data.frame(
    ID = final_sample,
    X = X_des[final_sample,2],
    strata = 1, # Assume same strata here
    psu = cluster_assignment[final_sample],
    weight = survey_weights[final_sample]
  )

  # Generate fixed effects and true Y.
  true_fixef <- X_des %*% beta_true
  pop = data.frame(id = 1:10000,
                   g = cluster_assignment, # cluster
                   X = X_des[,2]
  )
  Y_obs <- true_fixef + with(pop, rnorm(10, s = 1)[g]) # PSU level variability in Y
  Y_obs <- data.frame(Y_obs)
  colnames(Y_obs) <- paste0("Y", 1:L)

  sample_data <- cbind(dat.sim, Y_obs[final_sample,])

  return(sample_data)
}



source('R/wFUI.R')
################################################################################
## Do simulations on a local laptop (results in the manuscript are from JHPCE)
################################################################################
nsim <- 5
results <- list()
boot_types = 'Rao-Wu-Yue-Beaumont' # c('ref_design', 'Rao-Wu', 'ref_noDesign')
# 'Rao-Wu-Yue-Beaumont', 'Preston'
sim_local <- list() ## store simulation results
for (boot_type in boot_types){
  for(iter in 1:nsim){

    ################################################################################
    ## Generate simulated data
    ################################################################################
    set.seed(iter)
    data <- simulate_survey_fosr()

    ## pre-process simulated data
    Y_mat <- data[, grepl('Y.', colnames(data))]
    dat.fit <- data.frame(Y = Y_mat,
                          data[,!grepl('Y.', colnames(data))])
    formula.FUI = as.formula(paste0('Y~', 'X'))

    ################################################################################
    ## Implement different estimation methods
    ################################################################################
    ptm <- proc.time()
    mod.Output <- fui.weight(formula = as.formula('Y ~ X'),
                 analytic = FALSE,
                 argvals = seq(from = 1, to = dim(beta_true)[2], by = 1),
                 weighted = TRUE,
                 random = FALSE,
                 boot_type = boot_type,
                 num_boots = 100,
                 nknots_min = 5,
                 data = dat.fit)

    lfosr3stime <- (proc.time() - ptm)[3]

    ################################################################################
    ## Organize simulation results
    ################################################################################

    ## lfosr3s
    MISE_lfosr3s <- rowMeans((mod.Output$betaHat - beta_true)^2)
    cover_lfosr3s <- rep(NA, nrow(beta_true))
    cover_lfosr3s_pw <- matrix(FALSE, nrow(beta_true), L)
    for(p in 1:length(cover_lfosr3s)){
      cover_upper <- which(mod.Output$betaHat[p,]+mod.Output$qn[p]*sqrt(diag(mod.Output$betaHat.var[,,p])) > beta_true[p,])
      cover_lower <- which(mod.Output$betaHat[p,]-mod.Output$qn[p]*sqrt(diag(mod.Output$betaHat.var[,,p])) < beta_true[p,])
      cover_lfosr3s[p] <- length(intersect(cover_lower, cover_upper)) == L
      cover_upper_pw <- which(mod.Output$betaHat[p,]+1.96*sqrt(diag(mod.Output$betaHat.var[,,p])) > beta_true[p,])
      cover_lower_pw <- which(mod.Output$betaHat[p,]-1.96*sqrt(diag(mod.Output$betaHat.var[,,p])) < beta_true[p,])
      cover_lfosr3s_pw[p,intersect(cover_lower_pw, cover_upper_pw)] <- TRUE
    }

    ## results of a single simulation
    sim_result <- list(MISE_lfosr3s = MISE_lfosr3s,
                       cover_lfosr3s = cover_lfosr3s,
                       cover_lfosr3s_pw = cover_lfosr3s_pw,
                       time_lfosr3s = lfosr3stime,
                       I = I, L = L, SNR_sigma = SNR_sigma,
                       family = family, scenario = scenario)

    sim_local[[iter]] <- sim_result
    print(iter)
  }

  ################################################################################
  ## Obtain MISE, coverage, and computing time
  ################################################################################
  ## MISE
  MISE_lfosr3s <- lapply(sim_local, '[[', 1) %>% bind_rows()
  colMeans(MISE_lfosr3s)

  ## coverage
  ### joint
  cover_lfosr3s <- t(lapply(sim_local, '[[', 2) %>% bind_cols())
  ### pointwise
  cover_lfosr3s_pw <- array(NA, dim = c(nsim, L, p))
  for(i in 1:nsim){
    cover_lfosr3s_pw[i,,1] <- sim_local[[i]]$cover_lfosr3s_pw[2,] # intercept
    cover_lfosr3s_pw[i,,2] <- sim_local[[i]]$cover_lfosr3s_pw[2,] # fixed effects
  }
  ### print results
  print(boot_type)
  print(paste0("lfosr3s coverage (joint CI): ", colMeans(cover_lfosr3s)[2]))
  print(paste0("lfosr3s coverage (pointwise CI): ", apply(cover_lfosr3s_pw, 3, mean)[2]))

  ## time
  time_lfosr3s <- lapply(sim_local, '[[', 4) %>% bind_rows()
  # print(apply(time_lfosr3s, 2, median))

  results[[boot_type]] = data.frame(cover_ratio_jw = colMeans(cover_lfosr3s)[2],
                                    cover_ratio_pw = apply(cover_lfosr3s_pw, 3, mean)[2],
                                    time =  apply(time_lfosr3s, 2, median))
}

#
# ## Clean results
# ################################################################################
# ## Obtain MISE, coverage, and computing time
# ################################################################################
# ## MISE
# MISE_lfosr3s <- lapply(sim_local, '[[', 1) %>% bind_rows()
# colMeans(MISE_lfosr3s)
#
# ## coverage
# ### joint
# cover_lfosr3s <- t(lapply(sim_local, '[[', 2) %>% bind_cols())
# ### point-wise
# cover_lfosr3s_pw <- array(NA, dim = c(nsim, L, p))
# for(i in 1:nsim){
#   cover_lfosr3s_pw[i,,1] <- sim_local[[i]]$cover_lfosr3s_pw[2,] # intercept
#   cover_lfosr3s_pw[i,,2] <- sim_local[[i]]$cover_lfosr3s_pw[2,] # fixed effects
# }
# ### print results
# print(paste0("lfosr3s coverage (joint CI): ", colMeans(cover_lfosr3s)[2]))
# print(paste0("lfosr3s coverage (pointwise CI): ", apply(cover_lfosr3s_pw, 3, mean)[2]))
#
# ## time
# time_lfosr3s <- lapply(sim_local, '[[', 4) %>% bind_rows()
# print(apply(time_lfosr3s, 2, median))
