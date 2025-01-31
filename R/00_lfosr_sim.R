#' simulate functional response Y (matrix) and a scalar predictor X (vector)
#' single level data
library(tidyverse)
library(mvtnorm)
library(tidyfun)

# generate total population 10^7 ppl
lfosrsim_sl <- function(family = "gaussian",
                        I = 500, # n subs
                        L = 100, # number of grid points
                        SNR_sigma = 1,
                        psu_beta = c(0 , 0)){ # signal-to-noise ratio in Gaussian data generation

  # I <- 50 ## number of subjects
  # L <- 50 ## dimension of the functional domain
  # SNR_sigma <- 1 ## signal-to-noise ratio
  # family <- "gaussian" ## family of longitudinal functional data

  ## specify true fixed effects functions
  ## generate fixed effects
  # n <- sum(J_subj)
  n <- I
  psus <- rbinom(n, 1, 0.5) + 1 # PSUs are either 1 or 2 with equal probability
  # design matrix is 1, N(0, 4), and PSUs
  X_des = cbind(1, rnorm(n, 0, 2), if_else(psus == 1, psu_beta[1], psu_beta[2]))



  grid <- seq(0, 1, length = L)
  beta_true <- matrix(NA, 3, L)
  beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
  beta_true[2,] =  dnorm(grid, .6, .15)/20
  beta_true[3,] = 1
  rownames(beta_true) <- c("Intercept", "x1", "x2")


  fixef <- X_des %*% beta_true # fixed effect is xbeta

  fixef %>%
    as.data.frame() %>%
    mutate(id = row_number()) %>%
    filter(id <= 10) %>%
    mutate(psu = psus[1:10]) %>%
    pivot_longer(cols = -c(id, psu)) %>%
    mutate(x = as.numeric(sub(".*V", "", name))) %>%
    ggplot(aes(x = x, y = value, color = factor(psu), group = as.factor(id))) +
    geom_line()

  eta_true <- fixef
  ## generate longitudinal functional data
  Y_obs <- matrix(NA, n, L)
  Y_obs_psu <- matrix(NA, n, L)
  p_true <- plogis(eta_true)
  lam_true <- exp(eta_true)
  sd_signal <- sd(eta_true)

  ## trying to induce assocation between PSU and outcome
  for(i in 1:n){
    for(j in 1:L){
      if(family == "gaussian"){
        Y_obs[i, j] <- rnorm(1, mean = eta_true[i, j], sd = sd_signal/SNR_sigma)
        # try to induce dependence
      }else if(family == "binomial"){
        Y_obs[i, j] <- rbinom(1, 1, p_true[i, j])

      }else if(family == "poisson"){
        Y_obs[i, j] <- rpois(1, lam_true[i, j])
      }
    }
  }

  ## combine simulated data

  dat.sim <- data.frame(ID = 1:n, X1 = X_des[,2], PSU = psus, Y = I(Y_obs), Y_psu = I(Y_obs_psu), eta = I(eta_true))

  return(dat.sim)
}
I = 10^6 ## number of subjects
L = 50 ## dimension of the functional domain
SNR_sigma = 1 ## signal-to-noise ratio
family = "gaussian" ## family of longitudinal functional data

data_sim =
  lfosrsim_sl(family, I, L, SNR_sigma)

write_rds(data_sim, here::here("data", "sim", "superpop_nopsu.rds"), compress = "gz")

data_sim_psu = lfosrsim_sl(family, I, L, SNR_sigma, c(1.1, 1.2))

write_rds(data_sim_psu, here::here("data", "sim", "superpop_psu.rds"), compress = "gz")

### plot the simulated data

test = data_sim %>%
  slice(1:10^5) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50))

test %>%
  mutate(PSU = factor(PSU),
        y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(PSU) %>%
  mutate(mean_curve = mean(y_smooth)) %>%
  filter(ID <= 100) %>%
  ggplot(aes(y = y_smooth, color = PSU, group = ID)) +
  geom_spaghetti(lwd = .7, alpha = .1) + # Individual curves
  geom_spaghetti(aes(y = mean_curve), lwd = 1.2)

test_psu = data_sim_psu %>%
  slice(1:10^5) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50))

test_psu %>%
  mutate(PSU = factor(PSU),
         y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(PSU) %>%
  mutate(mean_curve = mean(y_smooth)) %>%
  filter(ID <= 100) %>%
  ggplot(aes(y = y_smooth, color = PSU, group = ID)) +
  geom_spaghetti(lwd = .7, alpha = .1) + # Individual curves
  geom_spaghetti(aes(y = mean_curve), lwd = 1.2)

## confirmed that it seems to work where we get different means by PSU
