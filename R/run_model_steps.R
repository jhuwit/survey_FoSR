rm(list = ls())
## fit models
library(tidyverse)
library(fastFMM)
library(haven)
library(dplyr)
library(survey)
library(progress)
library(lme4)
library(mgcv)
library(ggplot2)
library(furrr)
library(future)
library(gridExtra)
library(refund)
steps_data = read_rds(here::here("data", "steps_covariates.rds"))
source(here::here("R", "sim_functions_unwt.R"))
source(here::here("R", "utils.R"))

family = "poisson"
B = 100
len = 1440

steps_data =
  steps_data %>%
  rename_with(.cols = starts_with("min"), ~str_replace(., "min_", "Y"))

steps_data =
  steps_data %>%
  mutate(sex_bin = if_else(gender == "Female", 1, 0)) %>%
  mutate(age_cat = cut(age_in_years_at_screening, breaks=c(18, 30, 40, 50, 60, 70, 80), include.lowest = TRUE))

data = steps_data
model_formula = as.formula(paste0('Y~', 'age_cat + sex_bin'))
out_index <- grep(paste0("^", model_formula[2]), names(data)) # indices that start with the outcome name
if (length(out_index) != 1) {
  # functional observations stored in multiple columns
  L <- length(out_index)
} else{
  # observations stored as a matrix in one column using the I() function
  L <- ncol(data[, out_index])
}
argvals <- 1:L

unimm <- function(l) {
  data$Yl <- unclass(data[, out_index][, l]) %>% unlist() %>% unname() %>% round(., digits = 0)
  fit_uni <- suppressMessages(
    glm(
      Yl ~ age_cat + sex_bin,
      # formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
      family = family,
      control = glm.control(maxit = 5000),
      data = data
    )
  )
  betaTilde <- coef(fit_uni)
  return(list(betaTilde = betaTilde))
}

massmm <- lapply(argvals, unimm)


# Obtain betaTilde, fixed effects estimates
betaTilde <- t(do.call(rbind, lapply(massmm, '[[', 1)))
colnames(betaTilde) <- argvals

nknots <- round(L / 2)

betaHat <- t(apply(betaTilde, 1, function(x) mgcv::gam(x ~ s(argvals, bs = "tp",
                                                             k = (nknots + 1)),
                                                       method = "GCV.Cp")$fitted.values))

rownames(betaHat) <- rownames(betaTilde)
colnames(betaHat) <- 1:L

num_boots = 100
# array with dimensions num_coefficients x length functional domain x num_boots
betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), num_boots))

# initiliaze progress bar
pb <- progress_bar$new(format = "  bootstrapping across L [:bar] :percent eta: :eta",,
                       total = L)
# get bootstrap estimates at each point along functional domain

library(future.apply)
for(l in 1:L){
  pb$tick()
  data$Yl <- unclass(data[, out_index][, l]) %>% unlist() %>% unname() %>% round(., digits = 0)
  if(set_seed){
    set.seed(seed)
  }
  plan(multisession, workers = parallel::detectCores() - 1)  # Adjust workers if needed

  coefs <- do.call(
    cbind,
    future_lapply(1:num_boots, function(num) {
      set.seed(l * num)
      index <- sample(row_number(data), dim(data)[1], replace = TRUE)
      dat.tmp <- data[index, ]

      model <- glm(
        formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
        data = dat.tmp,
        family = family)
      model$coefficients
    }, future.seed = TRUE))
  betaTilde_boot[,l,] <- coefs
}
plan(sequential)
## start here

B = ncol(betaTilde_boot[1, , ])
nknots <- round(L / 2)
nknots_cov <- 35
nknots_fpca <- min(round(L / 2), 35)
betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), ncol(betaTilde_boot[1, , ])))
betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
# smooth


for (b in 1:B) {
  betaHat_boot[, , b] <- t(apply(betaTilde_boot[, , b], 1, function(x)
    gam(x ~ s(
      argvals, bs = "tp", k = (2 +  1)
    ), method = "GCV.Cp")$fitted.values))
}
smooth_for_variance = TRUE
for (r in 1:nrow(betaHat)) {
  if (smooth_for_variance) {
    betaHat.var[, , r] <- 1.2 * var(t(betaHat_boot[r, , ]))
  } else{
    betaHat.var[, , r] <- 1.2 * var(t(betaTilde_boot[r, , ]))
  }
}
smooth_for_ci = TRUE
N <- 10000
qn <- rep(0, length = nrow(betaHat))
set.seed(456)
for (i in 1:length(qn)) {
  if (smooth_for_ci) {
    est_bs <- t(betaHat_boot[i, , ])
  } else{
    est_bs <- t(betaTilde_boot[i, , ])
  }
  fit_fpca <- suppressWarnings(refund::fpca.face(est_bs, knots = nknots_fpca)) # suppress sqrt(Eigen$values) NaNs
  ## extract estimated eigenfunctions/eigenvalues
  phi <- fit_fpca$efunctions
  lambda <- fit_fpca$evalues
  K <- length(fit_fpca$evalues)

  ## simulate random coefficients
  theta <- matrix(stats::rnorm(N * K), nrow = N, ncol = K) # generate independent standard normals
  if (K == 1) {
    theta <- theta * sqrt(lambda) # scale to have appropriate variance
    X_new <- tcrossprod(theta, phi) # simulate new functions
  } else{
    theta <- theta %*% diag(sqrt(lambda)) # scale to have appropriate variance
    X_new <- tcrossprod(theta, phi) # simulate new functions
  }
  x_sample <- X_new + t(fit_fpca$mu %o% rep(1, N)) # add back in the mean function
  Sigma_sd <- Rfast::colVars(x_sample, std = TRUE, na.rm = FALSE) # standard deviation: apply(x_sample, 2, sd)
  x_mean <- colMeans(est_bs)
  un <- rep(NA, N)
  for (j in 1:N) {
    un[j] <- max(abs((x_sample[j, ] - x_mean) / Sigma_sd))
  }
  qn[i] <- stats::quantile(un, 0.95)
}

## plot!


Y_mat =
  steps_data %>%
  select(starts_with("Y")) %>%
  mutate(across(.cols = everything(), round)) %>%
  as.matrix()

pffr_df =
  steps_data %>%
  select(-starts_with("Y")) %>%
  mutate(Y_mat = I(Y_mat))  %>%
  mutate(SEQN = factor(SEQN))
model =
  pffr(
    Y_mat ~ sex_bin + age_cat,
    data = pffr_df,
    algorithm = "bam",
    discrete = TRUE,
    family = "poisson",
    method = "fREML",
    bs.yindex = list(bs = "ps", k = 15, m = c(2, 1))
  )

sim_res = list(betaHat = betaHat,
               betaHat.var = betaHat.var,
               qn = qn,
               pffr_mod = model)

save(sim_res, here::here("data", "steps_model.Rda"))
get_plot_df = function(sim_res) {
  num_var = nrow(sim_res$betaHat)
  map_dfr(
    .x = 1:num_var,
    .f = function(r) {
      data.frame(
        s = 1:length(sim_res$betaHat[r, ]),
        beta = sim_res$betaHat[r, ],
        lower = sim_res$betaHat[r, ] - 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
        upper = sim_res$betaHat[r, ] + 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
        lower.joint = sim_res$betaHat[r, ] - sim_res$qn[r] *
          sqrt(diag(sim_res$betaHat.var[, , r])),
        upper.joint = sim_res$betaHat[r, ] + sim_res$qn[r] * sqrt(diag(sim_res$betaHat.var[, , r]))
      ) %>%
        mutate(name = sim_res$betaHat %>% rownames() %>% .[r])
    }
  )
}
plt_df_boot = get_plot_df(sim_res)


plt_df_boot %>%
  ggplot(aes(x = s, y = beta)) +
  geom_line()  +
  facet_wrap(.~name, scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper), alpha = .3,
              fill = "gray10") +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint), alpha = .3,
              fill = "gray20") +
  labs(x = "L", y = "Estimate") +
  theme_classic() +
  geom_hline(aes(yintercept = 0))



# compare to pffr

Y_mat =
  steps_data %>%
  select(starts_with("Y")) %>%
  mutate(across(.cols = everything(), round)) %>%
  as.matrix()

pffr_df =
  steps_data %>%
  select(-starts_with("Y")) %>%
  mutate(Y_mat = I(Y_mat))  %>%
  mutate(SEQN = factor(SEQN))
model =
  pffr(
    Y_mat ~ sex_bin + age_cat,
    data = pffr_df,
    algorithm = "bam",
    discrete = TRUE,
    family = "poisson",
    method = "fREML",
    bs.yindex = list(bs = "ps", k = 15, m = c(2, 1))
  )


plot(model)





