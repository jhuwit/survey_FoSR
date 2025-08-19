source(here::here("R", "00_data_gen_function.R"))

L = 50
grid = seq(0, 1, length = L)
beta_true = matrix(NA, 2, L) # time-dependent beta-true trajectory
scenario = 1 ## indicate the scenario to use
if(scenario == 1){
  beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
  beta_true[2,] = dnorm(grid, .6, .15)/20
}else if(scenario == 2){
  beta_true[1,] = 0.53 + 0.06*sin(3*grid*pi) - 0.03*cos(6.5*grid*pi)
  beta_true[2,] = dnorm(grid, 0.2, .1)/60 + dnorm(grid, 0.35, .1)/200 -
    dnorm(grid, 0.65, .06)/250 + dnorm(grid, 1, .07)/60
}
rownames(beta_true) <- c("Intercept", "x")

I = 10e6 # superpopulation size
set.seed(4574)
X_des = cbind(1, rnorm(I, 0, 2)) # design matrix: Intercept + predicted variables
true_fixef <- X_des %*% beta_true # true fixed effect
SNR_sigma = 1 # signal to noise ratio

set.seed(4577) # generate Y obs for all data from normal rv
Y_obs <- matrix(
  rnorm(n = nrow(true_fixef) * L, mean = as.vector(t(true_fixef)), sd = SNR_sigma),
  nrow = nrow(true_fixef),
  ncol = L,
  byrow = TRUE
)

# save objects
write_rds(X_des, here::here("sim_data", "X_des.rds"), compress = "xz")
write_rds(Y_obs, here::here("sim_data", "Y_obs.rds"), compress = "xz")
write_rds(true_fixef, here::here("sim_data", "X_beta.rds"), compress = "xz")

scenario <- 2 ## indicate the scenario to use
if(scenario == 1){
  beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
  beta_true[2,] = dnorm(grid, .6, .15)/20
}else if(scenario == 2){
  beta_true[1,] <- 0.53 + 0.06*sin(3*grid*pi) - 0.03*cos(6.5*grid*pi)
  beta_true[2,] <- dnorm(grid, 0.2, .1)/60 + dnorm(grid, 0.35, .1)/200 -
    dnorm(grid, 0.65, .06)/250 + dnorm(grid, 1, .07)/60
}
rownames(beta_true) <- c("Intercept", "x")
true_fixef <- X_des %*% beta_true # true fixed effect

set.seed(4577) # generate Y obs for all data from normal rv
Y_obs <- matrix(
  rnorm(n = nrow(true_fixef) * L, mean = as.vector(t(true_fixef)), sd = SNR_sigma),
  nrow = nrow(true_fixef),
  ncol = L,
  byrow = TRUE
)
write_rds(Y_obs, here::here("sim_data", "Y_obs_s2.rds"), compress = "xz")
write_rds(true_fixef, here::here("sim_data", "X_beta_s2.rds"), compress = "xz")

######### stop here

set.seed(2)
d2 =  simulate_survey_fosr(X_des = X_des,
                           beta_true  = beta_true,
                           true_fixef = true_fixef,
                           L = 50,
                           I_n = 100,
                           alpha = 2,
                           num_strata = 30,
                           sampled_strata = 30,
                           Y_obs = Y_obs)

Y_mat = d2 %>% select(starts_with("Y")) %>%
  as.matrix()


d2 %>%
  mutate(psu = paste0("PSU: ", psu)) %>%
  select(-starts_with("Y")) %>%
  mutate(Y = I(Y_mat)) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50),
         Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  mutate(mean_curve = mean(Y_smooth)) %>%
  ungroup() %>%
  slice(1:10) %>%
  ggplot(aes(y = Y_tf)) +
  geom_spaghetti(aes(y = mean_curve), lwd = 1.1) +
  geom_spaghetti(lwd = .7, aes(color = factor(ID))) +
  theme(legend.position = "none") +
  labs(x = "L", y = "Y", title = "Mean curves for 10 subjects")


d2 %>%
  select(-starts_with("Y")) %>%
  mutate(Y = I(Y_mat)) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50),
         Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  mutate(mean_curve = mean(Y_smooth)) %>%
  ungroup() %>%
  slice(1:10) %>%
  ggplot(aes(y = Y_smooth)) +
  geom_spaghetti(aes(y = mean_curve), lwd = 1.1) +
  geom_spaghetti(lwd = .7, aes(color = factor(ID))) +
  theme(legend.position = "none") +
  labs(x = "L", y = "Y (smoothed)", title = "Mean curves and subjects by PSU")

# get betatilde/hat and plot

dat.fit <- data.frame(Y = Y_mat, d2[, !grepl('Y.', colnames(d2))])
formula.FUI = as.formula(paste0('Y~', 'X'))

weighted = TRUE
family = "gaussian"
unimm <- function(l){
  data$Yl <- unclass(data[,out_index][,l])

  svy_design <- svydesign(ids = ~ psu,
                          strata = ~ strata,
                          weights = ~ weight,
                          data = data,
                          nest = TRUE)
  # fit survey-weighted regression using svyglm
  if(weighted == TRUE){
    if(family == "gaussian"){
      fit_uni <- suppressMessages(svyglm(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                         design = svy_design,
                                         control = glm.control(maxit = 5000)))
    }else{
      fit_uni <- suppressMessages(svyglm(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                         design = svy_design,
                                         family = family,
                                         control = glm.control(maxit = 5000)))
    }} else{
      if(family == "gaussian"){
        fit_uni <- suppressMessages(lm(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                       data = data))
      }else{
        fit_uni <- suppressMessages(glm(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                        data = data,
                                        family = family) )
      }
    }

  ## fixed effects estimates
  if (weighted == TRUE) { ## this part seems redundant
    betaTilde <- coef(fit_uni)
  }else{
    betaTilde <- coef(fit_uni)
  }

  return(list(betaTilde = betaTilde))

}

model_formula <- formula.FUI
out_index <- grep(paste0("^", model_formula[2]), names(d2)) # indices that start with the outcome name
if(length(out_index) != 1){ # functional observations stored in multiple columns
  L <- length(out_index)
}else{ # observations stored as a matrix in one column using the I() function
  L <- ncol(data[,out_index])
}
argvals <- 1:L


massmm <- lapply(argvals, unimm)

# Obtain betaTilde, fixed effects estimates
betaTilde <- t(do.call(rbind, lapply(massmm, '[[', 1)))
colnames(betaTilde) <- argvals

nknots_min = NULL
nknots_min_cov = 35
nknots <- min(round(L / 2), nknots_min)
nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
nknots_fpca <- min(round(L / 2), 35)

betaHat <- t(apply(betaTilde, 1, function(x) mgcv::gam(x ~ s(argvals, bs = "tp", k = (nknots + 1)), method = "GCV.Cp")$fitted.values))

rownames(betaHat) <- rownames(betaTilde)
colnames(betaHat) <- 1:L



betaHat <- t(apply(betaTilde, 1, function(x) mgcv::gam(x ~ s(argvals, bs = "tp", k = (nknots + 1)), method = "GCV.Cp")$fitted.values))

rownames(betaHat) <- rownames(betaTilde)
colnames(betaHat) <- 1:L


bt_plot =
  betaHat %>%
  as_tibble() %>%
  mutate(beta = c("beta_0", "beta_1")) %>%
  pivot_longer(cols = -beta) %>%
  mutate(x = as.numeric(sub(".*V", "", name)))
btrue_plot =
  beta_true %>%
  as_tibble() %>%
  mutate(beta = c("beta_0", "beta_1")) %>%
  pivot_longer(cols = -beta) %>%
  mutate(x = as.numeric(sub(".*V", "", name)))

bt_plot %>% mutate(type = "estimate") %>%
  bind_rows(btrue_plot %>% mutate(type = "true")) %>%
  ggplot(aes(x = x, y = value, color = type)) +
  facet_wrap(.~beta) +
  geom_line() +
  labs(x = "L", y = "Value", title = "Coefficient Functions")
