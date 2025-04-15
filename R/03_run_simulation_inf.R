library(future)
library(furrr)
library(tidyverse)
force = TRUE
library(here)
library(devtools)
library(fastFMM)
library(haven)
library(dplyr)
library(survey)
library(progress)
library(lme4)
library(paletteer)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(tidyfun)
library(mvtnorm)
library(refund)
library(svrep)

source(here::here("R", "01_sim_functions.R"))
get_betatilde = function(data, family = "gaussian", type = "design",
                         model_formula = as.formula(paste0('Y~', 'X'))){
  if(!(type %in% c("design", "weighted", "none"))){
    stop("please specify boot_type as one of 'design', 'weighted', 'none'")
  }

  out_index <- grep(paste0("^", model_formula[2]), names(data)) # indices that start with the outcome name
  if(length(out_index) != 1){ # functional observations stored in multiple columns
    L <- length(out_index)
  }else{ # observations stored as a matrix in one column using the I() function
    L <- ncol(data[,out_index])
  }
  argvals <- 1:L
  unimm_design <- function(l) {
    data$Yl <- unclass(data[, out_index][, l])
    svy_design <- svydesign(
      ids = ~ psu,
      strata = ~ strata,
      weights = ~ weight,
      data = data,
      nest = TRUE
    )
    # fit survey-weighted regression using svyglm
    fit_uni <- suppressMessages(
      svyglm(
        formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
        design = svy_design,
        family = family,
        control = glm.control(maxit = 5000)
      )
    )
    betaTilde <- coef(fit_uni)
    return(list(betaTilde = betaTilde))
  }
  unimm_wt <- function(l) {
    data$Yl <- unclass(data[, out_index][, l])

    fit_uni <- suppressMessages(
      glm(
        formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
        family = family,
        control = glm.control(maxit = 5000),
        weights = weight,
        data = data
      )
    )
    betaTilde <- coef(fit_uni)
    return(list(betaTilde = betaTilde))
  }
  unimm <- function(l) {
    data$Yl <- unclass(data[, out_index][, l])
    fit_uni <- suppressMessages(glm(
      formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
      family = family,
      control = glm.control(maxit = 5000),
      data = data
    ))
    betaTilde <- coef(fit_uni)
    return(list(betaTilde = betaTilde))
  }
  if (type == "design") {
    massmm <- lapply(argvals, unimm_design)
  } else if (type == "weighted") {
    massmm <- lapply(argvals, unimm_wt)
  } else if (type == "none") {
    massmm <- lapply(argvals, unimm)
  }

  # Obtain betaTilde, fixed effects estimates
  betaTilde <- t(do.call(rbind, lapply(massmm, '[[', 1)))
  colnames(betaTilde) <- argvals
  return(betaTilde)
}

get_betahat = function(betaTilde, L = 50,
                       nknots_min = NULL){
  argvals = 1:L
  nknots <- min(round(L / 2), nknots_min)

  betaHat <- t(apply(betaTilde, 1, function(x) mgcv::gam(x ~ s(argvals, bs = "tp",
                                                               k = (nknots + 1)),
                                                         method = "GCV.Cp")$fitted.values))

  rownames(betaHat) <- rownames(betaTilde)
  colnames(betaHat) <- 1:L
  return(betaHat)
}

#
generate_superpopulation <- function(I = 10e6,
                                     L = 50,
                                     scenario = 1,
                                     family = "gaussian",
                                     seed = 4574,
                                     num_strata = 30,
                                     psu_per_strata = 2,
                                     SNR_sigma = 0.5,
                                     sigma_strata = 0.2,
                                     sigma_psu = 0.1) {
  stopifnot(scenario %in% c(1, 2))
  stopifnot(family %in% c("gaussian", "poisson", "binomial"))

  grid <- seq(0, 1, length = L)
  beta_fixed <- matrix(NA, 2, L)
  if (scenario == 1) {
    beta_fixed[1, ] <- -0.15 - 0.1 * sin(2 * pi * grid) - 0.1 * cos(2 * pi * grid)
    beta_fixed[2, ] <- dnorm(grid, 0.6, 0.15) / 20
  } else {
    beta_fixed[1, ] <- 0.53 + 0.06 * sin(3 * pi * grid) - 0.03 * cos(6.5 * pi * grid)
    beta_fixed[2, ] <- dnorm(grid, 0.2, 0.1) / 60 + dnorm(grid, 0.35, 0.1) / 200 -
      dnorm(grid, 0.65, 0.06) / 250 + dnorm(grid, 1, 0.07) / 60
  }

  rownames(beta_fixed) <- c("Intercept", "x")

  nbasis <- 5
  basis <- fda::create.bspline.basis(c(0, 1), nbasis)
  Phi <- fda::eval.basis(grid, basis)

  set.seed(seed)
  X_des <- cbind(1, rnorm(I, 0, 2))

  set.seed(seed)
  dirichlet_probs <- gtools::rdirichlet(1, rep(1, num_strata))
  stratum_assignments <- sample(1:num_strata, I, replace = TRUE, prob = dirichlet_probs)

  set.seed(seed)
  psu_scores <- matrix(rnorm(num_strata * psu_per_strata * nbasis, 0, sigma_psu), num_strata * psu_per_strata, nbasis)
  psu_assignments <- rep(NA, I)
  for (s in 1:num_strata) {
    psus_in_stratum <- sample(1:psu_per_strata, sum(stratum_assignments == s), replace = TRUE)
    psu_assignments[stratum_assignments == s] <- paste0(s, "_", psus_in_stratum)
  }

  strata_scores <- matrix(rnorm(num_strata * nbasis, 0, sigma_strata), num_strata, nbasis)
  strata_random_effects <- strata_scores %*% t(Phi)
  psu_random_effects <- psu_scores %*% t(Phi)

  strata_effects_indiv <- strata_random_effects[stratum_assignments, ]
  psu_effects_indiv <- psu_random_effects[as.numeric(factor(psu_assignments)), ]
  random_effects <- strata_effects_indiv + psu_effects_indiv

  ## 🚨 Add stratum-specific slope modifications
  set.seed(seed)
  hetero_slopes <- matrix(rnorm(num_strata * L, mean = 0, sd = 0.05), nrow = num_strata, ncol = L)
  beta1_by_stratum <- hetero_slopes[stratum_assignments, ] + matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)

  # Construct fixed effect signal
  fixef_signal <- matrix(rep(beta_fixed[1, ], I), nrow = I, byrow = TRUE) +
    X_des[, 2] * beta1_by_stratum

  # Generate outcomes
  set.seed(seed)
  if (family == "gaussian") {
    Y_obs <- matrix(
      rnorm(n = I * L, mean = as.vector(t(fixef_signal)) + as.vector(t(random_effects)), sd = SNR_sigma),
      nrow = I,
      ncol = L,
      byrow = TRUE
    )
  } else {
    stop("Only Gaussian implemented for now.")
  }

  return(list(
    Y_obs = Y_obs,
    X_des = X_des,
    stratum_assignments = stratum_assignments,
    psu_assignments = psu_assignments,
    dirichlet_probs = dirichlet_probs,
    beta_true = beta_fixed
  ))
}


simulate_survey_fosr_inf <- function(X_des, Y_obs, fpca_Y,
                                     I_n = 100,
                                     num_strata = 30,
                                     sampled_strata = 30,
                                     stratum_assignments,
                                     psu_assignments,
                                     dirichlet_probs,
                                     L = 50) {
  score1 = fpca_Y$scores[,1]
  score2 = fpca_Y$scores[,2]
  score3 = fpca_Y$scores[,3]
  X1 = X_des[, 2]
  I = nrow(X_des)

  if (sampled_strata == num_strata) {
    selected_strata = 1:num_strata
  } else {
    selected_strata = sample(1:num_strata, sampled_strata, replace = FALSE, prob = dirichlet_probs)
  }

  p_strata_design <- dirichlet_probs
  p_strata_design[!(1:num_strata %in% selected_strata)] <- 0
  p_strata_design <- p_strata_design / sum(p_strata_design)

  final_sample <- c()
  p_i_given_strata <- rep(NA, I)
  psus <- c()

  for (strata in selected_strata) {
    inds <- which(stratum_assignments == strata)

    # Induce informative sampling based on X and PC scores
    incl_score <- -0.6 * score1[inds] + 0.3 * score2[inds] + 0.4 * score3[inds] + 0.6 * X1[inds]
    inclusion_probs <- plogis(incl_score)
    inclusion_probs <- inclusion_probs / sum(inclusion_probs)
    inclusion_probs <- inclusion_probs * I_n * 2
    inclusion_probs[inclusion_probs > 1] <- 1

    sampled_in_strata <- inds[rbinom(length(inds), 1, inclusion_probs) == 1]
    final_sample <- c(final_sample, sampled_in_strata)
    psus <- c(psus, psu_assignments[sampled_in_strata])
    p_i_given_strata[inds] <- inclusion_probs
  }

  p_i <- p_strata_design[stratum_assignments] * p_i_given_strata
  survey_weights <- 1 / p_i

  dat.sim <- data.frame(
    ID = final_sample,
    X = X1[final_sample],
    strata = stratum_assignments[final_sample],
    psu = sub(".*\\_", "", psus),
    weight = survey_weights[final_sample]
  ) %>%
    filter(!is.na(weight)) %>%
    mutate(weight = weight / mean(weight))

  Y_sample <- data.frame(Y_obs[final_sample,])
  colnames(Y_sample) <- paste0("Y", 1:L)

  return(cbind(dat.sim, Y_sample))
}



source(here::here("R", "utils.R"))

nsim = 200
parallel = TRUE
ncores = parallelly::availableCores() - 1
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'Preston', 'weighted', 'none')
fit_types = c('design', 'design', 'design', 'weighted', 'none')
options(survey.lonely.psu = "adjust")



scenario = 1
SNR_sigma = 0.5
family = "gaussian"
B = 500
I_n = 50
sampled_strata = 30

L = 50
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

I = 10e6 # superpopulation size
lst = generate_superpopulation(scenario = scenario,
                               family = family,
                               I = 10e5)
fpcay = refund::fpca.face(lst$Y_obs, center = TRUE)

I_n = 50


sim_res <- list() ## store simulation results
for (iter in 1:nsim) {
  set.seed(iter)
  x = try({
    data <- simulate_survey_fosr_inf(
      X_des = lst$X_des,
      Y_obs = lst$Y_obs,
      fpca_Y = fpcay,
      ## dimension of the functional domain
      I_n = I_n,
      # subjects in each psu-strata combination
      # number PSU
      num_strata = 30,
      sampled_strata = 30,
      stratum_assignments = lst$stratum_assignments,
      psu_assignments = lst$psu_assignments,
      dirichlet_probs = lst$dirichlet_probs
    )


    betaTilde = map(
      .x = fit_types,
      .f = get_betatilde,
      data = data,
      family = family
    )

    betaHat = map(.x = betaTilde, get_betahat)

    if (parallel) {
      plan(multisession, workers = ncores)
      res = future_map2(
        .x = boot_types,
        .y = betaHat,
        .f = run_boots,
        data = data,
        family = family,
        num_boots = B,
        set_seed = TRUE,
        seed = 2025,
        L = 50,
        .options = furrr_options(seed = TRUE)
      )

      cis = future_map2(
        .x = res,
        .y = betaHat,
        .f = get_cis,
        smooth_for_ci = TRUE,
        .options = furrr_options(seed = TRUE)
      )
      plan(sequential)
    } else{
      res = map2(
        .x = boot_types,
        .y = betaHat,
        .f = run_boots,
        data = data,
        family = family,
        num_boots = B,
        set_seed = TRUE,
        seed = 2025,
        L = 50
      )

      cis = map2(
        .x = res,
        .y = betaHat,
        .f = get_cis,
        smooth_for_ci = TRUE
      )
    }



    # get_plot_df = function(sim_res, tname){
    #   num_var = nrow(sim_res$betaHat)
    #   map_dfr(.x = 1:num_var,
    #           .f = function(r){
    #             data.frame(s = 1:length(sim_res$betaHat[r,]),
    #                        beta = sim_res$betaHat[r,],
    #                        lower = sim_res$betaHat[r, ] - 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
    #                        upper = sim_res$betaHat[r, ] + 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
    #                        lower.joint = sim_res$betaHat[r, ] - sim_res$qn[r] *
    #                          sqrt(diag(sim_res$betaHat.var[, , r])),
    #                        upper.joint = sim_res$betaHat[r, ] + sim_res$qn[r] * sqrt(diag(sim_res$betaHat.var[, , r]))) %>%
    #               mutate(name = sim_res$betaHat %>% rownames() %>% .[r])
    #           }) %>%
    #     mutate(varname = tname)
    # }
    # plot_objs = map2(.x = cis, .y = boot_types, .f = get_plot_df) %>% bind_rows()
    #
    # btrue =
    #   beta_true %>%
    #   as_tibble() %>%
    #   mutate(beta = c("beta_0", "beta_1")) %>%
    #   pivot_longer(cols = -beta) %>%
    #   mutate(x = as.numeric(sub(".*V", "", name)),
    #          type = "betaTrue")
    #
    # btrue_df =
    #   btrue %>%
    #   mutate(name = if_else(beta == "beta_0", "(Intercept)", "X")) %>%
    #   select(-beta) %>%
    #   rename(s = x, beta = value)
    # plot_objs %>%
    #   ggplot(aes(x = s, y = beta)) +
    #   geom_line()  +
    #   facet_grid(varname~name) +
    #   geom_ribbon(aes(x = s, ymin = lower, ymax = upper), alpha = .3,
    #               fill = "gray10") +
    #   geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint), alpha = .3,
    #               fill = "gray20") +
    #   labs(x = "L", y = "Estimate") +
    #   theme_classic() +
    #   geom_line(data = btrue_df, aes(x = s, y = beta), color = "red")
    #
    #
    # plot_objs %>%
    #   left_join(btrue_df %>% rename(truth = beta), by = c("name", "s")) %>%
    #   mutate(cover = if_else(between(truth, lower, upper), 1, 0),
    #          cover_joint = if_else(between(truth, lower.joint, upper.joint), 1, 0)) %>%
    #   ggplot(aes(x = s, y = beta)) +
    #   # geom_point(aes(color = factor(cover))) +
    #   geom_line(aes(color = factor(cover), group = 1))  +
    #   facet_grid(varname~name) +
    #   geom_ribbon(aes(x = s, ymin = lower, ymax = upper), alpha = .3,
    #               fill = "gray10") +
    #   geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint), alpha = .3,
    #               fill = "gray20") +
    #   labs(x = "L", y = "Estimate") +
    #   theme_classic() +
    #   scale_color_discrete(name = "Pointwise coverage") +
    #   theme(legend.position = "bottom")
    # # geom_line(data = btrue_df, aes(x = s, y = beta), color = "red")

    stats = map2(
      .x = cis,
      .y = boot_types,
      .f = get_coverage_stats,
      beta_true = beta_true
    ) %>%
      list_rbind() %>%
      mutate(n = nrow(data), n_boot = B)

    sim_res[[iter]] <- stats
    sim_all[[iter]] <- cis
  })
  rm(x)
  print(iter)
}

result =
  sim_res %>%
  keep(., is.data.frame) %>%
  list_rbind(names_to = "id")




write_rds(result, here::here("results", "simulations", "test1.rds"))
write_rds(sim_all, here::here("results", "simulations", "test1_all.rds"))



result %>%
  ggplot(aes(x = boot_type, y = cover_pw)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  facet_grid(.~var) +
  geom_hline(aes(yintercept = 50 * .95))

