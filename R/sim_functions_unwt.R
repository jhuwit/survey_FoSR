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


generate_population = function(n = 500, # size of superpopulation
                               L = 50, # length of functional domain
                               scenario = 1, # shape of function
                               family = "gaussian",
                               seed = 4574,
                               SNR_sigma = 1
){
  stopifnot("scenario must be either 1 or 2" = scenario %in% c(1, 2))
  stopifnot("family must be either 'gaussian', 'poisson', or 'binomial'" = family %in% c("gaussian", "poisson", "binomial"))

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
  set.seed(seed)
  X_des = cbind(1, rnorm(n, 0, 2)) # design matrix: Intercept + predicted variables
  true_fixef <- X_des %*% beta_true

  set.seed(seed + 1)
  if (family == "gaussian") {
    Y_obs <- matrix(
      rnorm(
        n = nrow(true_fixef) * L,
        mean = as.vector(t(true_fixef)),
        sd = SNR_sigma
      ),
      nrow = nrow(true_fixef),
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "binomial") {
    p_true = plogis(true_fixef)
    Y_obs <- matrix(
      rbinom(
        n = nrow(p_true) * L,
        size = 1,
        prob = as.vector(t(p_true))
      ),
      nrow = nrow(p_true),
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "poisson") {
    lam_true <- exp(true_fixef)
    Y_obs <- matrix(
      rpois(n = nrow(lam_true) * L, lambda = as.vector(t(lam_true))),
      nrow = nrow(lam_true),
      ncol = L,
      byrow = TRUE
    )
  }

  dat.sim <- data.frame(
    ID = 1:n,
    X = X_des[,2]
  )

  # Generate fixed effects and true Y.
  Y_df <- data.frame(Y_obs)
  colnames(Y_df) <- paste0("Y", 1:L)

  sample_data <- cbind(dat.sim, Y_df)

  return(list(sample_data = sample_data, beta_true = beta_true))
}

# get beta hat
get_betatilde = function(data, family = "gaussian",
                         model_formula = as.formula(paste0('Y~', 'X'))){


  out_index <- grep(paste0("^", model_formula[2]), names(data)) # indices that start with the outcome name
  if(length(out_index) != 1){ # functional observations stored in multiple columns
    L <- length(out_index)
  }else{ # observations stored as a matrix in one column using the I() function
    L <- ncol(data[,out_index])
  }
  argvals <- 1:L

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

  massmm <- lapply(argvals, unimm)


  # Obtain betaTilde, fixed effects estimates
  betaTilde <- t(do.call(rbind, lapply(massmm, '[[', 1)))
  colnames(betaTilde) <- argvals
  return(betaTilde)
}

get_betahat = function(betaTilde, L = 50,
                       nknots_min = 35){
  argvals = 1:L
  nknots <- min(round(L/4), nknots_min)

  betaHat <- t(apply(betaTilde, 1, function(x) mgcv::gam(x ~ s(argvals, bs = "tp",
                                                               k = (nknots + 1)),
                                                         method = "GCV.Cp")$fitted.values))

  rownames(betaHat) <- rownames(betaTilde)
  colnames(betaHat) <- 1:L
  return(betaHat)
}


get_betahat_analytic = function(betaTilde, L = 50, splines = "tp",
                                smooth_method = "GCV.Cp",
                                nknots_min = 35,
                                nknots_min_cov = NULL
                                ){
    p <- nrow(betaTilde) # Number of fixed effects parameters
    betaHat <- matrix(NA, nrow = p, ncol = L)
    lambda <- rep(NA, p)
    argvals = 1:L
    nknots <- min(nknots_min,  round(L/4))
    # Number of knots for covariance matrix
    nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
    nknots_fpca <- min(round(L / 2), 35)
    for (r in 1:p) {
      fit_smooth <- gam(
        betaTilde[r,] ~ s(argvals, bs = splines, k = (nknots + 1)),
        method = smooth_method
      )
      betaHat[r,] <- fit_smooth$fitted.values
      lambda[r] <- fit_smooth$sp # Smoothing parameter
    }

    sm <- smoothCon(
      s(argvals, bs = splines, k = (nknots + 1)),
      data = data.frame(argvals = argvals),
      absorb.cons = TRUE
    )
    S <- sm[[1]]$S[[1]] # Penalty matrix
    B <- sm[[1]]$X # Basis functions
    rm(fit_smooth, sm)
  return(betaHat)
}
# b_hat = get_betahat_analytic(betaTilde)
# betaTilde = get_betatilde(data)
# betaHat = get_betahat(betaTilde)

run_boots = function(data, betaHat, family = "gaussian",
                     num_boots = 500, set_seed = TRUE, seed = 2025, L = 50,
                     model_formula = as.formula(paste0('Y~', 'X')),
                     parallel = FALSE, ncores_inner = 4){
  argvals = 1:L

  out_index <- grep(paste0("^", model_formula[2]), names(data))
  # array with dimensions num_coefficients x length functional domain x num_boots
  betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), num_boots))
  # get BRR replicates if type is BRR - only need to do this once

  # initiliaze progress bar
  pb <- progress_bar$new(format = "  bootstrapping across L [:bar] :percent eta: :eta",,
                         total = L)
  # get bootstrap estimates at each point along functional domain
  if(parallel){
    plan(multisession, workers = ncores_inner)
    for(l in 1:L) {
      pb$tick()
      data$Yl <- unclass(data[, out_index][, argvals[l]])

      if (set_seed) {
        set.seed(seed)  # Change seed per L
      }

      coefs <- do.call(
        cbind,
        future_lapply(1:num_boots, function(num) {
          set.seed(l * num)
          index <- sample(nrow(data), nrow(data), replace = TRUE)
          dat.tmp <- data[index, ]

          model <- glm(
            formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
            data = dat.tmp,
            family = family
          )
          model$coefficients
        }, future.seed = TRUE,
        future.scheduling = 1,
        future.globals = c("data", "model_formula", "family"))
      )

      betaTilde_boot[, l, ] <- coefs
    }
    plan(sequential)
  } else {
    for(l in 1:L){
    pb$tick()
    data$Yl <- unclass(data[,out_index][,argvals[l]])
    if(set_seed){
      set.seed(seed)
    }
      coefs <- do.call(
        cbind,
        lapply(1:num_boots, function(num) {
          set.seed(l * num)
          index <- sample(row_number(data), dim(data)[1], replace = TRUE)
          dat.tmp <- data[index, ]

          model <- glm(
            formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
            data = dat.tmp,
            family = family)
          model$coefficients
        }))
      betaTilde_boot[,l,] <- coefs
    }
  }
  return(betaTilde_boot)
}

eigenval_trim <- function(V) {
  edcomp <- eigen(V, symmetric = TRUE)
  eigen.positive <- which(edcomp$values > 0)
  q <- ncol(V)

  if (length(eigen.positive) == q) {
    return(V)
  } else if (length(eigen.positive) == 0) {
    return(tcrossprod(edcomp$vectors[, 1]) * edcomp$values[1])
  } else if (length(eigen.positive) == 1) {
    return(tcrossprod(as.vector(edcomp$vectors[, 1])) * as.numeric(edcomp$values[1]))
  } else {
    return(matrix(
      edcomp$vectors[, eigen.positive] %*%
        tcrossprod(diag(edcomp$values[eigen.positive]), edcomp$vectors[, eigen.positive]),
      ncol = q
    ))
  }
}


library(zoo)  # For rolling mean (moving average)
library(Matrix)  # For handling covariance matrices

# Function to apply moving average (Gaussian-weighted)
sandwich_smooth <- function(mat, kernel_width = 5) {
  n <- nrow(mat)
  smoothed_mat <- mat

  # Apply moving average row-wise
  for (i in 1:n) {
    smoothed_mat[i, ] <- rollapply(mat[i, ], width = kernel_width, FUN = mean, fill = NA, align = "center", partial = TRUE)
  }

  # Apply moving average column-wise
  for (j in 1:n) {
    smoothed_mat[, j] <- rollapply(smoothed_mat[, j], width = kernel_width, FUN = mean, fill = NA, align = "center", partial = TRUE)
  }

  return(smoothed_mat)
}

get_cis_analytic = function(betaHat, betaTilde, L = 50, data){
  p = nrow(betaHat) # number of beta coefs
  designmat = cbind(rep(1, nrow(data)), data$X) # design matrix X
  Y_mat = data %>% # Y matrix
    select(starts_with("Y")) %>%
    as.matrix()

  Y_pred = designmat %*% betaTilde # y = x beta
  resid_mat = Y_mat - Y_pred # residual matrix
  n = nrow(Y_mat)

  R_hat = matrix(NA, nrow = L, ncol = L) # covariance matrix of Ys
  for (l1 in 1:L) {
    for (l2 in 1:L) {
      if (l1 == l2) {
        R_hat[l1, l2] = (1 / (n - p)) * sum((resid_mat[,l1]^2)) # if l1 = l2, use residual variance as estimate
      } else{
        R_hat[l1, l2] = mean(resid_mat[,l1] * resid_mat[,l2]) # method of moments estimator
      }
    }
  }


  XTX_inv <- solve(t(designmat) %*% designmat) # precompute
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  for (l1 in 1:len) {
    for (l2 in 1:len) {
      betaHat.var[l1, l2, ] <- diag((XTX_inv %*% t(designmat) * R_hat[l1, l2]) %*% designmat %*% XTX_inv)
    }
  }
  # smooth
  betaHat.var.smooth <- array(NA, dim = dim(betaHat.var))
  for (r in 1:nrow(betaHat)) {
    betaHat.var.smooth[, , r] <- sandwich_smooth(betaHat.var[, , r])
    # Ensure positive definite
    betaHat.var[, , r] <- pmax(0, betaHat.var[, , r])
    betaHat.var.smooth[, , r] <- eigenval_trim(betaHat.var.smooth[, , r])
  }

  return(list(
    betaHat = betaHat,
    betaHat.var = betaHat.var.smooth
  ))
}



get_cis_boot = function(betaTilde_boot,
                   betaHat,
                   smooth_for_ci = TRUE,
                   smooth_for_variance = TRUE,
                   L = 50,
                   nknots_min = NULL,
                   nknots_min_cov = 35) {
  argvals = 1:L
  B = ncol(betaTilde_boot[1, , ])
  nknots <- min(round(L / 2), nknots_min)
  nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
  nknots_fpca <- min(round(L / 2), 35)
  betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), ncol(betaTilde_boot[1, , ])))
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  # smooth


  for (b in 1:B) {
    betaHat_boot[, , b] <- t(apply(betaTilde_boot[, , b], 1, function(x)
      gam(x ~ s(
        argvals, bs = "tp", k = (nknots + 1)
      ), method = "GCV.Cp")$fitted.values))
  }
  for (r in 1:nrow(betaHat)) {
    if (smooth_for_variance) {
      betaHat.var[, , r] <- 1.2 * var(t(betaHat_boot[r, , ]))
    } else{
      betaHat.var[, , r] <- 1.2 * var(t(betaTilde_boot[r, , ]))
    }
  }

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

  return(list(
    betaHat = betaHat,
    betaHat.var = betaHat.var,
    qn = qn
  ))

}
#
# mod_output = get_cis(res[[1]], betaHat = betaHat)
#

get_coverage_stats_fui = function(mod_output, beta_true, L, type = "boot"){
  if(type == "boot"){
    MISE <- rowMeans((mod_output$betaHat - beta_true)^2)
    mean_pw_se = apply(mod_output$betaHat.var, 3, function(mat) mean(sqrt(diag(mat))))
    mean_joint_se <- numeric(nrow(beta_true))  # Store mean CI width

    cover_joint <- logical(nrow(beta_true))  # Using logical instead of NA-filled vector
    cover_pw <- matrix(FALSE, nrow(beta_true), L)

    # Extract the diagonal elements of each variance matrix upfront
    sqrt_diag_var <- apply(mod_output$betaHat.var, 3, function(mat) sqrt(diag(mat)))

    # Loop efficiently
    for (p in seq_along(cover_joint)) {
      true <- beta_true[p,]
      sqrt_diag_p <- sqrt_diag_var[, p]  # Extract precomputed diagonal elements

      # Compute upper and lower bounds in one step
      margin <- mod_output$qn[p] * sqrt_diag_p
      mean_joint_se[p] <- mean(margin)

      upper <- mod_output$betaHat[p,] + margin
      lower <- mod_output$betaHat[p,] - margin

      margin_pw <- 1.96 * sqrt_diag_p
      upper_pw <- mod_output$betaHat[p,] + margin_pw
      lower_pw <- mod_output$betaHat[p,] - margin_pw

      # Compute joint and pointwise coverage efficiently
      covered <- (lower < true) & (upper > true)
      covered_pw <- (lower_pw < true) & (upper_pw > true)

      cover_joint[p] <- all(covered)  # Joint coverage check
      cover_pw[p, ] <- covered_pw     # Assign logical vector directly
    }
    cover_pw_int = which(cover_pw[1, ])
    cover_pw_x = which(cover_pw[2, ])
    result = tibble(
      MISE = MISE %>% unname(),
      mean_joint_se = mean_joint_se,
      mean_pw_se = mean_pw_se,
      cover_joint = cover_joint,
      cover_pw = rowSums(cover_pw) / L,
      var = rownames(beta_true)) %>%
      mutate(vector_col = map(var, ~ ifelse(.x  == "Intercept", list(cover_pw_int), list(cover_pw_x))))
  } else if (type == "analytic"){
    MISE <- rowMeans((mod_output$betaHat - beta_true)^2)
    mean_pw_se = apply(mod_output$betaHat.var, 3, function(mat) mean(sqrt(diag(mat))))

    cover_joint <- rep(NA, 2) # Using logical instead of NA-filled vector
    cover_pw <- matrix(FALSE, nrow(beta_true), L)

    # Extract the diagonal elements of each variance matrix upfront
    sqrt_diag_var <- apply(mod_output$betaHat.var, 3, function(mat) sqrt(diag(mat)))

    # Loop efficiently
    for (p in seq_along(cover_joint)) {
      true <- beta_true[p,]
      sqrt_diag_p <- sqrt_diag_var[, p]  # Extract precomputed diagonal elements


      margin_pw <- 1.96 * sqrt_diag_p
      upper_pw <- mod_output$betaHat[p,] + margin_pw
      lower_pw <- mod_output$betaHat[p,] - margin_pw

      covered_pw <- (lower_pw < true) & (upper_pw > true)

      cover_pw[p, ] <- covered_pw     # Assign logical vector directly
    }
    cover_pw_int = which(cover_pw[1, ])
    cover_pw_x = which(cover_pw[2, ])
    result = tibble(
      MISE = MISE %>% unname(),
      mean_joint_se = NA_real_,
      mean_pw_se = mean_pw_se,
      cover_joint = NA_real_,
      cover_pw = rowSums(cover_pw) / L,
      var = rownames(beta_true)) %>%
      mutate(vector_col = map(var, ~ ifelse(.x  == "Intercept", list(cover_pw_int), list(cover_pw_x))))
  }
  return(result)

}


get_coverage_stats_bayes = function(res, beta_true, L = 50){
  beta_hat_bayes = res$beta.hat
  lb = res$beta.LB
  ub = res$beta.UB
  width = ub - lb
  se = width / 2 / 1.96
  MISE = rowMeans((beta_hat_bayes - beta_true)^2)
  mean_pw_se = apply(se, 1, mean)

  cover_pw = matrix(FALSE, nrow(beta_true), L)
  for(p in seq_along(MISE)){
    cover_upper = which(ub[p,] > beta_true[p,])
    cover_lower = which(lb[p,] < beta_true[p,])
    cover_pw[p, intersect(cover_lower, cover_upper)] = TRUE
  }

  cover_pw_int = which(cover_pw[1, ])
  cover_pw_x = which(cover_pw[2, ])

  tibble(
    MISE = MISE %>% unname(),
    mean_pw_se = mean_pw_se,
    cover_pw = rowSums(cover_pw) / L,
    var = rownames(beta_true)) %>%
    mutate(vector_col = map(var, ~ ifelse(.x  == "Intercept", list(cover_pw_int), list(cover_pw_x))))
}

get_coverage_stats_famm = function(res, beta_true, L = 50){
  coef_pffr = coef(res, n1 = L)
  betaHat_pffr <- betaHat_pffr.se <- matrix(NA, nrow = 2, ncol = L)
  betaHat_pffr[1,] = as.vector(coef_pffr$smterms$`Intercept(yindex)`$value) + unname(res$coef[1])
  betaHat_pffr[2,] = as.vector(coef_pffr$smterms$`X(yindex)`$value)
  betaHat_pffr.se[1,] = as.vector(coef_pffr$smterms$`Intercept(yindex)`$se)
  betaHat_pffr.se[2,] = as.vector(coef_pffr$smterms$`X(yindex)`$se)
  MISE = rowMeans((betaHat_pffr - beta_true)^2)
  mean_pw_se = apply(betaHat_pffr.se, 1, mean)

  cover_pw = matrix(FALSE, nrow(beta_true), L)
  for(p in seq_along(MISE)){
    cover_upper = which((betaHat_pffr[p,]+(1.96*betaHat_pffr.se[p,])) > beta_true[p,])
    cover_lower = which((betaHat_pffr[p,]-(1.96*betaHat_pffr.se[p,])) < beta_true[p,])
    cover_pw[p, intersect(cover_lower, cover_upper)] = TRUE
  }

  cover_pw_int = which(cover_pw[1, ])
  cover_pw_x = which(cover_pw[2, ])


  tibble(
    MISE = MISE %>% unname(),
    mean_pw_se = mean_pw_se,
    cover_pw = rowSums(cover_pw) / L,
    var = rownames(beta_true)) %>%
    mutate(vector_col = map(var, ~ ifelse(.x  == "Intercept", list(cover_pw_int), list(cover_pw_x))))
}

pffrcoef_to_tf = function(coef) {

  coef =
    coef %>%
    mutate(
      id = "ID"
    ) %>%
    tf_nest(value:se, .id = id, .arg = yindex.vec) %>%
    select(coef = value, se)

  coef

}
