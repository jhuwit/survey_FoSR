get_cis_revise0606 = function(betaTilde_boot,
                              betaHat,
                              smooth_for_ci = TRUE,
                              smooth_for_variance = TRUE,
                              L = 50,
                              scaling_factor = 1.3,
                              nknots_min = L/2,
                              nknots_fpca_min = NULL,
                              nknots_min_cov = 35) {
  argvals = 1:L
  B = ncol(betaTilde_boot[1, , ])
  nknots <- nknots_min
  nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
  nknots_fpca <- nknots_fpca_min
  # old: nknots_fpca <- min(round(L / 2), 35)
  betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), ncol(betaTilde_boot[1, , ])))
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  
  # smooth
  for (b in 1:B) {
    betaHat_boot[, , b] <- t(apply(betaTilde_boot[, , b], 1, function(x)
      gam(x ~ s(
        argvals, bs = "tp", k = (nknots + 1)
      ), method = "GCV.Cp")$fitted.values))
  } # smooth each bootstrapped estimate (functional domain)
  
  
  for (r in 1:nrow(betaHat)) {
    if (smooth_for_variance) {
      betaHat.var[, , r] <- scaling_factor * var(t(betaHat_boot[r, , ]))
    } else{
      betaHat.var[, , r] <- scaling_factor * var(t(betaTilde_boot[r, , ]))
    }
  }
  
  N <- 10000
  qn <- rep(0, length = nrow(betaHat))
  set.seed(456)
  for (i in 1:length(qn)) { # the number of beta estimates
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
    
    # Sigma_sd <- Rfast::colVars(x_sample,
    #                            std = TRUE, 
    #                            na.rm = FALSE) # standard deviation: apply(x_sample, 2, sd)
    
    # use bootstaped based std(t)
    Sigma_sd <- Rfast::colVars(est_bs,
                               std = TRUE, 
                               na.rm = FALSE) # standard deviation: apply(x_sample, 2, sd)
    
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


get_cis_revise0505 = function(betaTilde_boot,
                              betaHat,
                              smooth_for_ci = TRUE,
                              smooth_for_variance = TRUE,
                              L = 50,
                              nknots_min = NULL,
                              nknots_fpca_min = NULL,
                              nknots_min_cov = 35) {
  argvals = 1:L
  B = ncol(betaTilde_boot[1, , ])
  nknots <- min(round(L / 2), nknots_min)
  nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
  nknots_fpca <- nknots_fpca_min
  # old: nknots_fpca <- min(round(L / 2), 35)
  betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), ncol(betaTilde_boot[1, , ])))
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  
  # smooth
  for (b in 1:B) {
    betaHat_boot[, , b] <- t(apply(betaTilde_boot[, , b], 1, function(x)
      gam(x ~ s(
        argvals, bs = "tp", k = (nknots + 1)
      ), method = "GCV.Cp")$fitted.values))
  } # smooth each bootstrapped estimate (functional domain)
  
  
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
  for (i in 1:length(qn)) { # the number of beta estimates
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
    
    # Sigma_sd <- Rfast::colVars(x_sample,
    #                            std = TRUE, 
    #                            na.rm = FALSE) # standard deviation: apply(x_sample, 2, sd)
    
    # use bootstaped based std(t)
    Sigma_sd <- Rfast::colVars(est_bs,
                               std = TRUE, 
                               na.rm = FALSE) # standard deviation: apply(x_sample, 2, sd)
    
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



get_cis_revise0422 = function(betaTilde_boot,
                              betaHat,
                              smooth_for_ci = TRUE,
                              smooth_for_variance = TRUE,
                              L = 50,
                              nknots_min = NULL,
                              nknots_fpca_min = NULL,
                              nknots_min_cov = 35) {
  argvals = 1:L
  B = ncol(betaTilde_boot[1, , ])
  nknots <- min(round(L / 2), nknots_min)
  nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
  nknots_fpca <- nknots_fpca_min
  # old: nknots_fpca <- min(round(L / 2), 35)
  betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), ncol(betaTilde_boot[1, , ])))
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  # smooth
  
  for (b in 1:B) {
    betaHat_boot[, , b] <- t(apply(betaTilde_boot[, , b], 1, function(x)
      gam(x ~ s(
        argvals, bs = "tp", k = (nknots + 1)
      ), method = "GCV.Cp")$fitted.values))
  } # smooth each bootstrapped estimate (functional domain)
  
  for (r in 1:nrow(betaHat)) {
    if (smooth_for_variance) {
      betaHat.var[, , r] <- 1.2 * var(t(betaHat_boot[r, , ])) # use smooth beta-estimate
    } else{
      betaHat.var[, , r] <- 1.2 * var(t(betaTilde_boot[r, , ])) # not smooth beta-estimate
    }
  }
  
  N <- 10000
  qn <- rep(0, length = nrow(betaHat))
  set.seed(456)
  for (i in 1:length(qn)) {
    if (smooth_for_ci) { # not consider scaling factor here
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


get_cis_revise0415 = function(betaTilde_boot,
                   betaHat,
                   smooth_for_ci = TRUE,
                   smooth_for_variance = TRUE,
                   L = 50,
                   nknots_min = NULL,
                   nknots_fpca_min = NULL,
                   nknots_min_cov = 35) {
  argvals = 1:L
  B = ncol(betaTilde_boot[1, , ])
  nknots <- min(round(L / 2), nknots_min)
  nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
  nknots_fpca <- min(round(L/2), nknots_fpca_min)
    # old: nknots_fpca <- min(round(L / 2), 35)
  betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), ncol(betaTilde_boot[1, , ])))
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  # smooth
  
  
  for (b in 1:B) {
    betaHat_boot[, , b] <- t(apply(betaTilde_boot[, , b], 1, function(x)
      gam(x ~ s(
        argvals, bs = "tp", k = (nknots + 1)
      ), method = "GCV.Cp")$fitted.values))
  } # smooth each bootstrapped estimate (functional domain)
  
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



visualize_notCover <- function(beta_true = NA,
                               mod_output = NA,
                               subtitle = NA,
                               uncovered_pw = FALSE,
                               L_domain = c(1:10)){
  
  results = data.frame()
  sqrt_diag_var <- apply(mod_output$betaHat.var, 3, function(mat) sqrt(diag(mat)))
  for (p in 1:length(var_name)){
    # Extract beta and SE estimate
    true <- beta_true[p,]
    sqrt_diag_p <- sqrt_diag_var[, p]  # Extract precomputed diagonal elements
    
    # Compute upper and lower bounds in one step
    ## joint-wise
    margin <- mod_output$qn[p] * sqrt_diag_p #joint
    upper <- mod_output$betaHat[p,] + margin
    lower <- mod_output$betaHat[p,] - margin
    
    ## point-wise
    margin_pw <- 1.96 * sqrt_diag_p
    upper_pw <- mod_output$betaHat[p,] + margin_pw
    lower_pw <- mod_output$betaHat[p,] - margin_pw
    
    # return data for plots
    tmp = data.frame(type = var_name[p],
                     beta_estimate = mod_output$betaHat[p,],
                     beta_true = true,
                     beta_pw_lower = lower_pw,
                     beta_pw_upper = upper_pw,
                     beta_jw_lower = lower,
                     beta_jw_upper = upper)
    
    results = rbind(results, tmp)
  }
  
  p1 = results |>
    mutate(x = rep(L_domain, 2)) |>
    ggplot(aes(x = x)) +
    
    ## joint-wise CI
    geom_ribbon(aes(ymin = beta_pw_lower, ymax = beta_pw_upper), fill = "blue", alpha = 0.3) +
    ## point-wise CI
    geom_ribbon(aes(ymin = beta_jw_lower, ymax = beta_jw_upper), fill = "red", alpha = 0.3) +
    ## true values
    geom_line(aes(y = beta_true, color = "True value"), size = 0.5) +
    ## estimated values
    geom_line(aes(y = beta_estimate, color = "Estimated value"), size = 0.5) +
    
    ## Add non-covered jw-points
    geom_point(data = results |> 
                 mutate(x = rep(L_domain, 2)) |> 
                 mutate(cover_jw = ifelse(beta_true >= beta_jw_lower & beta_true <= beta_jw_upper, 1, 0)) |> 
                 filter(cover_jw == 0),
               aes(x = x, y = beta_true), 
               shape = 1, size = 2, color = "green") +
    
    facet_wrap(~type, scales = "free_y") + 
    labs(x = "X", 
         y = "Y", 
         subtitle = subtitle,
         title = "True vs Estimated with Confidence Intervals", color = "Line Type") +
    
    scale_color_manual(values = c("True value" = "black", "Estimated value" = "blue"))
  
  ## add uncovered pw-point
  if (uncovered_pw == TRUE){
    p1 = p1 + geom_point(data = results |> 
                           mutate(x = rep(L_domain, 2)) |> 
                           mutate(cover_pw = ifelse(beta_true >= beta_pw_lower & beta_true <= beta_pw_upper, 1, 0)) |> 
                           filter(cover_pw == 0),
                         aes(x = x, y = beta_true), 
                         shape = 1, size = 2, color = "orange")
  } else {p1 = p1}
  
  return(p1)
}
