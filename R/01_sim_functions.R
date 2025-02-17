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


# get beta hat
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

# betaTilde = get_betatilde(data)
# betaHat = get_betahat(betaTilde)

run_boots = function(data, boot_type, betaHat, family = "gaussian",
                     num_boots = 500, set_seed = TRUE, seed = 2025, L = 50,
                     model_formula = as.formula(paste0('Y~', 'X'))){
  argvals = 1:L
  if(!(boot_type %in% c("BRR", "Rao-Wu-Yue-Beaumont", "Preston", "no_design_bootweight", "weighted", "none"))){
    stop("please specify boot_type as one of 'BRR', 'Rao-Wu-Yue-Beaumont', 'Preston', 'weighted', 'none'")
  }
  out_index <- grep(paste0("^", model_formula[2]), names(data))
  # array with dimensions num_coefficients x length functional domain x num_boots
  betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), num_boots))
  # get BRR replicates if type is BRR - only need to do this once
  if(boot_type == "BRR"){
    sample_data = function(df){
      df %>%
        group_by(strata) %>%
        summarise(
          psu = sample(psu, size = 1)
        )
    }
    set.seed(seed)
    indices = replicate(num_boots, sample_data(data), simplify = FALSE)
  }

  # initiliaze progress bar
  pb <- progress_bar$new(format = "  bootstrapping across L [:bar] :percent eta: :eta",,
                         total = L)
  # get bootstrap estimates at each point along functional domain
  for(l in 1:L){
    pb$tick()
    data$Yl <- unclass(data[,out_index][,argvals[l]])
    svy_design <- svydesign(ids = ~ psu,
                            strata = ~ strata,
                            weights = ~ weight,
                            data = data,
                            nest = TRUE)
    if(set_seed){
      set.seed(seed)
    }
    if(boot_type %in% c("Rao-Wu-Yue-Beaumont", "Preston")){
      bs_design <- suppressWarnings(svy_design %>%
        as_bootstrap_design(
          type = boot_type,
          replicates = num_boots
        ))
      coef.output = svyglm(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
               design = bs_design,
               family = family,
               return.replicates = TRUE,
               control = glm.control(maxit = 5000))

      betaTilde_boot[,l,] <- t(as.matrix(coef.output$replicates[,1:nrow(betaHat)])) # store bootstrap estimates
    } else if(boot_type == "BRR"){
      coefs <- do.call(
        cbind,
        lapply(1:num_boots, function(r) {
          dat_tmp <- data %>%
            inner_join(indices[[r]], by = c("psu", "strata")) %>%
            mutate(weight = weight * 2) # double weights via BRR method

          model <- glm(
            formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
            data = dat_tmp,
            weights = weight,
            family = family,
            control = glm.control(maxit = 5000)
          )

          model$coefficients
        })
      )
      betaTilde_boot[,l,] <- coefs
    } else if(boot_type == "weighted"){
      coefs <- do.call(
        cbind,
        lapply(1:num_boots, function(num) {
          set.seed(l * num)
          index <- sample(row_number(data), dim(data)[1], replace = TRUE, prob = data$weight / sum(data$weight))
          dat.tmp <- data[index, ]

          model <- glm(
            formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
            data = dat.tmp,
            weights = weight,
            family = family)
          model$coefficients
        }))
      betaTilde_boot[,l,] <- coefs
    } else if (boot_type == "none"){
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

# res = c("BRR", "Rao-Wu-Yue-Beaumont", "Preston", "no_design") %>%
#   map(\(x) run_boot_sim(data, betaHat = betaHat, boot_type = x, num_boots = 50, set_seed = TRUE, seed = 2025, L = 50))

# res = run_boots(data, betaHat = betaHat,
#                    boot_type = "BRR", num_boots = 50)

get_cis = function(betaTilde_boot,
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

get_coverage_stats = function(mod_output, beta_true, name){
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

  tibble(
    MISE = MISE %>% unname(),
    mean_joint_se = mean_joint_se,
    mean_pw_se = mean_pw_se,
    cover_joint = cover_joint,
    cover_pw = rowSums(cover_pw),
    var = rownames(beta_true),
    boot_type = name)
}

