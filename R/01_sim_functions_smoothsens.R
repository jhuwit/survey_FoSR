library(dplyr)

source(here::here("R", "create_survey_settings.R"))
settings_small =
  settings |>
  dplyr::filter(
    (family == "gaussian" & strata_scale == .125 & strata_sigma == 0.05 &
       snr_b == 0.5 & snr_eps == 1 & inf_level == 10 & In == 100) |
      (family %in% c("binomial", "poisson") & strata_scale == .125 &
         strata_sigma == 0.05 & snr_b == 0.5 & inf_level == 10 & In == 100)
  )


sensitivity_grid = tidyr::expand_grid(
  nknots_spec = c("conservative", "current", "flexible"),
  splines = c("tp", "cr", "ps"),
  method = c("GCV.Cp", "REML")
)


full_grid = tidyr::expand_grid(
  nknots_spec = c("conservative", "current", "flexible"),
  splines = c("tp", "cr", "ps"),
  method = c("GCV.Cp", "REML"),
  family = c("gaussian", "binomial", "poisson"),
  strata_scale = 0.125,
  strata_sigma = 0.05,
  inf_level = 10,
  In = 100,
  snr_b = 0.5,
  snr_eps = 1,
  len = c(50, 100, 1440))

full_grid =
  full_grid |>
  mutate(snr_eps = if_else(family == "gaussian", snr_eps, NA_real_))


smaller_grid =
  full_grid |>
  filter((splines == "tp" & method == "GCV.Cp") |
           (nknots_spec == "current" & splines == "cr" & method == "GCV.Cp") |
           (nknots_spec == "current" & splines == "tp" & method == "REML") |
           (nknots_spec == "current" & splines == "ps" & method == "REML")
  ) |>
  mutate(fold = row_number())
# key_scenarios <- list(
#   baseline = list(nknots = "current", splines = "tp", method = "GCV.Cp"),
#   fewer_knots = list(nknots = "conservative", splines = "tp", method = "GCV.Cp"),
#   more_knots = list(nknots = "flexible", splines = "tp", method = "GCV.Cp"),
#   alt_basis = list(nknots = "current", splines = "cr", method = "GCV.Cp"),
#   alt_method = list(nknots = "current", splines = "tp", method = "REML"),
#   penalized = list(nknots = "current", splines = "ps", method = "REML")
# )



get_betahat_sensitivity = function(betaTilde, L,
                                   nknots_spec = "current",
                                   splines = "tp",
                                   method = "GCV.Cp") {
  argvals = 1:L

  # Knot specifications
  nknots_options <- list(
    "conservative" = min(round(L/4), 20),
    "current" = min(round(L/2), 35),
    "flexible" = min(round(3*L/4), 50),
    "very_flexible" = min(L - 5, 70)
  )

  nknots <- nknots_options[[nknots_spec]]

  betaHat <- t(apply(betaTilde, 1, function(x) {
    fit <- mgcv::gam(x ~ s(argvals, bs = splines, k = (nknots + 1)),
                     method = method)
    fit$fitted.values
  }))

  rownames(betaHat) <- rownames(betaTilde)
  colnames(betaHat) <- 1:L
  return(betaHat)
}

get_cis_sensitivity = function(betaTilde_boot,
                   betaHat,
                   smooth_for_ci = TRUE,
                   smooth_for_variance = TRUE,
                   L,
                   nknots_spec = "current",
                   splines = "tp",
                   method = "GCV.Cp",
                   # nknots_min = 35,
                   nknots_min_cov = 35,
                   nknots_fpca = 35,
                   mult_fac = 1.2) {
  argvals = 1:L
  B = ncol(betaTilde_boot[1, , ])
  # nknots <- min(round(L / 2), nknots_min)
  nknots_options <- list(
    "conservative" = min(round(L/4), 20),
    "current" = min(round(L/2), 35),
    "flexible" = min(round(3*L/4), 50),
    "very_flexible" = min(L - 5, 70)
  )

  nknots <- nknots_options[[nknots_spec]]

  nknots_fpca <- min(round(L / 2), nknots_fpca)
  betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), ncol(betaTilde_boot[1, , ])))
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  # smooth


  # for (b in 1:B) {
  #   betaHat_boot[, , b] <- t(apply(betaTilde_boot[, , b], 1, function(x)
  #     gam(x ~ s(
  #       argvals, bs = "tp", k = (nknots + 1)
  #     ), method = "GCV.Cp")$fitted.values))
  # }

  for (b in 1:B) {
    betaHat_boot[, , b] <- t(apply(betaTilde_boot[, , b], 1, function(x)
      bam(x ~ s(
        argvals, bs = splines, k = (nknots + 1)
      ), method = method)$fitted.values))
  }

  for (r in 1:nrow(betaHat)) {
    if (smooth_for_variance) {
      betaHat.var[, , r] <- mult_fac * var(t(betaHat_boot[r, , ]))
    } else{
      betaHat.var[, , r] <- mult_fac * var(t(betaTilde_boot[r, , ]))
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
    # x_sample: N x p matrix
    # x_mean: length-p vector
    # Sigma_sd: length-p vector

    # Vectorized computation
    z <- sweep(x_sample, 2, x_mean, "-")     # center
    z <- sweep(z, 2, Sigma_sd, "/")          # standardize
    un <- apply(abs(z), 1, max)              # row-wise max of absolute values

    qn[i] <- stats::quantile(un, 0.95)
  }

  return(list(
    betaHat = betaHat,
    betaHat.var = betaHat.var,
    qn = qn
  ))

}
