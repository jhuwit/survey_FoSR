set.seed(123)

# future direction - could write fn to plot this over varying n, p, L
# Simulation parameters
n <- 1000   # number of observations
p <- 10    # predictors
L <- 100    # outcomes

# Design matrix
X <- matrix(rnorm(n * p), n, p)

# True coefficients
B_true <- matrix(rnorm(p * L, sd = 0.5), p, L)

# Outcomes (Gaussian for now)
Y <- X %*% B_true + matrix(rnorm(n * L), n, L)

# Weights
w <- runif(n, 0.5, 1.5)


# ---------------------------
# Naive: loop over outcomes
# ---------------------------
naive_fit <- function(X, Y, w) {
  coefs <- matrix(NA, ncol = ncol(Y), nrow = ncol(X))
  for (l in 1:ncol(Y)) {
    fit <- glm.fit(
      x = X,
      y = Y[, l],
      weights = w,
      family = gaussian()
    )
    coefs[, l] <- fit$coefficients
  }
  coefs
}

naive_fit_qr <- function(X, Y, w) {
  coefs <- matrix(NA, ncol = ncol(Y), nrow = ncol(X))
  for (l in 1:ncol(Y)) {
    y = Y[, l]
    fit <- lm_ls_multi(X, y, w)
    coefs[, l] <- fit
  }
  coefs
}


# ---------------------------
# Batched: one QR for all Y
# ---------------------------
lm_wls_multi <- function(X, Y, w) {
  W_half <- sqrt(w)
  Xw <- X * W_half
  Yw <- Y * W_half
  qr_decomp <- qr(Xw)
  qr.coef(qr_decomp, Yw)  # (p × L) matrix
}

lm_ls_multi <- function(X, y, w) {
  W_half <- sqrt(w)
  Xw <- X * W_half
  Yw <- y * W_half
  qr_decomp <- qr(Xw)
  qr.coef(qr_decomp, Yw)  # (p × L) matrix
}


# ---------------------------
# Benchmark
# ---------------------------
library(microbenchmark)

res <- microbenchmark(
  naive = naive_fit(X, Y, w),
  batched = lm_wls_multi(X, Y, w),
  naive_qr = naive_fit_qr(X, Y, w),
  times = 10
)
print(res)

# Accuracy check
coef_naive  <- naive_fit(X, Y, w)
coef_batch  <- lm_wls_multi(X, Y, w)

max(abs(coef_naive - coef_batch))  # should be near 1e-12


set.seed(123)

# Simulation parameters
n <- 1000   # number of observations
p <- 10    # predictors
L <- 100    # outcomes
X <- matrix(rnorm(n * p), n, p)

# True coefficients
B_true <- matrix(rnorm(p * L, sd = 0.3), p, L)

# Linear predictor and probabilities
Eta <- X %*% B_true
P <- 1 / (1 + exp(-Eta))

# Outcomes (Bernoulli)
Y <- matrix(rbinom(n * L, size = 1, prob = as.vector(P)), n, L)

# Weights
w <- runif(n, 0.5, 1.5)


# ---------------------------
# Naive: loop over outcomes
# ---------------------------
naive_fit <- function(X, Y, w) {
  coefs <- matrix(NA, nrow = ncol(X) + 1, ncol = ncol(Y))
  for (l in 1:ncol(Y)) {
    fit <- glm.fit(
      x = cbind(Intercept = 1, X),
      y = Y[, l],
      weights = w,
      family = binomial()
    )
    coefs[, l] <- fit$coefficients
  }
  coefs
}


# ---------------------------
# Batched: your glm_batch_multiY
# ---------------------------
glm_batch_multiY <- function(
    X, y_mat, w = NULL,
    family,
    add_intercept = TRUE,
    offset = NULL,
    start = NULL,
    maxit = 50, tol = 1e-8, ridge = 1e-8,
    return_se = FALSE, estimate_phi = FALSE, verbose = FALSE
) {
  stopifnot(is.matrix(X), is.matrix(y_mat))
  n <- nrow(X); p <- ncol(X); B <- ncol(y_mat)
  if (is.null(w)) w <- rep(1, n)

  if (add_intercept) X <- cbind(Intercept = 1, X)
  p <- ncol(X)

  if (is.null(offset)) offset <- matrix(0, n, B)

  linkinv    <- family$linkinv
  mu_eta_fun <- family$mu.eta
  variance   <- family$variance
  eps <- .Machine$double.eps^(2/3)

  Beta <- if (is.null(start)) matrix(0, p, B) else matrix(start, p, B)
  Xt <- t(X)
  converged <- rep(FALSE, B)

  for (it in seq_len(maxit)) {
    Eta <- X %*% Beta + offset
    Eta <- pmin(pmax(Eta, -35), 35)
    Mu  <- linkinv(Eta)

    mu_eta <- mu_eta_fun(Eta); mu_eta[abs(mu_eta) < eps] <- eps
    VarMu  <- variance(Mu);    VarMu[VarMu  < eps] <- eps

    Ww <- (mu_eta^2 / VarMu) * matrix(w, n, B)
    Z  <- Eta + (y_mat - Mu) / mu_eta
    RHS <- Xt %*% (Ww * Z)

    step_max <- rep(0, B)
    active <- which(!converged)
    if (!length(active)) break

    s_all <- (mu_eta^2 / VarMu)
    for (b in active) {
      wb <- w * s_all[, b]
      Xw <- X * wb
      H  <- Xt %*% Xw
      diag(H) <- diag(H) + ridge

      beta_new <- tryCatch({
        Rchol <- chol(H)
        backsolve(Rchol, forwardsolve(t(Rchol), RHS[, b]))
      }, error = function(e) solve(H, RHS[, b], tol = 1e-12))

      step_max[b] <- max(abs(beta_new - Beta[, b]))
      Beta[, b]   <- beta_new
    }

    newly_conv <- (!converged) & (step_max < tol)
    converged[newly_conv] <- TRUE
    if (all(converged)) break
  }

  list(coef = Beta, iter = it, converged = converged)
}


# ---------------------------
# Benchmark
# ---------------------------
library(microbenchmark)

res <- microbenchmark(
  naive = naive_fit(X, Y, w),
  batched = glm_batch_multiY(X, Y, w, family = binomial())$coef,
  times = 10
)
print(res)

# Accuracy check (up to tolerance)
coef_naive  <- naive_fit(X, Y, w)
coef_batch  <- glm_batch_multiY(X, Y, w, family = binomial())$coef
max(abs(coef_naive - coef_batch), na.rm = TRUE)

