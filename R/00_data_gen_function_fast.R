# Generate 'structure' only (no Y_obs)
get_p_i = function(i, probs) probs[i] * (1 + sum((probs[-i]) / (1-probs[-i])))

generate_population_structure <- function(I = 1e6, L = 50, family = "gaussian",
                                          seed = 4574, num_strata = 30,
                                          strata_sigma = 0.05, psu_factor = 0.5,
                                          strata_scale = 0.125, snr_b = 1, snr_eps = 1) {

  stopifnot(family %in% c("gaussian","poisson","binomial"))

  set.seed(seed)

  # scalar covariate (we keep this vector; it's small vs I x L)
  X_des <- cbind(1, rnorm(I, 0, 2))

  # beta fixed
  grid <- seq(0, 1, length = L)
  beta_fixed <- matrix(NA, 2, L)
  beta_fixed[1, ] <- -0.15 - 0.1 * sin(2 * pi * grid) - 0.1 * cos(2 * pi * grid)
  beta_fixed[2, ] <- dnorm(grid, 0.6, 0.15) / 20

  rownames(beta_fixed) <- c("Intercept", "x")

  # strata sizes (draw exactly as original): dirichlet then sample assignments
  set.seed(seed)
  dirichlet_probs <- gtools::rdirichlet(1, rep(4, num_strata))
  set.seed(seed)
  stratum_assignments <- sample(1:num_strata, I, replace = TRUE, prob = dirichlet_probs)

  # build PSUs per stratum exactly like original and collect their labels
  psu_assignments <- rep(NA, I)
  num_psu_by_stratum <- integer(num_strata)
  psu_label_counter <- 0L
  psu_label_map <- list() # to map (stratum, local_psu) -> global psu index

  for (s in seq_len(num_strata)) {
    set.seed(seed + s)
    num_in_strata <- sum(stratum_assignments == s)
    set.seed(seed + s)
    num_psu <- round(runif(1, 75, 125))
    num_psu_by_stratum[s] <- num_psu

    set.seed(seed + s)
    dps <- gtools::rdirichlet(1, rep(10, num_psu))

    set.seed(seed + s)
    psu_in_stratum <- sample(1:num_psu, num_in_strata, replace = TRUE, prob = dps)

    # assign labels like "s_local"
    psu_assignments[stratum_assignments == s] <- paste0(s, "_", psu_in_stratum)

    # record mapping of local->global indices
    local_ids <- unique(psu_in_stratum)
    # we'll map later in the same order as unique(psu_assignments) if needed
    psu_label_map[[s]] <- seq_len(num_psu)  # simply 1:num_psu
  }

  # Build basis and Phi for functional random effects
  nbasis <- 5
  basis <- fda::create.bspline.basis(c(0,1), nbasis)
  Phi <- fda::eval.basis(grid, basis)

  # strata-specific scores and random effects
  set.seed(seed)
  strata_scores <- matrix(rnorm(num_strata * nbasis, 0, strata_sigma), num_strata, nbasis)
  strata_random_effects <- strata_scores %*% t(Phi)   # num_strata x L

  # Now compute total_psu (same as original)
  total_psu <- length(unique(psu_assignments))

  # psu-specific scores and random effects (same seed usage)
  psu_sigma <- sqrt(strata_sigma^2 * psu_factor)
  set.seed(seed)
  psu_scores <- matrix(rnorm(total_psu * nbasis, 0, psu_sigma), total_psu, nbasis)
  psu_random_effects <- psu_scores %*% t(Phi)   # total_psu x L

  # stratum-specific slope scaling
  set.seed(seed)
  stratum_scaling <- rnorm(num_strata, mean = 1, sd = strata_scale)

  # produce beta1_by_stratum (num_strata x L) and beta1_by_indiv (I x L) if you want
  beta1_by_stratum <- matrix(rep(stratum_scaling, each = L), nrow = num_strata, byrow = TRUE) *
    matrix(rep(beta_fixed[2, ], times = num_strata), nrow = num_strata, byrow = TRUE)

  # return the structure needed for streaming
  list(
    I = I, L = L, grid = grid, beta_fixed = beta_fixed,
    X_des = X_des, dirichlet_probs = dirichlet_probs,
    stratum_assignments = stratum_assignments,
    psu_assignments = psu_assignments,
    num_psu_by_stratum = num_psu_by_stratum,
    strata_random_effects = strata_random_effects,
    psu_random_effects = psu_random_effects,
    stratum_scaling = stratum_scaling,
    beta1_by_stratum = beta1_by_stratum,
    Phi = Phi
  )
}


# Stream-sampling by stratum using the stored structure. This creates the final sampled data frame
sample_by_strata_stream <- function(structure,
                                    I_n = 500,
                                    num_selected_psu = 2,
                                    seed = 1,
                                    inf_level = 1,
                                    compression = 2,
                                    family = "gaussian") {
  # unpack
  I <- structure$I; L <- structure$L
  X_des <- structure$X_des
  stratum_assignments <- structure$stratum_assignments
  psu_assignments <- structure$psu_assignments
  dirichlet_probs <- structure$dirichlet_probs
  beta_fixed <- structure$beta_fixed
  strata_random_effects <- structure$strata_random_effects
  psu_random_effects <- structure$psu_random_effects
  stratum_scaling <- structure$stratum_scaling

  final_sample <- integer(0)
  psus <- character(0)

  p1 <- rep(NA_real_, I)
  p2 <- rep(NA_real_, I)
  p_overall <- rep(NA_real_, I)

  num_strata <- nrow(strata_random_effects)
  psu_ids_all <- unique(psu_assignments)

  # helper to compute stage-1 psu selection probability under PPS WOR approx:
  psu_prob_selected_all <- vector("numeric", length = length(psu_ids_all))
  names(psu_prob_selected_all) <- psu_ids_all
  # compute psu sizes by stratum once
  for (s in seq_len(num_strata)) {
    inds_in_stratum <- which(stratum_assignments == s)
    psu_sizes <- table(psu_assignments[inds_in_stratum])
    psu_ids <- names(psu_sizes)
    psu_probs <- as.numeric(psu_sizes) / sum(psu_sizes)
    # inclusion prob for sampling without replacement of num_selected_psu PSUs:

    psu_prob_selected <- 1 - (1 - psu_probs) ^ num_selected_psu
    names(psu_prob_selected) <- psu_ids
    psu_prob_selected_all[psu_ids] <- psu_prob_selected
  }

  # iterate strata, generate Y rows for members of that strata on the fly
  for (strata in 1:num_strata) {
    inds_in_stratum <- which(stratum_assignments == strata)

    # get psu ids inside this stratum
    psu_sizes <- table(psu_assignments[inds_in_stratum])
    psu_ids <- names(psu_sizes)

    # select PSUs (keep same sampling behavior as your original: sample(..., prob = psu_sizes))
    set.seed(strata + seed)
    selected_psus <- sample(psu_ids, size = num_selected_psu, replace = FALSE, prob = psu_sizes)

    # for each selected PSU, we need to evaluate inclusion_probs and actually sample units
    for (psu in selected_psus) {
      inds_in_psu <- which(psu_assignments == psu & stratum_assignments == strata)

      # compute per-individual fixed effect (matrix of 1 x L rows by individuals)
      # fixef_signal_i = intercept + X * beta2(t)
      # we don't build whole matrix; compute rowMeans(Y) requires the entire row of Y,
      # but we can compute lin_pred (vector length L) per individual and then simulate Y row.
      # form fixef_signal per individual: replicate beta_fixed[1,] + X * beta_fixed[2,]
      X_vals <- X_des[inds_in_psu, 2]

      # get strata-level random effect for this stratum (length L)
      strata_re <- strata_random_effects[strata, , drop = TRUE]

      # map psu to global index for psu_random_effects: factor(as.numeric(factor(psu_assignments)))
      # create mapping of global psu index: we replicate original ordering by using factor levels from unique(psu_assignments)
      psu_factor_levels <- levels(factor(psu_assignments))
      global_psu_index <- which(psu_factor_levels == psu)
      psu_re <- psu_random_effects[global_psu_index, , drop = TRUE]

      # per-individual slope modification:
      slope_re_indiv <- (stratum_scaling[strata] - 1) * beta_fixed[2, ]

      # Combined random effect vector for individuals in this psu:
      # ranef_row = strata_re + psu_re + slope_re_indiv
      ranef_row_common <- strata_re + psu_re + slope_re_indiv  # length L

      # Now build lin_pred for each individual (do sequentially to limit memory)
      n_psu <- length(inds_in_psu)
      # compute inclusion base scores (we will need rowMeans over Y, so first generate Y rows)
      # We'll generate Y_i (length L) per individual, compute rowMean, then compute inclusion probs
      # But we must use same sigma generation as original: sigma = sd(as.vector(fixef_signal + ranef)) / snr_eps
      # To compute sigma, we need sd of all lin_preds in the stratum/psu or globally. Original computed sd_lp based on whole population.
      # To maintain same sigma as original, we approximate by computing sigma per stratum using fixef_signal + ranef_row_common.
      # To match exactly original you'd compute sd_lp across entire population; here we approximate per-stratum.
      # If you want exact match to previous code's sigma, compute sd_lp once globally before streaming (requires computing fixef_signal + ranef for all individuals,
      # which defeats streaming). So you must accept sigma computed per-stratum or pre-compute global sd_lp outside streaming.
      # ----> We'll compute global sd_lp once before streaming using a light-weight approach:
      #
      # (We assume the structure object included a precomputed global sd_lp; if not, you can compute it earlier with a cheap pass.)
      stopifnot(!is.null(structure$global_sd_lp))  # require that generate_population_structure computed this
      sigma <- structure$global_sd_lp / (ifelse(is.null(structure$snr_eps), 1, structure$snr_eps))

      # Now generate Y rows sequentially and compute inclusion probs
      # For speed, precompute the matrix for beta_fixed[2, ] to multiply by X_vals
      beta2_mat <- matrix(rep(beta_fixed[2, ], n_psu), nrow = n_psu, byrow = TRUE)
      intercept_mat <- matrix(rep(beta_fixed[1, ], n_psu), nrow = n_psu, byrow = TRUE)
      Xcol_mat <- matrix(X_vals, nrow = n_psu, ncol = L, byrow = FALSE)
      # compute fixef_signal for each individual: intercept + X * beta2(t)
      fixef_signal_mat <- intercept_mat + Xcol_mat * beta2_mat

      # add ranef common to all individuals (broadcast)
      lin_preds_mat <- fixef_signal_mat + matrix(rep(ranef_row_common, each = n_psu), nrow = n_psu, byrow = TRUE)

      # Simulate Y rows now (Gaussian / binomial / poisson)
      if (family == "gaussian") {
        set.seed(strata + seed + which(selected_psus == psu)) # keep reproducible
        Y_psu <- matrix(rnorm(n = n_psu * L, mean = as.vector(t(lin_preds_mat)), sd = sigma),
                        nrow = n_psu, ncol = L, byrow = TRUE)
        # compute row means
        y_mean <- rowMeans(Y_psu)
      } else if (family == "binomial") {
        p_true <- plogis(as.vector(t(lin_preds_mat)))
        set.seed(strata + seed + which(selected_psus == psu))
        Y_psu <- matrix(rbinom(n = n_psu * L, size = 1, prob = p_true),
                        nrow = n_psu, ncol = L, byrow = TRUE)
        y_mean <- rowMeans(Y_psu)
      } else if (family == "poisson") {
        lam_true <- exp(as.vector(t(lin_preds_mat)))
        set.seed(strata + seed + which(selected_psus == psu))
        Y_psu <- matrix(rpois(n = n_psu * L, lambda = lam_true),
                        nrow = n_psu, ncol = L, byrow = TRUE)
        y_mean <- rowMeans(Y_psu)
      }

      # Now compute inclusion scores & probs (clamp for binomial)
      if (inf_level == 0) {
        inclusion_probs <- rep(1 / n_psu, n_psu)
      } else {
        if (family == "gaussian") {
          incl_score <- y_mean * inf_level
        } else if (family == "poisson") {
          incl_score <- log(y_mean) * inf_level
        } else if (family == "binomial") {
          m <- pmin(pmax(y_mean, 1e-6), 1 - 1e-6)
          incl_score <- qlogis(m) * inf_level
        }
        score_compressed <- pmax(pmin(incl_score, compression), -compression)
        inclusion_probs <- plogis(score_compressed)
      }

      inclusion_probs_adj <- inclusion_probs / sum(inclusion_probs) * I_n
      inclusion_probs_adj[inclusion_probs_adj > 1] <- 1

      # sample via Bernoulli as before
      set.seed(strata + seed + which(selected_psus == psu))
      draw <- rbinom(length(inclusion_probs_adj), 1, inclusion_probs_adj)
      sampled_units_local_idx <- which(draw == 1)
      sampled_units_global_idx <- inds_in_psu[sampled_units_local_idx]

      # append sampled ids and psu labels
      final_sample <- c(final_sample, sampled_units_global_idx)
      psus <- c(psus, rep(psu, length(sampled_units_local_idx)))

      # stage probabilities: p1 is psu-level pi computed earlier; p2 is the actual Bernoulli prob used
      p_psu <- psu_prob_selected_all[psu]
      p1[inds_in_psu] <- p_psu
      p2[inds_in_psu] <- inclusion_probs_adj
      p_overall[inds_in_psu] <- p_psu * inclusion_probs_adj

      # discard Y_psu, lin_preds_mat, etc — loop continues
    }
  }

  # Build final sampled data frame
  survey_weights <- 1 / p_overall

  dat.sim <- data.frame(
    ID = final_sample,
    X = X_des[final_sample,2],
    strata = stratum_assignments[final_sample],
    psu = sub(".*\\_", "", psus),
    weight = survey_weights[final_sample],
    p_stage1 = p1[final_sample],
    p_stage2 = p2[final_sample]
  )

  Y_sample <- data.frame(matrix(NA_real_, nrow = length(final_sample), ncol = L))
  colnames(Y_sample) <- paste0("Y", 1:L)
  # fill Y_sample by regenerating the same Y rows for the sampled indices (or store above when sampled)
  # In the loop we generated Y_psu; to avoid storing them earlier, regenerate Y rows for final_sample here
  # This requires re-computing lin_pred for each sampled unit — cheap since final_sample << I usually.

  for (k in seq_along(final_sample)) {
    iidx <- final_sample[k]
    # compute lin_pred for unit iidx (same formula as used above)
    s <- stratum_assignments[iidx]
    # intercept + x * beta2 + strata_random + psu_random + slope_re
    xval <- X_des[iidx, 2]
    # find global_psu_index for psu_assignments[iidx]
    psu_label <- psu_assignments[iidx]
    psu_factor_levels <- levels(factor(psu_assignments))
    gp_idx <- which(psu_factor_levels == psu_label)
    fixef_vec <- beta_fixed[1, ] + xval * beta_fixed[2, ]
    ranef_vec <- strata_random_effects[s, ] + psu_random_effects[gp_idx, ] + (stratum_scaling[s] - 1) * beta_fixed[2, ]
    linp <- fixef_vec + ranef_vec

    if (family == "gaussian") {
      set.seed(seed + iidx) # seed scheme to regenerate deterministic row (you can choose same as earlier)
      Yrow <- rnorm(L, mean = linp, sd = structure$global_sd_lp / (ifelse(is.null(structure$snr_eps), 1, structure$snr_eps)))
    } else if (family == "binomial") {
      ptrue <- plogis(linp)
      set.seed(seed + iidx)
      Yrow <- rbinom(L, 1, ptrue)
    } else {
      lam <- exp(linp)
      set.seed(seed + iidx)
      Yrow <- rpois(L, lam)
    }
    Y_sample[k, ] <- Yrow
  }

  data_final <- cbind(dat.sim, Y_sample)
  data_final
}

sample_by_strata_stream_fast <- function(structure,
                                         I_n = 500,
                                         num_selected_psu = 2,
                                         seed = 1,
                                         inf_level = 1,
                                         compression = 2,
                                         family = "gaussian") {
  # unpack
  I <- structure$I; L <- structure$L
  X_des <- structure$X_des
  stratum_assignments <- structure$stratum_assignments
  psu_assignments <- structure$psu_assignments
  dirichlet_probs <- structure$dirichlet_probs
  beta_fixed <- structure$beta_fixed
  strata_random_effects <- structure$strata_random_effects
  psu_random_effects <- structure$psu_random_effects
  stratum_scaling <- structure$stratum_scaling

  # checks
  stopifnot(!is.null(structure$global_sd_lp))
  sigma <- structure$global_sd_lp / (ifelse(is.null(structure$snr_eps), 1, structure$snr_eps))

  # Precompute maps / levels once (very cheap)
  psu_levels <- unique(psu_assignments)          # order as in your structure
  psu_map_all <- match(psu_assignments, psu_levels)  # length-I mapping to psu row index
  # Precompute constants
  mean_beta1 <- mean(beta_fixed[1, ])
  mean_beta2 <- mean(beta_fixed[2, ])

  # outputs
  final_sample <- integer(0)
  psus <- character(0)
  p1 <- rep(NA_real_, I)
  p2 <- rep(NA_real_, I)
  p_overall <- rep(NA_real_, I)

  num_strata <- nrow(strata_random_effects)

  # compute psu-level first-stage selection probabilities per stratum
  psu_prob_selected_all <- numeric(length(psu_levels))
  names(psu_prob_selected_all) <- psu_levels
  for (s in seq_len(num_strata)) {
    inds_in_stratum <- which(stratum_assignments == s)
    psu_sizes <- table(psu_assignments[inds_in_stratum])
    psu_ids <- names(psu_sizes)
    psu_probs <- as.numeric(psu_sizes) / sum(psu_sizes)
    # If num_selected_psu == 2 and you have get_p_i, use it; else use exact routine if you want
    if (num_selected_psu == 2 && exists("get_p_i")) {
      # compute pi_i for WOR with n=2 using your get_p_i
      pi_vec <- sapply(seq_along(psu_probs), function(ii) get_p_i(ii, psu_probs))
      names(pi_vec) <- psu_ids
    } else {
      # fallback approx (or use sampling::inclusionprobabilities for exact)
      pi_vec <- pmin(1, num_selected_psu * psu_probs)
      names(pi_vec) <- psu_ids
    }
    psu_prob_selected_all[psu_ids] <- pi_vec
  }

  # iterate strata
  for (strata in 1:num_strata) {
    inds_in_stratum <- which(stratum_assignments == strata)
    psu_sizes <- table(psu_assignments[inds_in_stratum])
    psu_ids <- names(psu_sizes)

    set.seed(strata + seed)
    selected_psus <- sample(psu_ids, size = num_selected_psu, replace = FALSE, prob = psu_sizes)

    # process each selected PSU
    for (psu in selected_psus) {
      inds_in_psu <- which(psu_assignments == psu & stratum_assignments == strata)
      n_psu <- length(inds_in_psu)
      X_vals <- X_des[inds_in_psu, 2]

      # prepare ranef for this psu (length L)
      # map psu to global index in psu_random_effects
      global_psu_idx <- match(psu, psu_levels)
      strata_re <- strata_random_effects[strata, ]
      psu_re <- psu_random_effects[global_psu_idx, ]
      slope_re_common <- (stratum_scaling[strata] - 1) * beta_fixed[2, ]
      ranef_row_common <- strata_re + psu_re + slope_re_common   # length L

      # For each individual in PSU compute the mean of linpred across grid cheaply:
      # mean_lin_i = mean(beta1) + x_i * mean(beta2) + mean(ranef_row_common)
      mean_ranef_common <- mean(ranef_row_common)
      mean_lin_pred_vec <- mean_beta1 + X_vals * mean_beta2 + mean_ranef_common  # length n_psu

      if (inf_level == 0) {
        # uniform sampling inside PSU -> equal inclusion probs
        inclusion_probs <- rep(1 / n_psu, n_psu)
      } else {
        # We need y_mean for each unit. Use distributional shortcuts:
        if (family == "gaussian") {
          # row mean of L normals with means lin_pred_j and sd sigma has distribution:
          # mean = mean_lin_pred_i, sd = sigma / sqrt(L)
          sd_mean <- sigma / sqrt(L)
          # simulate the observed row means (same RNG call structure as before)
          set.seed(strata + seed + which(selected_psus == psu))
          y_mean <- stats::rnorm(n_psu, mean = mean_lin_pred_vec, sd = sd_mean)
        } else if (family == "binomial") {
          # approximate mean of L Bernoulli(p_j) by Normal(mean(p_j), var = sum p_j(1-p_j) / L^2)
          # need p_j = plogis(lin_pred_j) per individual (we must compute lin_pred_j vector)
          y_mean <- numeric(n_psu)
          set.seed(strata + seed + which(selected_psus == psu))
          for (ii in seq_len(n_psu)) {
            xval <- X_vals[ii]
            lin_pred_vec <- beta_fixed[1, ] + xval * beta_fixed[2, ] + ranef_row_common  # length L
            p_vec <- plogis(lin_pred_vec)
            mean_p <- mean(p_vec)
            var_mean <- sum(p_vec * (1 - p_vec)) / (L^2)
            if (var_mean <= 0) var_mean <- 1e-12
            y_mean[ii] <- rnorm(1, mean = mean_p, sd = sqrt(var_mean))
          }
        } else if (family == "poisson") {
          # approximate mean of L Poissons by Normal(mean=lambda_mean, var = sum lambda_j / L^2)
          y_mean <- numeric(n_psu)
          set.seed(strata + seed + which(selected_psus == psu))
          for (ii in seq_len(n_psu)) {
            xval <- X_vals[ii]
            lin_pred_vec <- beta_fixed[1, ] + xval * beta_fixed[2, ] + ranef_row_common
            lam_vec <- exp(lin_pred_vec)
            mean_lam <- mean(lam_vec)
            var_mean <- sum(lam_vec) / (L^2)
            if (var_mean <= 0) var_mean <- 1e-12
            y_mean[ii] <- rnorm(1, mean = mean_lam, sd = sqrt(var_mean))
          }
        } else stop("unknown family")
        # compute inclusion score using y_mean
        if (family == "gaussian") {
          incl_score <- y_mean * inf_level
        } else if (family == "poisson") {
          incl_score <- log(pmax(y_mean, 1e-8)) * inf_level
        } else { # binomial
          # ensure in (0,1)
          m <- pmin(pmax(y_mean, 1e-8), 1 - 1e-8)
          incl_score <- qlogis(m) * inf_level
        }

        score_compressed <- pmax(pmin(incl_score, compression), -compression)
        inclusion_probs <- plogis(score_compressed)
      }

      # adjust to target expected sample size I_n (your original semantics)
      inclusion_probs_adj <- inclusion_probs / sum(inclusion_probs) * I_n
      inclusion_probs_adj[inclusion_probs_adj > 1] <- 1

      # Bernoulli draws to decide selection (same RNG call style as original)
      set.seed(strata + seed + which(selected_psus == psu))
      draw <- rbinom(length(inclusion_probs_adj), 1, inclusion_probs_adj)
      sampled_local_idx <- which(draw == 1)
      sampled_global_idx <- inds_in_psu[sampled_local_idx]

      # store stage probs and overall
      p_psu <- psu_prob_selected_all[psu]
      p1[inds_in_psu] <- p_psu
      p2[inds_in_psu] <- inclusion_probs_adj
      p_overall[inds_in_psu] <- p_psu * inclusion_probs_adj

      # append sample ids and psu labels
      final_sample <- c(final_sample, sampled_global_idx)
      psus <- c(psus, rep(psu, length(sampled_global_idx)))

      # We do NOT create the full Y_psu here; we will regenerate full Y rows only for final_sample later.
    } # end loop psu
  } # end loop strata

  # Build dat.sim
  survey_weights <- 1 / p_overall
  dat.sim <- data.frame(
    ID = final_sample,
    X = X_des[final_sample, 2],
    strata = stratum_assignments[final_sample],
    psu = sub(".*\\_", "", psus),
    weight = survey_weights[final_sample],
    p_stage1 = p1[final_sample],
    p_stage2 = p2[final_sample]
  )

  # Now regenerate full Y rows only for selected units (cheap)
  n_selected <- length(final_sample)
  Y_sample <- matrix(NA_real_, nrow = n_selected, ncol = L)
  psu_factor_levels <- psu_levels  # reuse
  for (k in seq_len(n_selected)) {
    iidx <- final_sample[k]
    s <- stratum_assignments[iidx]
    xval <- X_des[iidx, 2]
    psu_label <- psu_assignments[iidx]
    gp_idx <- match(psu_label, psu_factor_levels)
    fixef_vec <- beta_fixed[1, ] + xval * beta_fixed[2, ]
    ranef_vec <- strata_random_effects[s, ] + psu_random_effects[gp_idx, ] + (stratum_scaling[s] - 1) * beta_fixed[2, ]
    linp <- fixef_vec + ranef_vec

    set.seed(seed + iidx) # deterministic per sampled unit
    if (family == "gaussian") {
      Y_sample[k, ] <- rnorm(L, mean = linp, sd = sigma)
    } else if (family == "binomial") {
      ptrue <- plogis(linp)
      Y_sample[k, ] <- rbinom(L, 1, ptrue)
    } else {
      lam <- exp(linp)
      Y_sample[k, ] <- rpois(L, lam)
    }
  }

  Y_sample_df <- as.data.frame(Y_sample)
  colnames(Y_sample_df) <- paste0("Y", seq_len(L))
  data_final <- cbind(dat.sim, Y_sample_df)
  data_final
}


# You **must** compute a global sd_lp for consistency with your old code.
# If you want exact equivalence to previous sd_lp which used sd(as.vector(fixef_signal)+ranef) across the whole population,
# you need to compute that once; but doing so requires a single pass through all individuals to compute their lin_pred (cheap: no Y).
# Here's a quick pass to compute global sd_lp (no Y retained):

compute_global_sd_lp <- function(structure) {
  I <- structure$I; L <- structure$L
  beta_fixed <- structure$beta_fixed
  X_des <- structure$X_des
  stratum_assignments <- structure$stratum_assignments
  psu_assignments <- structure$psu_assignments
  strata_re <- structure$strata_random_effects
  psu_re <- structure$psu_random_effects
  stratum_scaling <- structure$stratum_scaling

  # compute a vector of lin_pred values in the same order as original by streaming small chunks:
  # We'll compute sd of all linpred values without storing them: accumulate sum and sumsq
  sum1 <- 0; sumsq <- 0; nvals <- 0
  chunk_size <- 10000
  idx_seq <- seq_len(I)
  for (start in seq(1, I, by = chunk_size)) {
    end <- min(start + chunk_size - 1, I)
    ids <- start:end
    # compute lin_pred for these ids
    svec <- stratum_assignments[ids]
    pvec <- psu_assignments[ids]
    xvec <- X_des[ids, 2]
    # mapping psu -> global index
    psu_levels <- levels(factor(psu_assignments))
    gp_idx <- match(pvec, psu_levels)
    # vectorized: compute fixef + ranef per id
    fixef_mat <- matrix(rep(beta_fixed[1, ], length(ids)), nrow = length(ids), byrow = TRUE) +
      (xvec %o% beta_fixed[2, ])
    ranef_mat <- strata_re[svec, , drop = FALSE] + psu_re[gp_idx, , drop = FALSE] +
      (matrix(rep((stratum_scaling[svec] - 1), each = L), nrow = length(ids)) * matrix(rep(beta_fixed[2, ], length(ids)), nrow = length(ids), byrow = TRUE))
    lin_mat <- fixef_mat + ranef_mat
    vec <- as.vector(t(lin_mat))
    sum1 <- sum1 + sum(vec)
    sumsq <- sumsq + sum(vec^2)
    nvals <- nvals + length(vec)
  }
  mu <- sum1 / nvals
  sd_lp <- sqrt(max((sumsq - nvals * mu^2) / (nvals - 1), 0))
  sd_lp
}
compute_global_sd_lp_fast <- function(structure, chunk_size = 50000) {
  I <- structure$I; L <- structure$L
  beta_fixed <- structure$beta_fixed
  X_des <- structure$X_des[, 2]  # take just covariate column
  stratum_assignments <- structure$stratum_assignments
  psu_assignments <- structure$psu_assignments
  strata_re <- structure$strata_random_effects
  psu_re <- structure$psu_random_effects
  stratum_scaling <- structure$stratum_scaling

  # Precompute mapping psu -> row index in psu_re
  psu_levels <- unique(psu_assignments)   # much faster than factor() every loop
  psu_map <- match(psu_assignments, psu_levels)

  sum1 <- 0; sumsq <- 0; nvals <- 0

  for (start in seq(1, I, by = chunk_size)) {
    end <- min(start + chunk_size - 1, I)
    ids <- start:end

    svec <- stratum_assignments[ids]
    gp_idx <- psu_map[ids]
    xvec <- X_des[ids]

    # Compute linear predictor pieces efficiently
    # Intercept + slope*X (broadcasted manually)
    fixef_mat <- matrix(beta_fixed[1, ] + outer(xvec, beta_fixed[2, ]), nrow = length(ids))

    ranef_mat <- strata_re[svec, , drop = FALSE] +
      psu_re[gp_idx, , drop = FALSE] +
      (stratum_scaling[svec] - 1) %o% beta_fixed[2, ]

    lin_mat <- fixef_mat + ranef_mat

    vec <- as.vector(lin_mat)
    sum1 <- sum1 + sum(vec)
    sumsq <- sumsq + sum(vec * vec)
    nvals <- nvals + length(vec)
  }

  mu <- sum1 / nvals
  sqrt(max((sumsq - nvals * mu^2) / (nvals - 1), 0))
}

# structure$global_sd_lp <- compute_global_sd_lp(structure)
# structure$global_sd_lp2 <- compute_global_sd_lp_fast(structure)
#
# structure$snr_eps <- 1  # match input
#
# structure <- generate_population_structure(I = 10e6, L = 1440, scenario = 1, family = "gaussian", seed = 4574, num_strata = 30,
#                                            strata_sigma = 0.05, psu_factor = 0.5, strata_scale = 0.125, snr_b = 1, snr_eps = 1)
# structure$global_sd_lp <- compute_global_sd_lp(structure)
# structure$snr_eps <- 1  # match input
#
# sampled_data <- sample_by_strata_stream(structure,
#                                         I_n = 500,
#                                         num_selected_psu = 2,
#                                         seed = 1,
#                                         inf_level = 1,
#                                         compression = 3,
#                                         family = "gaussian")
#
