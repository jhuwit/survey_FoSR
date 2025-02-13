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
get_betatilde = function(data, family = "gaussian", formula = Y ~ X + (1 | ID),
                         argvals = NULL){
  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)

  L <- ncol(data[,model_formula[2]]) ## number of observations on the functional domain
  if(is.null(argvals)) argvals <- 1:L


  unimm <- function(l){
    data$Yl <- unclass(data[,model_formula[2]][,l])
      options(survey.lonely.psu = "adjust") # process only one observation in the strata
      data$psu = as.character(data$psu)
      data$wid = 1
      data$visit = as.character(data$visit)

      # design <- svydesign(id = ~ psu + ID,
      #                     strata = ~ strata,
      #                     weights = ~ weight + wid,
      #                     data = data,
      #                     nest = TRUE)

      design2 <- svydesign(id = ~ psu + visit,
                          strata = ~ strata,
                          weights = ~ weight + wid,
                          data = data,
                          nest = TRUE)


      fit_uni <- svy2lme(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                         design = design2,
                         sterr = FALSE,
                         method = "general",
                         return.devfun = FALSE)

      # fit_nodes = lmer(formula = as.formula(paste0("Yl ~ ", model_formula[3])),
      #      data = data, control = lmerControl(optimizer = "bobyqa"))

    # betaTilde <- lme4::fixef(fit_uni)
      betaTilde <- coef(fit_uni)
    # return(list(betaTilde = betaTilde, group = as.data.frame(VarCorr(fit_uni))[1,1]))
      return(list(betaTilde = betaTilde))
  }

  num_cores = parallel::detectCores() - 1
  if(parallel == TRUE){
    massmm <- parallel::mclapply(argvals, unimm, mc.cores = num_cores)
  }else{
    massmm <- lapply(argvals, unimm)
  }
  betaTilde <- t(lapply(massmm, '[[', 1) %>% bind_rows())


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

### START HERE

run_boots = function(data, betaHat, family = "gaussian", boot_type, num_boots = 500, set_seed = TRUE, seed = 2025, L = 50,
                     formula = Y ~ X + (1 | ID)){
  data$psu <- as.character(data$psu)
  model_formula <- as.character(formula)
  argvals = 1:L
  set_seed = TRUE
  seed  = 2025
  boot_type = "bootstrap"
  if(!(boot_type %in% c("BRR", "bootstrap", "JKn", "no_design_bootweight", "no_design", "none"))){
    stop("please specify boot_type as one of 'BRR', 'bootstrap', 'JKn', 'no_design_bootweight', 'no_design', or 'none'")
  }
  # array with dimensions num_coefficients x length functional domain x num_boots
  betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), num_boots))
  # get BRR replicates if type is BRR - only need to do this once
  if(boot_type == "BRR"){
    sample_data = function(df){
      df %>%
        select(strata, psu, ID) %>%
        distinct() %>%
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
    data$Yl <- unclass(data[,model_formula[2]][,l])
    data$wid <- 1
    data$visit = as.character(data$visit)
    data_subj =
      data %>%
      select(ID, weight) %>%
      distinct()
    svy_design <- svydesign(ids = ~ psu + visit,
                            strata = ~ strata,
                            weights = ~ weight + wid,
                            data = data,
                            nest = TRUE)
    if(set_seed){
      set.seed(seed)
    }
    if(boot_type %in% c("bootstrap", "JKn")){

      model <- svy2lme(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                         design = svy_design,
                         sterr = FALSE,
                         method = "general",
                         return.devfun = TRUE)

      # env <- environment(model$devfun)
      # env$ctrl$rhoend <- 1e-4

      bs_design <- svy_design %>%
        as.svrepdesign(
          type = boot_type,
          replicates = num_boots
        )
      b <- try({boot2lme(model, bs_design, verbose = TRUE)})

      if(!inherits(b, 'try-error')){
        betaTilde_boot[,l,] <- t(b$beta) # store bootstrap estimates
      } else{
        betaTilde_boot[,l,] <- matrix(NA, nrow = nrow(betaHat), ncol = num_boots)
      }

    }
    if(boot_type == "BRR"){
      coefs <- do.call(
        cbind,
        lapply(1:num_boots, function(r) {
          dat_tmp <- data %>%
            inner_join(indices[[r]], by = c("psu", "strata")) %>%
            mutate(weight = weight * 2) # double weights via BRR method
          dat_tmp$w_scaled <- dat_tmp$weight / mean(dat_tmp$weight)

          model <- lmer(
            formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
            data = dat_tmp,
            weights = w_scaled,
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
          )

          lme4::fixef(model)
        })
      )
      betaTilde_boot[,l,] <- coefs
    } else if(boot_type == "no_design_bootweight"){
      coefs <- do.call(
        cbind,
        lapply(1:num_boots, function(num) {
          set.seed(l * num)
          sampled_subs <- sample(data_subj$ID, nrow(data_subj), replace = TRUE, prob = data_subj$weight / sum(data_subj$weight))
          boot_list <- lapply(seq_along(sampled_subs), function(i) {
            subj <- sampled_subs[i]
            df_subj <- data[data$ID == subj, ]
            df_subj$boot_id <- i
            # Optionally assign a new subject ID for the bootstrapped dataset.
            # This distinguishes between multiple copies of the same subject.
            return(df_subj)
          })

          # Combine the bootstrapped subjects into one data frame
          boot_data <- do.call(rbind, boot_list)
          boot_data$w_scaled = boot_data$weight / mean(boot_data$weight)
          model <- lmer(
            formula = stats::as.formula(paste0("Yl ~ ", "X + (1 | boot_id)")),
            data = boot_data,
            weights = w_scaled,
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
          lme4::fixef(model)
        }))
      betaTilde_boot[,l,] <- coefs
    } else if (boot_type == "no_design"){
      coefs <- do.call(
        cbind,
        lapply(1:num_boots, function(num) {
          set.seed(l * num)
          sampled_subs <- sample(data_subj$ID, nrow(data_subj), replace = TRUE)
          boot_list <- lapply(seq_along(sampled_subs), function(i) {
            subj <- sampled_subs[i]
            df_subj <- data[data$ID == subj, ]
            df_subj$boot_id <- i
            # Optionally assign a new subject ID for the bootstrapped dataset.
            # This distinguishes between multiple copies of the same subject.
            return(df_subj)
          })

          # Combine the bootstrapped subjects into one data frame
          boot_data <- do.call(rbind, boot_list)
          boot_data$w_scaled = boot_data$weight / mean(boot_data$weight)
          model <- lmer(
            formula = stats::as.formula(paste0("Yl ~ ", "X + (1 | boot_id)")),
            data = boot_data,
            weights = w_scaled,
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
          lme4::fixef(model)
        }))
      betaTilde_boot[,l,] <- coefs
    } else if (boot_type == "none"){
      coefs <- do.call(
        cbind,
        lapply(1:num_boots, function(num) {
          set.seed(l * num)
          sampled_subs <- sample(data_subj$ID, nrow(data_subj), replace = TRUE, prob = data_subj$weight / sum(data_subj$weight))
          boot_list <- lapply(seq_along(sampled_subs), function(i) {
            subj <- sampled_subs[i]
            df_subj <- data[data$ID == subj, ]
            df_subj$boot_id <- i
            # Optionally assign a new subject ID for the bootstrapped dataset.
            # This distinguishes between multiple copies of the same subject.
            return(df_subj)
          })

          # Combine the bootstrapped subjects into one data frame
          boot_data <- do.call(rbind, boot_list)
          model <- lmer(
            formula = stats::as.formula(paste0("Yl ~ ", "X + (1 | boot_id)")),
            data = boot_data,
            control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
          lme4::fixef(model)
        }))
      betaTilde_boot[,l,] <- coefs
    }
  }
  return(betaTilde_boot)
}

# run_boots(data, betaHat = betaHat, boot_type = "BRR", num_boots = 10,
#              set_seed = TRUE, seed = 2025)
# res = c("BRR", "Rao-Wu-Yue-Beaumont", "Preston", "no_design") %>%
#   map(\(x) run_boot_sim(data, betaHat = betaHat, boot_type = x, num_boots = 50, set_seed = TRUE, seed = 2025, L = 50))

get_cis = function(betaTilde_boot,
                   smooth_for_ci = FALSE,
                   smooth_for_variance = TRUE,
                   L = 50,
                   betaHat,
                   nknots_min = NULL,
                   nknots_min_cov = 35) {
  argvals = 1:L
  B = ncol(betaTilde_boot[1, , ])
  nknots <- min(round(L / 2), nknots_min)
  nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov)
  nknots_fpca <- min(round(L / 2), 35)
  # betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), ncol(betaTilde_boot[1, , ])))
  betaHat.var <- array(NA, dim = c(L, L, nrow(betaHat)))
  # approx NAs in betaTilde with linear interpolation
  smooth_nas = function(x){
    if (any(is.na(x))) {
      x <- zoo::na.approx(x, x = argvals, rule = 2)
    }
    x
  }
  # impute NAs in betaTilde_boot
  betaTilde_boot <- aperm(apply(betaTilde_boot, c(1, 3), smooth_nas), c(2, 1, 3))

  # smooth betaTilde_boot
  temp <- aperm(betaTilde_boot, c(1, 3, 2))  # dims: 2 x B x L

  # Apply the smoothing function over subject and replicate.
  # For each (subject, replicate) pair, x is a vector of length T.
  res <- apply(temp, c(1, 2), function(x) {
    gam(x ~ s(argvals, bs = "tp", k = (nknots + 1)), method = "GCV.Cp")$fitted.values
  }) # dim L x 2 x B

  # repermute dimensions
  betaHat_boot <- aperm(res, c(2, 1, 3)) # dim 2 x B x L
  rm(temp); rm(res)

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
  cover_joint = rep(NA, nrow(beta_true))
  cover_pw = matrix(FALSE, nrow(beta_true), L)

  for(p in 1:length(cover_joint)){
    true = beta_true[p,]
    upper = mod_output$betaHat[p,] + mod_output$qn[p] * sqrt(diag(mod_output$betaHat.var[,,p]))
    lower = mod_output$betaHat[p,] - mod_output$qn[p] * sqrt(diag(mod_output$betaHat.var[,,p]))
    upper_pw = mod_output$betaHat[p,] + 1.96 * sqrt(diag(mod_output$betaHat.var[,,p]))
    lower_pw = mod_output$betaHat[p,] - 1.96 * sqrt(diag(mod_output$betaHat.var[,,p]))
    cover_upper = which(upper > true)
    cover_lower = which(lower < true)
    cover_joint[p] <- length(intersect(cover_lower, cover_upper)) == L

    cover_upper_pw <- which(upper_pw > true)
    cover_lower_pw <- which(lower_pw < true)
    cover_pw[p,intersect(cover_lower_pw, cover_upper_pw)] <- TRUE
  }
  tibble(
    MISE = MISE %>% unname(),
    cover_joint = cover_joint,
    cover_pw = rowSums(cover_pw),
    var = rownames(beta_true),
    boot_type = name)
}

