fui.weight <- function(formula,
                       data,
                       random = TRUE,
                       weighted = TRUE,
                       family = "gaussian",
                       var = TRUE,
                       analytic = FALSE,
                       parallel = FALSE,
                       silent = FALSE,
                       argvals = NULL,
                       nknots_min = NULL,
                       nknots_min_cov = 35,
                       boot_type = 'FPB',
                       smooth_method = "GCV.Cp",
                       splines = "tp",
                       design_mat = FALSE,
                       residuals = FALSE,
                       G_return = FALSE,
                       num_boots = 500,
                       seed = 1,
                       subj_ID = NULL,
                       num_cores = 1,
                       caic = FALSE,
                       REs = FALSE,
                       non_neg = 0,
                       MoM = 2){

  # If doing parallel computing, set up the number of cores
  if(parallel & !is.integer(num_cores) ) num_cores <- parallel::detectCores() - 1

  # For non-Gaussian family, only do bootstrap inference
  if(family != "gaussian") analytic <- FALSE

  # Organize the input from the model formula
  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)

  # stop function if there are column names with "." to avoid issues with covariance G() and H() calculations below
  dep_str <- deparse(model_formula[3])
  if(grepl(".", dep_str, fixed = TRUE)){
    # make sure it isn't just a call to all covariates with "Y ~. "
    dep_str_rm <- substr(dep_str, 3, nchar(dep_str)) # remove first character of parsed formula string and check
    if(grepl(".", dep_str_rm, fixed = TRUE)){
      stop('Remove the character "." from all non-functional covariate names and rerun fui() function
           -For example, change the name "X.1" to "X_1"
           -The string "." *should* be kept in the functional outcome names (e.g., "Y.1" *is* proper naming).')
    }
  }
  rm(dep_str)

  ##############################################################################
  ## Step 1
  ##############################################################################
  if(silent == FALSE) print("Step 1: Fit Massively Univariate Mixed Models")

  # Obtain the dimension of the functional domain
  out_index <- grep(paste0("^", model_formula[2]), names(data)) # indices that start with the outcome name
  if(length(out_index) != 1){ # functional observations stored in multiple columns
    L <- length(out_index)
  }else{ # observations stored as a matrix in one column using the I() function
    L <- ncol(data[,out_index])
  }
  if(analytic & !is.null(argvals) & var)  message("'argvals' argument is currently only supported for bootstrap. `argvals' ignored: fitting model to ALL points on functional domain")
  if(is.null(argvals) | analytic){
    # Set up the functional domain when not specified
    argvals <- 1:L
  }else{
    # when only using a subset of the functional domain
    if(max(argvals) > L) stop("Maximum index specified in argvals is greater than the total number of columns for the functional outcome")
    L <- length(argvals)
  }
  if(family == "gaussian" & analytic & L > 400 & var)   message("Your functional data is dense! Consider subsampling along the functional domain (i.e., reduce columns of outcome matrix) or using bootstrap inference.")

  # Create a matrix to store AICs
  AIC_mat <- matrix(NA, nrow = L, ncol = 2)

  # # when consider weighted functions, build a survey design object
  # if (weighted == TRUE){
  #   nhats_design <- svydesign(ids = ~ psu,
  #                             strata = ~ strata,
  #                             weights = ~ weight,
  #                             data = data,
  #                             nest=TRUE)
  # }

  ### Create a function that fit a mixed model at location l
  ### Input: l:location of the functional domain
  ### Output: A list containing point estimates, variance estimates, etc.
  unimm <- function(l){
    data$Yl <- unclass(data[,out_index][,l])

    if (random == FALSE){
      if (weighted == TRUE){ # Determine if weights need to be taken into account
        options(survey.lonely.psu = "adjust") # process only one observation in the strata
        nhats_design <- svydesign(ids = ~ psu,
                                  strata = ~ strata,
                                  weights = ~ weight,
                                  data = data,
                                  nest=TRUE)

        if(family == "gaussian"){
          fit_uni <- suppressMessages(svyglm(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                             design = nhats_design,
                                             control = glm.control(maxit = 5000)))
        }else{
          fit_uni <- suppressMessages(svyglm(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                             design = nhats_design,
                                             family = family,
                                             control = glm.control(maxit = 5000)))
        }
      }else{
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
      if (weighted == TRUE) {
        betaTilde <- coef(fit_uni)
      }else{
        betaTilde <- coef(fit_uni)
      }

      re_df <- aic_met <- resids <- NA
      if(residuals) resids <- as.numeric(residuals(fit_uni)) # these are residuals from including the random effects (i.e., with BLUPs): not JUST from fixed effects -- can verify by comparing with nlme::lme() and seeing 2 columns of residuals in lme: mod$residuals
      if(caic) aic_met <- as.numeric(cAIC4::cAIC(fit_uni)$caic)[1]
      if(REs) re_df <- ranef(fit_uni) ## random effects
    }
    if (random == TRUE){
      if (weighted == TRUE){ # Determine if weights need to be taken into account
        options(survey.lonely.psu = "adjust") # process only one observation in the strata
        data$psu = as.character(data$psu)

        nhats_design <- svydesign(id = ~ psu + subjectid,
                                  strata = ~ strata,
                                  weights = ~ weight + wid,
                                  data = data,
                                  nest = TRUE)

        fit_uni <- suppressWarnings( svy2lme(formula = stats::as.formula(paste0("Y.Y.1 ~ ", model_formula[3])),
                                             design = nhats_design,
                                             return.devfun = TRUE) )

      }else{
        if(family == "gaussian"){
          fit_uni <- suppressMessages(lmer(formula = as.formula(paste0("Yl ~ ", model_formula[3])),
                                           data = data, control = lmerControl(optimizer = "bobyqa")))
        }else{
          fit_uni <- suppressMessages(glmer(formula = as.formula(paste0("Yl ~ ", model_formula[3])),
                                            data = data, family = family, control = glmerControl(optimizer = "bobyqa")))
        }
        betaTilde <- lme4::fixef(fit_uni)
      }

      ## fixed effects estimates
      if (weighted == TRUE) {
        betaTilde <- coef(fit_uni)
      }else{
        betaTilde <- lme4::fixef(fit_uni)
      }

      re_df <- aic_met <- resids <- NA
      if(residuals) resids <- as.numeric(residuals(fit_uni)) # these are residuals from including the random effects (i.e., with BLUPs): not JUST from fixed effects -- can verify by comparing with nlme::lme() and seeing 2 columns of residuals in lme: mod$residuals
      if(caic) aic_met <- as.numeric(cAIC4::cAIC(fit_uni)$caic)[1]
      if(REs) re_df <- ranef(fit_uni) ## random effects
    }

    if (analytic == TRUE){
      return(print("This session haven not been developed"))
    }
    if (analytic == FALSE && random == TRUE && weighted == FALSE) {
      return(list(betaTilde = betaTilde,
                  group = as.data.frame(VarCorr(fit_uni))[1,1],
                  aic = stats::AIC(fit_uni),
                  bic = stats::BIC(fit_uni),
                  residuals = resids,
                  caic = aic_met,
                  re_df = re_df))
    } else {
      return(list(betaTilde = betaTilde,
                  #group = as.data.frame(VarCorr(fit_uni))[1,1],
                  #aic = stats::AIC(fit_uni),
                  #bic = stats::BIC(fit_uni),
                  residuals = resids,
                  caic = aic_met,
                  re_df = re_df))
    }
  }

  # Fit massively univariate mixed models
  if(parallel == TRUE){
    massmm <- mclapply(argvals, unimm, mc.cores = num_cores)
  }else{
    massmm <- lapply(argvals, unimm)
  }

  # Obtain betaTilde, fixed effects estimates
  betaTilde <- t(do.call(rbind, lapply(massmm, '[[', 1)))
  colnames(betaTilde) <- argvals


  ##############################################################################
  ## Step 2
  ##############################################################################
  if(silent == FALSE) print("Step 2: Smoothing")

  # Penalized splines smoothing and extract components (analytic)
  nknots <- min(round(L/2), nknots_min) ## number of knots for regression coefficients
  nknots_cov <- ifelse(is.null(nknots_min_cov), 35, nknots_min_cov) ## number of knots for covariance matrix
  nknots_fpca <- min(round(L/2), 35)

  ## smoothing parameter, spline basis, penalty matrix (analytic)
  if(analytic == TRUE){
    return(print("This session haven not been developed"))
  }else{
    betaHat <- t(apply(betaTilde, 1, function(x) mgcv::gam(x ~ s(argvals, bs = splines, k = (nknots + 1)), method = smooth_method)$fitted.values))
  }
  rownames(betaHat) <- rownames(betaTilde)
  colnames(betaHat) <- 1:L


  ##############################################################################
  ## Step 3
  ##############################################################################
  if(var == TRUE){ ## skip the step when var = FALSE
    if(analytic == FALSE){

      ##########################################################################
      ## Bootstrap Inference
      ##########################################################################
      if(silent == FALSE) print("Step 3: Inference (Bootstrap)")

      # Check to see if group contains ":" which indicates hierarchical structure and group needs to be specified
      if(weighted == FALSE & random == TRUE){
        group <- massmm[[1]]$group
        if(grepl(":", group, fixed = TRUE)){
          if(is.null(subj_ID)){
            group <- str_remove(group, ".*:") # assumes the ID name is to the right of the ":"
          }else if(!is.null(subj_ID)){
            group <- subj_ID # use user specified if it exists
          }else{
            message("You must specify the argument: ID")
          }
        }
        ID.number <- unique(data[,group])
      } else {
        ID.number <- unique(data$subjectid)
      }

      idx_perm <- t( replicate(num_boots, sample.int(length(ID.number), length(ID.number), replace = TRUE)) )
      # generate Resampling index matrix
      B <- num_boots
      betaHat_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))

      if(family != "gaussian" & boot_type %in% c("wild", "reb")){
        stop('Non-gaussian outcomes only supported for some bootstrap procedures. \n Set argument `boot_type` to one of the following: "parametric", "semiparametric", "cluster", "case", "residual"')
      }
      message(paste("Bootstrapping Procedure:", as.character(boot_type)))

      ######### Different weighted bootstrap methods to solve complex survey design
      if (boot_type %in% c('Rao-Wu-Yue-Beaumont', 'Preston', 'Rao-Wu') & random == TRUE & weighted == TRUE){
        B <- num_boots
        betaHat_boot <- betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))
        if(!silent)   print("Step 3.1: Bootstrap resampling-") #
        pb <- progress_bar$new(total = L)
        for(l in 1:L){
          pb$tick()
          data$Yl <- unclass(data[,out_index][,argvals[l]])

          ## Select single subjects with multiple observations
          nhats_design <- svydesign(ids = ~ 1,
                                    # strata = ~ strata,
                                    weights = ~ weight,
                                    data = data,
                                    nest=TRUE)

          m.tmp = suppressWarnings( svy2lme(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                          design = nhats_design,
                          return.devfun=TRUE) )

          rwyb_bootstrap <- nhats_design |> as_bootstrap_design(
            type = boot_type, replicates = B
          )

          a = suppressWarnings( boot2lme(m.tmp, rwyb_bootstrap) )

          betaTilde_boot[,l,] = t(a$beta)
        }
      }

      ## Finite Population Bootstrap (FPB)
      #### Not consider random effects (FPB)
      if (boot_type == 'FPB' & random == FALSE & weighted == TRUE){
        B <- num_boots
        betaHat_boot <- betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))
        if(!silent)   print("Step 3.1: Bootstrap resampling-") # , as.character(boot_type)
        pb <- progress_bar$new(total = L)
        for(l in 1:L){
          pb$tick()
          data$Yl <- unclass(data[,out_index][,argvals[l]])
          nhats_design <- svydesign(ids = ~ psu,
                                    strata = ~ strata,
                                    weights = ~ weight,
                                    data = data,
                                    nest=TRUE)

          betaTilde_boot[,l,] = suppressWarnings( t(fpb_bootstrap(nhats_design,
                                                R = B,
                                                formula = as.formula(paste0("Yl ~ ", model_formula[3])),
                                                random = FALSE)) )
        }
      }

      #### Consider random effects (FPB)
      if (boot_type == 'FPB' & random == TRUE & weighted == TRUE){
        B <- num_boots
        betaHat_boot <- betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))
        if(!silent)   print("Step 3.1: Bootstrap resampling-") # , as.character(boot_type)
        pb <- progress_bar$new(total = L)
        for(l in 1:L){
          pb$tick()
          data$Yl <- unclass(data[,out_index][,argvals[l]])
          nhats_design <- svydesign(ids = ~ psu,
                                    strata = ~ strata,
                                    weights = ~ weight,
                                    data = data,
                                    nest=TRUE)

          betaTilde_boot[,l,] = suppressWarnings( t(fpb_bootstrap(nhats_design,
                                                                  R = B,
                                                                  formula = as.formula(paste0("Yl ~ ", model_formula[3])),
                                                                  random = FALSE)) )
        }
      }

      if (boot_type %in% c('reb', "residual", "semiparametric",
                           "wild", "reb", "case") & random == TRUE & weighted == FALSE){

        B <- num_boots
        betaHat_boot <- betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))
        if(!silent)   print("Step 3.1: Bootstrap resampling-") # , as.character(boot_type)
        pb <- progress_bar$new(total = L)
        for(l in 1:L) {
          pb$tick()
          data$Yl <- unclass(data[,out_index][,argvals[l]])
          fit_uni <- suppressMessages(
            lmer(
              formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
              data = data,
              control = lmerControl(
                optimizer = "bobyqa", optCtrl = list(maxfun = 5000))
            )
          )

          # set seed to make sure bootstrap replicate (draws) are correlated
          # across functional domains
          set.seed(seed)

          fit_uni <- suppressMessages(
            lmer(
              formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
              data = data,
              control = lmerControl(optimizer = "bobyqa")
            ) )

          if (boot_type == "residual") {
            # for residual bootstrap to avoid singularity problems
            boot_sample <- lmeresampler::bootstrap(
              model = fit_uni,
              B = B,
              type = boot_type,
              rbootnoise = 0.0001
            )$replicates
            betaTilde_boot[,l,] <- t(as.matrix(boot_sample[,1:nrow(betaHat)]))
          }else if (boot_type %in% c("wild", "reb", "case") ) {
            # for case
            flist <- lme4::getME(fit_uni, "flist")
            re_names <- names(flist)
            clusters_vec <- c(rev(re_names), ".id")
            # for case bootstrap, we only resample at first (subject level)
            # because doesn't make sense to resample within-cluster for
            # longitudinal data
            resample_vec <- c(TRUE, rep(FALSE, length(clusters_vec) - 1))
            boot_sample <- lmeresampler::bootstrap(
              model = fit_uni,
              B = B,
              type = boot_type,
              resample = resample_vec, # only matters for type = "case"
              hccme = "hc2", # wild bootstrap
              aux.dist = "mammen", # wild bootstrap
              reb_type = 0
            )$replicates # for reb bootstrap only
            betaTilde_boot[, l, ] <- t(as.matrix(boot_sample[,1:nrow(betaHat)]))
          } else {
            use.u <- ifelse(boot_type == "semiparametric", TRUE, FALSE)
            betaTilde_boot[,l,] <- t(
              lme4::bootMer(
                x = fit_uni, FUN = function(.) {fixef(.)},
                nsim = B,
                seed = seed,
                type = boot_type,
                use.u = use.u
              )$t
            )
          }
        }
      }

      ## 'Rao-Wu-Yue-Beaumont'; 'Preston'
      if (boot_type %in% c('Rao-Wu-Yue-Beaumont', 'Preston') & random == FALSE & weighted == TRUE){
        B <- num_boots
        betaHat_boot <- betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))
        if(!silent)   print("Step 3.1: Bootstrap resampling-") #
        pb <- progress_bar$new(total = L)
        for(l in 1:L){
          pb$tick()
          data$Yl <- unclass(data[,out_index][,argvals[l]])
          nhats_design <- svydesign(ids = ~ psu,
                                    strata = ~ strata,
                                    weights = ~ weight,
                                    data = data,
                                    nest=TRUE)

          rwyb_bootstrap <- nhats_design |> as_bootstrap_design(
            type = boot_type, replicates = B
          )

          coef.output = suppressWarnings( svyglm(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                                 design = rwyb_bootstrap,
                                                 return.replicates=TRUE,
                                                 control = glm.control(maxit = 5000)) )
          betaTilde_boot[,l,] = coef.output$replicates[,2]
        }
      }

      ## Bootstrap reference (consider complex survey design)
      if (boot_type == 'ref_design' & random == FALSE){
        B <- num_boots
        betaHat_boot <- betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))
        if(!silent)   print("Step 3.1: Bootstrap resampling-") # , as.character(boot_type)
        pb <- progress_bar$new(total = L)
        for(l in 1:L){
          pb$tick()
          data$Yl <- unclass(data[,out_index][,argvals[l]])

          for (num in 1:B){
            set.seed(num)
            index = sample(row_number(data), dim(data)[1], replace = TRUE)
            dat.tmp = data[index, ]

            new_design <- svydesign(ids = ~ psu,
                                    strata = ~ strata,
                                    weights = ~ weight,
                                    data = dat.tmp,
                                    nest=TRUE)

            coef.tmp = suppressWarnings( coef(svyglm(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                                     design = new_design,
                                                     control = glm.control(maxit = 5000))) )

            betaTilde_boot[,l,num] = coef.tmp
          }
        }
      }

      ## Bootstrap reference (Not consider complex survey design)
      if (boot_type == 'ref_noDesign' & random == FALSE){
        B <- num_boots
        betaHat_boot <- betaTilde_boot <- array(NA, dim = c(nrow(betaHat), ncol(betaHat), B))
        if(!silent)   print("Step 3.1: Bootstrap resampling-") # , as.character(boot_type)
        pb <- progress_bar$new(total = L)
        for(l in 1:L){
          pb$tick()
          data$Yl <- unclass(data[,out_index][,argvals[l]])

          for (num in 1:B){
            set.seed(l*num)
            index = sample(row_number(data), dim(data)[1], replace = TRUE)
            dat.tmp = data[index, ]
            coef.tmp = suppressWarnings( coef(lm(formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
                                                 data=dat.tmp)) )

            betaTilde_boot[,l,num] = coef.tmp
          }
        }
      }

      # smooth across functional domain
      if(!silent)   print("Step 3.2: Smooth Bootstrap estimates")
      for(b in 1:B){
        betaHat_boot[,,b] <- t(apply(betaTilde_boot[,,b], 1, function(x) gam(x ~ s(argvals, bs = splines, k = (nknots + 1)), method = smooth_method)$fitted.values))
      }

      # Obtain bootstrap variance
      betaHat.var <- array(NA, dim = c(L,L,nrow(betaHat)))
      for(r in 1:nrow(betaHat)){
        betaHat.var[,,r] <- 1.2*var(t(betaHat_boot[r,,])) ## account for within-subject correlation
      }

      # Obtain qn to construct joint CI using the fast approach
      if(!silent)   print("Step 3.3")
      qn <- rep(0, length = nrow(betaHat))
      N <- 10000 ## sample size in simulation-based approach
      set.seed(seed) # set seed to make sure bootstrap replicate (draws) are correlated across functional domains
      for(i in 1:length(qn)){
        est_bs <- t(betaHat_boot[i,,])
        fit_fpca <- suppressWarnings( refund::fpca.face(est_bs, knots = nknots_fpca) ) # suppress sqrt(Eigen$values) NaNs
        ## extract estimated eigenfunctions/eigenvalues
        phi <- fit_fpca$efunctions
        lambda <- fit_fpca$evalues
        K <- length(fit_fpca$evalues)

        ## simulate random coefficients
        theta <- matrix(stats::rnorm(N*K), nrow=N, ncol=K) # generate independent standard normals
        if(K == 1){
          theta <- theta * sqrt(lambda) # scale to have appropriate variance
          X_new <- tcrossprod(theta, phi) # simulate new functions
        }else{
          theta <- theta %*% diag(sqrt(lambda)) # scale to have appropriate variance
          X_new <- tcrossprod(theta, phi) # simulate new functions
        }
        x_sample <- X_new + t(fit_fpca$mu %o% rep(1,N)) # add back in the mean function
        Sigma_sd <- Rfast::colVars(x_sample, std = TRUE, na.rm = FALSE) # standard deviation: apply(x_sample, 2, sd)
        x_mean <- colMeans(est_bs)
        un <- rep(NA, N)
        for(j in 1:N){
          un[j] <- max(abs((x_sample[j,] - x_mean)/Sigma_sd))
        }
        qn[i] <- stats::quantile(un, 0.95)
      }

      if(!silent)  message("Complete! \n -Use plot_fui() function to plot estimates \n -For more information, run the command:  ?plot_fui")

      resids = NA
      return(list(betaHat = betaHat, betaHat.var = betaHat.var, qn = qn,
                  aic = AIC_mat, residuals = resids, bootstrap_samps = B,
                  argvals = argvals))
    }
  }
}
