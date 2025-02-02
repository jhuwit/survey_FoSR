fpb_bootstrap <- function(design,
                          R,
                          random = FALSE,
                          formula = formula) {
  N <- nrow(design$cluster)
  n <- design$fpc$sampsize

  replicates <- lapply(1:R, function(i) {
    # generate FPB weights
    new_weights <- rep(0, N)
    for (h in unique(design$strata)) {
      Nh <- sum(design$strata == h)
      nh <- sum(design$strata == h & design$prob > 0)
      mh <- rbinom(1, Nh, nh/Nh)
      new_weights[design$strata == h] <- rmultinom(1, mh, rep(1/nh, nh))
    }

    # update design objects
    new_design <- design
    new_design$prob <- 1 / (new_weights * design$prob)
    new_design$pweights <- 1 / new_design$prob

    # fit model and return coefficients
    if (random == FALSE){
      model <- svyglm(formula,
                      design = new_design, family = gaussian())
    } else {
      nhats_design <- svydesign(ids = ~ psu,
                                strata = ~ strata,
                                weights = ~ weight,
                                data = data,
                                nest=TRUE)

      model <- svy2lme(formula,
                      design = nhats_design,
                      method = "nested")
    }

    return(coef(model))
  })

  return(do.call(rbind, replicates))
}


boot.simple <- function(B = 500){
  sim_data = generate_nhats_like_data(n_samples = 100)
  beta_boot_point = array(NA, dim = c(2, B))
  for (num in 1:B){
    set.seed(num)
    index = sample(row_number(sim_data), dim(sim_data)[1], replace = TRUE)
    dat.tmp = sim_data[index,]

    coef.tmp = suppressWarnings( coef(lm(formula = stats::as.formula('Y ~ X'),
                                         data = dat.tmp)) )
    beta_boot_point[, num] = coef.tmp
  }
  return(beta_boot_point)
}

boot.simple.design <- function(B = 500){
  sim_data = generate_nhats_like_data(n_samples = 100)
  beta_boot_point = array(NA, dim = c(2, B))
  for (num in 1:B){
    set.seed(num)
    index = sample(row_number(sim_data), dim(sim_data)[1], replace = TRUE)
    dat.tmp = sim_data[index,]

    sim_design <- svydesign(ids = ~ psu,
                            strata = ~ strata,
                            weights = ~ weight,
                            data = dat.tmp,
                            nest=TRUE)

    coef.tmp = suppressWarnings( coef(svyglm(formula = stats::as.formula('Y ~ X'),
                                             design = sim_design,
                                             control = glm.control(maxit = 5000))) )
    beta_boot_point[, num] = coef.tmp
  }
  return(beta_boot_point)
}


boot.svydesign <- function(B = 100, boot_type = 'Rao-Wu'){
  sim_data = generate_nhats_like_data()
  beta_boot_point = array(NA, dim = c(2, B))

  ## simulated results
  sim_design <- svydesign(
    ids = ~ psu,
    strata = ~ strata,
    weights = ~ weight,
    data = sim_data,
    nest = TRUE
  )

  rwyb_bootstrap <- sim_design |> as_bootstrap_design(
    type = boot_type, replicates = B)

  coef.tmp = suppressWarnings( svyglm(formula = stats::as.formula('Y~X'),
                                      design = rwyb_bootstrap,
                                      return.replicates=TRUE,
                                      control = glm.control(maxit = 5000)) )
  beta_boot_point[,1:B] = coef.tmp$replicates[,2]
  return(beta_boot_point)
}

