required_packages <- c("svrep", "purrr", "progress", "svylme", 'dplyr', 'readr')

install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}

invisible(lapply(required_packages, install_and_load))

# Generate Population Data
set.seed(209)

# One simulated sample and model fitting
one.simeff_fixed <- function() {
  ## A total of 10 strata, each containing 100 PSU, each PSU has 20 individuals
  pop <- data.frame(
    x = rnorm(10000, s=1),
    z = rgamma(10000, 2),
    id = 1:10000,
    g = rep(1:1000, each = 10), # Define a total of 1000 PSU, each with 10 individuals
    strat = rep(1:10, each = 1000) # Define 10 strata, each containing 100 PSU
  )
  # pop$x <- with(pop, x + rnorm(1000, mean = 1, s = 2)[g])
  pop$y <- with(pop, 5*x + z
                + rnorm(10000, s = 1.5) # random error
                + rnorm(1000, s = 2)[g]
                + rnorm(10, s = 3)[strat]
  ) # random slope
  m0 <- lm(y ~ x + z, data = pop)  # Overall effect estimation

  # First stage sampling: select 40 PSU from each stratum
  v <- rep(40, 10); names(v) <- 1:10
  # v <- c(rep(60,5), rep(20,5)); names(v) <- 1:10
  psus <- survey::stratsample(rep(1:10, each = 100), v)
  # Extracts a specified number of PSUs (Primary Sampling Units) from each stratum.
  # Generate a vector of length 1000, which corresponds to 10 strata (strat 1 to 10), each containing 100 PSUs
  # rep(100, 10), means extract 10 PSUs from each stratum.
  # When rep(100, 10) is equivalent to selecting all, since we don't secondary sample individuals
  stage1 <- subset(pop, g %in% psus)
  # When sampling with equal probability, denotes the total number of individuals within the PSU,
  # which does not need to be explicitly defined, since it is assumed that individuals within the PSU are sampled with equal probability

  # Second stage sampling: further sample individuals within the selected groups
  n <- as.vector(sapply(unique(stage1$g), function(gi) ifelse(gi < 500, 5, 5)))
  in2 <- mapply(function(gi, ni) sample(which(stage1$g == gi), ni), unique(stage1$g), n)
  stage2 <- stage1[unlist(in2), ]

  # Generate weights
  stage2$w1 = with(stage2, 100 / 40) # 100 PSU per stratum, 40 sampled
  # stage2$w1 = with(stage2, ifelse(strat <= 5, 100/60, 100/20))
  stage2$w2 <- with(stage2, ifelse(g < 500, 10 / 5, 10 / 5))
  stage2$w <- stage2$w1 * stage2$w2

  stage2$n1 <- rep(100, nrow(stage2)) # Total number of PSU in each stratum
  stage2$n2 <- rep(10, nrow(stage2)) # Total number of individuals per PSU

  # Fit least squares regression
  m1 <- lm(y ~ x + z, data = stage2)
  se1 <- summary(m1)$coefficients["x", "Std. Error"]

  # Fit survey design weighted regression model
  des <- survey::svydesign(id = ~g+id,
                           strata = ~strat,
                           weights = ~w1+w2,
                           # fpc = ~n1 + n2, # fpc adjustment
                           data = stage2,
                           nest = TRUE)
  m2 <- survey::svyglm(y ~ x + z, design = des)
  se2 <- summary(m2)$coefficients["x", "Std. Error"]

  boot.num = 200
  # Non-parametric Bootstrap
  nonparametric_bootstrap <- function(stage2, boot.num = 50) {
    # replicate(boot.num, {
    #   boot_sample <- stage2 %>%
    #     group_by(strat, g) %>%
    #     slice_sample(prop = 1, replace = TRUE)
    #
    #   boot_des <- svydesign(id = ~g+id,
    #                         strata = ~strat,
    #                         weights = ~w1+w2,
    #                         # fpc = ~n1 + n2,
    #                         data = boot_sample, nest = TRUE)
    #   coef(suppressWarnings(svyglm(y ~ x + z, design = boot_des)))["x"]
    # })

    replicate(boot.num, {
      boot_id <- sample(1:dim(stage2)[1], dim(stage2)[1], replace = TRUE)
      boot_sample <- stage2[boot_id,]
      coef(suppressWarnings(lm(y ~ x + z,
                               weights = w,
                               data = boot_sample)))["x"]
    })
  }

  # Rao-Wu/Preston Bootstrap
  bootstrap_fit <- function(des, type, boot.num = 50) {
    boot_design <- des |> as_bootstrap_design(type = type, replicates = boot.num)
    sim_data <- as_data_frame_with_weights(boot_design, full_wgt_name = "FULL_SAMPLE_WGT", rep_wgt_prefix = "REP_WGT_")
    sapply(1:boot.num, function(num) {
      dat.tmp <- stage2
      dat.tmp$weight <- sim_data[[paste0('REP_WGT_', num)]]
      new_design <- svydesign(ids = ~g+id,
                              strata = ~strat,
                              weights = ~weight,
                              # fpc = ~n1 + n2,
                              data = dat.tmp, nest = TRUE)
      coef(svyglm(y ~ x + z, design = new_design))["x"]
    })
  }

  # Bootstrap Results
  np_results <- nonparametric_bootstrap(stage2, boot.num = boot.num)
  rwyb_results <- bootstrap_fit(des, "Rao-Wu-Yue-Beaumont", boot.num = boot.num)
  preston_results <- bootstrap_fit(des, "Preston", boot.num = boot.num)
  ## Warnings: observations with zero weight not used for calculating dispersion

  # Return results as list
  res = list(
    sample = c(coef(m1), se = se1,
               lower = confint(m1)[2,1],
               upper = confint(m1)[2,2]),
    svy = c(coef(m2), se = se2,
            lower = confint(m2)[2,1],
            upper = confint(m2)[2,2]),
    svy_np = c(coef(m2), se = sd(np_results),
               lower = unname(quantile(np_results, 0.025)),
               upper = unname(quantile(np_results, 0.975))),
    svy_rwyb = c(coef(m2), se = sd(rwyb_results),
                 lower = unname(quantile(rwyb_results, 0.025)),
                 upper = unname(quantile(rwyb_results, 0.975))),
    svy_preston = c(coef(m2), se = sd(preston_results),
                    lower = unname(quantile(preston_results, 0.025)),
                    upper = unname(quantile(preston_results, 0.975)))
  )

  # Return results as list
  methods <- c("sample", "svy", "svy_rwyb", "svy_preston", "svy_np")

  res_metrics <- lapply(methods, function(method) {
    lower <- res[[method]]["lower"]
    upper <- res[[method]]["upper"]
    ci_width <- unname(upper - lower)
    # coverage <- unname(lower <= lme4::fixef(m0)[2] & upper >= lme4::fixef(m0)[2])
    coverage <- unname(lower <= 5 & upper >= 5)
    return(c(coverage = coverage, ci_width = ci_width))
  })
  names(res_metrics) <- methods

  return(list(results = res, metrics = res_metrics))
}

# 3.Save Simulation results
sim.num <- 500
set.seed(208)
pb <- progress_bar$new(
  format = "Simulation [:bar] :current/:total (:percent) eta: :eta",
  total = sim.num,
  clear = FALSE,
  width = 60
)

rval_fixed <- rerun(sim.num, {
  pb$tick()
  one.simeff_fixed()
})

results_df <- do.call(rbind, lapply(rval_fixed, function(res) {
  df <- as.data.frame(do.call(rbind, res$metrics))
  df$Method <- rownames(df)
  return(df)
}))

average_results <- aggregate(. ~ Method, data = results_df, mean)
print(knitr::kable(average_results, digits = 3))

# 4. Confiden Interval
methods <- c("sample", "svy", "svy_rwyb", "svy_preston", "svy_np")

calculate_coverage <- function(results, method) {
  coverage <- sapply(results, function(res) {
    lower <- res[['results']][[method]]["x"] - 1.96*res[['results']][[method]]["se"]
    upper <- res[['results']][[method]]["x"] + 1.96*res[['results']][[method]]["se"]
    ci_width <- upper - lower
    return(c( (lower <= 5 & upper >= 5),
              CI_width <- ci_width) )
  })
  return(list( coverage_rate = mean(coverage[1,]),
               mean_width = mean(coverage[2,]),
               sd_width = sd(coverage[2,])) )
}

detailed_results <- lapply(methods, function(m) {
  stats <- calculate_coverage(rval_fixed, m)
  data.frame(
    Method = m,
    Coverage_Rate = stats$coverage_rate,
    Mean_CI_Width = stats$mean_width,
    SD_CI_Width = stats$sd_width
  )
})

results_df <- do.call(rbind, detailed_results)

readr::write_rds(rval_fixed, 'results/measurement_raw/single-fixed-S2.rds')
readr::write_csv(average_results, 'results/metrix_summary/single-fixed-S2.csv')
