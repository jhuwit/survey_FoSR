required_packages <- c("svrep", "purrr", "progress", "svylme", 'dplyr', 'readr')

install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package)
    library(package, character.only = TRUE)
  }
}
invisible(lapply(required_packages, install_and_load))
## When the PSU or individual sampling ratio is extreme
## the performance of svy is poor, which can be a breakthrough point.


# Generate Population Data
set.seed(209)
## A total of 10 strata, each stratum has 100 PSUs, each PSU has 10 individuals.

options(survey.lonely.psu = "adjust")

# One simulated sample and model fitting
one.simeff_fixed <- function(pop) {

  pop <- data.frame(
    x = rnorm(10000, s = 1),
    z = rep(rgamma(2000, 2), each = 5), # Covariate should be the same for each subject
    id = rep(1:2000, each = 5),
    g = rep(1:100, each = 100), # Define a total of 100 PSUs, each with 20 individuals, each individual has 5 observations
    strat = rep(1:10, each = 1000) # Define 10 strata, each stratum has 10 PSUs, each PSU has 20 individuals, each individual has 5 observations
  )

  # A total of 10 strata, each stratum has 100 PSUs, each PSU has 10 individuals.
  pop$u <- rnorm(2000, s = 0.7)[pop$id] # Random effect at the individual level
  # pop$x <- with(pop, x + rnorm(100, mean = 1, s = 2)[g])
  pop$y <- with(pop, 5*x + z
                + rnorm(10000, s = 1.5) # random error
                # + rnorm(100, s = 2)[g]
                + rnorm(10, s = 3)[strat]
                + rnorm(2000, s = 0.7)[id] # random intercept
                # + u * x
                ) # random slope

  m0 <- lme4::lmer(y ~ (1|id) + x + z,
                   data = pop)

  # First-stage sampling: select 40 PSUs from each stratum
  # v <- rep(5, 10); names(v) <- 1:10
  v <- c(rep(4,5), rep(4,5)); names(v) <- 1:10
  psus <- survey::stratsample(rep(1:10, each = 10), v)
  # Draw a specified number of PSUs (Primary Sampling Units) from each stratum.
  # Generate a vector of length 100, corresponding to 10 strata (strat 1 to 10), each containing 10 PSUs.
  # rep(4, 10) means drawing 4 PSUs from each of the 10 strata.

  stage1 <- subset(pop, g %in% psus)
  # When using equal probability sampling, the total number of individuals within a PSU does not need to be explicitly defined,
  # because it is assumed that individuals within a PSU are sampled with equal probability.

  # Second-stage sampling: further sample individuals within the selected PSUs.
  # Define the number of individuals drawn from each PSU.
  n <- as.vector(sapply(unique(stage1$g), function(gi) ifelse(gi < 50, 10, 10)))

  # Draw complete individuals from each PSU.
  in2 <- mapply(function(gi, ni) {
    ids <- unique(stage1$id[stage1$g == gi])  # Select all individual IDs within the current PSU.
    sampled_ids <- sample(ids, ni)            # Draw ni individual IDs.
    which(stage1$id %in% sampled_ids)         # Return all observation rows of these individuals.
  }, unique(stage1$g), n)

  # Merge sampling results and retain complete observations.
  stage2 <- stage1[unlist(in2), ]

  # Generate weights
  # stage2$w1 = with(stage2, 10 / 4) # There are 100 PSUs in a stratum, and 40 are sampled.
  stage2$w1 = with(stage2, ifelse(strat <= 5, 10/4, 10/4))
  stage2$w2 <- with(stage2, ifelse(g < 50, 20 / 10, 20 / 10))
  stage2$w <- stage2$w1 * stage2$w2
  stage2$wid <- 1

  stage2$n1 <- rep(10, nrow(stage2)) # Total number of PSUs in each stratum (strat).
  stage2$n2 <- rep(10, nrow(stage2)) # Total number of individuals in each PSU.

  # Fit LMER
  m1 <- lme4:::lmer(y~(1|id)+x+z,data=stage2)
  m1.coef = lme4::fixef(m1)[2]
  m1.se = confint(m1, method = "Wald")["x",]; names(m1.se) = c('lower', 'upper')

  # Fit survey design weighted regression model
  des <- survey::svydesign(id = ~g+id,
                           strata = ~strat,
                           weights = ~w1+w2,
                           # fpc = ~n1 + n2,
                           data = stage2)

  m2 <- svylme::svy2lme(y ~ (1|id) + x + z,
                        design = des,
                        method = 'nested',
                        return.devfun=TRUE)
  m2.coef = m2$beta[2]
  m2.se <- confint(m2)['x',]; names(m2.se) = c('lower', 'upper')

  boot.num = 100
  # Non-parametric Bootstrap
  nonparametric_bootstrap <- function(stage2, boot.num = 100) {
    # replicate(boot.num, {
    #   boot_sample <- stage2 %>%
    #     group_by(strat, g) %>%
    #     slice_sample(prop = 1, replace = TRUE)
    #   boot_des <- svydesign(id = ~g+id,
    #                         strata = ~strat,
    #                         weights = ~w1+w2,
    #                         # fpc = ~n1 + n2,
    #                         data = boot_sample, nest = TRUE)
    #   m.tmp = suppressWarnings( svylme::svy2lme(y ~ (1|id) + x + z,
    #                                             method = 'nested',
    #                                             design= boot_des) )
    #   m.tmp$beta[2]
    # })

    replicate(boot.num, {
      boot_id = sample(unique(stage2$id), length(unique(stage2$id)), replace = TRUE)
      # row_num = c()
      # for (i in 1:length(boot_id)) {
      #   num = which(stage2$id == boot_id[i])
      #   row_num = c(row_num, num)
      # }
      row_num <- unlist(lapply(boot_id, function(id) which(stage2$id == id)))
      boot_sample = stage2[row_num,]
      boot_des <- svydesign(id = ~g+id,
                            strata = ~strat,
                            weights = ~w1+w2,
                            # fpc = ~n1 + n2,
                            data = boot_sample, nest = TRUE)
      m.tmp = suppressWarnings( svylme::svy2lme(y ~ (1|id) + x + z,
                                                method = 'nested',
                                                design= boot_des) )
      m.tmp$beta[2]
    })
  }

  # Rao-Wu/Preston Bootstrap
  bootstrap_fit <- function(des, type, boot.num = 100) {
    boot_design <- des |> as_bootstrap_design(type = type,
                                              replicates = boot.num)
    m3.coef = boot2lme(m2, boot_design)
    m3.coef$beta[,2]
  }

  # Bootstrap Results
  np_results <- nonparametric_bootstrap(stage2,
                                        boot.num = boot.num)
  rwyb_results <- bootstrap_fit(des, "Rao-Wu-Yue-Beaumont",
                                boot.num = boot.num)
  preston_results <- bootstrap_fit(des, "Preston",
                                   boot.num = boot.num)

  # Return results as list
  res = list(
    sample = c(m1.coef,
               upper = unname(m1.se[2]),
               lower = unname(m1.se[1]) ),
    svy = c(coef(m2),
            upper = unname(m2.se[2]),
            lower = unname(m2.se[1]) ),
    svy_np = c(coef(m2), se = sd(np_results),
               lower = unname( quantile(np_results, 0.025) ),
               upper = unname(  quantile(np_results, 0.975)) ),
    svy_rwyb = c(coef(m2), se = sd(rwyb_results),
                 lower = unname( quantile(rwyb_results, 0.025) ),
                 upper = unname( quantile(rwyb_results, 0.975)) ),
    svy_preston = c(coef(m2), se = sd(preston_results),
                    lower = unname( quantile(preston_results, 0.025) ),
                    upper = unname( quantile(preston_results, 0.975) ))
  )

  # Return simulated metrics
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

# 3.Simulation Results
sim.num <- 500
set.seed(209)
pb <- progress_bar$new(
  format = "Simulation [:bar] :current/:total (:percent) eta: :eta",
  total = sim.num,
  clear = FALSE,
  width = 60
)

rval_fixed <- map(1:sim.num, ~ {
  pb$tick()
  one.simeff_fixed(pop)
})

results_df <- do.call(rbind, lapply(rval_fixed, function(res) {
  df <- as.data.frame(do.call(rbind, res$metrics))
  df$Method <- rownames(df)
  return(df)
}))

average_results <- aggregate(. ~ Method, data = results_df, mean)
print(knitr::kable(average_results, digits = 3))

# 4.Confidence Interval
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

readr::write_rds(rval_fixed, 'results/measurement_raw/single-random-S0.rds')
readr::write_csv(average_results, 'results/metrix_summary/single-random-S0.csv')
