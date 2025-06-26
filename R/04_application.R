## work on real world data result
library(tidyverse)
steps_df = read_rds(here::here("data", "steps_covariates.rds")) %>%
  mutate(SEQN = as.character(SEQN))
source(here::here("R", "01_sim_functions.R"))

steps_df_small =
 steps_df %>%
   pivot_longer(cols = starts_with("min")) %>%
  mutate(num = as.numeric(sub(".*min\\_", "", name)),
         hr = floor(num / 30)) %>%
  group_by(hr, SEQN, data_release_cycle, gender, age_in_years_at_screening,
           full_sample_2_year_mec_exam_weight,
           race_hispanic_origin,masked_variance_pseudo_stratum,
           masked_variance_pseudo_psu
           ) %>%
  summarize(name = first(name),
            value = mean(value, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(race_cat =
           case_when(race_hispanic_origin == "Non-Hispanic White" ~ "NH White",
                     race_hispanic_origin == "Non-Hispanic Black" ~ "NH Black",
                     race_hispanic_origin %in% c("Other Hispanic", "Mexican American") ~ "Hispanic",
                     .default = "Other"))

steps_df_small =
  steps_df_small %>%
  mutate(race_cat = factor(race_cat, levels = c("NH White", "NH Black", "Hispanic", "Other")))


fit_types = c('weighted', 'weighted', 'weighted', 'unweighted')
boot_types = c('BRR', 'weighted', 'unweighted', 'unweighted')

data =
  steps_df_small %>%
  mutate(sex_bin = if_else(gender == "Female", 1, 0)) %>%
  mutate(age_cat = cut(age_in_years_at_screening, breaks=c(18, 30, 40, 50, 60, 70, 80), include.lowest = TRUE)) %>%
  rename(weight = full_sample_2_year_mec_exam_weight,
         psu = masked_variance_pseudo_psu,
         strata = masked_variance_pseudo_stratum) %>%
  rename_with(.cols = starts_with("min"), ~str_replace(., "min_", "Y"))


# betaTilde = map(.x = fit_types, .f = get_betatilde, data = data, family = "poisson")


get_betatilde = function(data, family = "gaussian", type = "weighted",
                         model_formula = as.formula(paste0('Y~', 'age_cat + sex_bin + race_cat'))){
  if(!(type %in% c("weighted", "unweighted"))){
    stop("please specify boot_type as one of 'weighted', 'unweighted'")
  }

  out_index <- grep(paste0("^", model_formula[2]), names(data)) # indices that start with the outcome name
  if(length(out_index) != 1){ # functional observations stored in multiple columns
    L <- length(out_index)
  }else{ # observations stored as a matrix in one column using the I() function
    L <- ncol(data[,out_index])
  }
  argvals <- 1:L

  unimm_wt <- function(l) {
    data$Yl <- unclass(data[, out_index][, l]) %>%  unlist() %>% unname() %>% round(., digits = 0)

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
    data$Yl <- unclass(data[, out_index][, l]) %>%  unlist() %>% unname() %>% round(., digits = 0)
    fit_uni <- suppressMessages(glm(
      formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
      family = family,
      control = glm.control(maxit = 5000),
      data = data
    ))
    betaTilde <- coef(fit_uni)
    return(list(betaTilde = betaTilde))
  }
  if (type == "weighted") {
    massmm <- lapply(argvals, unimm_wt)
  } else if (type == "unweighted") {
    massmm <- lapply(argvals, unimm)
  }

  # Obtain betaTilde, fixed effects estimates
  betaTilde <- t(do.call(rbind, lapply(massmm, '[[', 1)))
  colnames(betaTilde) <- argvals
  return(betaTilde)
}
betaTilde = map(.x = fit_types, .f = get_betatilde, data = data, family = "poisson",
                model_formula = as.formula(paste0('Y~', 'age_cat + sex_bin + race_cat')))
betaHat = map(.x = betaTilde, .f = get_betahat, L = 49)
# get_betahat(betaTilde[[1]], L = 49)


run_boots = function(data, boot_type, betaHat, family = "gaussian",
                     num_boots = 500, set_seed = TRUE, seed = 2025, L = 50,
                     samp_stages = NULL,
                     model_formula = as.formula(paste0('Y~', 'X'))){
  argvals = 1:L
  if(!(boot_type %in% c("BRR", "Rao-Wu-Yue-Beaumont", "weighted", "unweighted"))){
    stop("please specify boot_type as one of 'BRR', 'Rao-Wu-Yue-Beaumont', 'weighted', 'unweighted'")
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
    data$Yl <- unclass(data[, out_index][, l]) %>%  unlist() %>% unname() %>% round(., digits = 0)
    svy_design <- svydesign(ids = ~ psu,
                            strata = ~ strata,
                            weights = ~ weight,
                            data = data,
                            nest = TRUE)
    if(set_seed){
      set.seed(seed)
    }
    if(boot_type %in% c("Rao-Wu-Yue-Beaumont")){
      rwyb_wts = make_rwyb_bootstrap_weights(
        num_replicates = num_boots,
        samp_unit_ids = svy_design$cluster,
        strata_ids = svy_design$strata,
        samp_unit_sel_probs = matrix(c(data$p_stage1, data$p_stage2), byrow = FALSE, ncol = 2),
        samp_method_by_stage = samp_stages,
        allow_final_stage_singletons = TRUE,
        output = "weights"
      )

      bs_design <- svrepdesign(
        weights = ~weight,
        repweights = rwyb_wts,
        rscales = rep(1, num_boots),
        scale = 1,
        type = "other",
        data = data,
        combined.weights = TRUE
      )
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
    } else if (boot_type == "weighted"){
      coefs <- do.call(
        cbind,
        lapply(1:num_boots, function(num) {
          set.seed(l * num)
          index <- sample(seq_len(nrow(data)), dim(data)[1], replace = TRUE, prob = data$weight / sum(data$weight))
          dat.tmp <- data[index, ]

          model <- glm(
            formula = stats::as.formula(paste0("Yl ~ ", model_formula[3])),
            data = dat.tmp,
            weights = weight,
            family = family)
          model$coefficients
        }))
      betaTilde_boot[,l,] <- coefs
    } else if (boot_type == "unweighted"){
      coefs <- do.call(
        cbind,
        lapply(1:num_boots, function(num) {
          set.seed(l * num)
          index <- sample(seq_len(nrow(data)), dim(data)[1], replace = TRUE)
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

data = data %>%
  mutate(ID = SEQN)

library(future)
library(furrr)
cores = parallelly::availableCores() - 1
plan(multisession, workers = cores)

res = future_map2(
  .x = boot_types,
  .y = betaHat,
  .f = run_boots,
  data = data,
  family = "poisson",
  num_boots = 500,
  set_seed = TRUE,
  seed = 2025,
  L = 49,
  model_formula = as.formula(paste0('Y~', 'age_cat + sex_bin + race_cat')),
  .options = furrr_options(seed = TRUE)
)
cis = future_map2(
  .x = res,
  .y = betaHat,
  .f = get_cis,
  smooth_for_ci = TRUE,
  .options = furrr_options(seed = TRUE)
)

plan(sequential)

write_rds(cis, here::here("results", "application.rds"))
# res = map2(
#   .x = boot_types,
#   .y = betaHat,
#   .f = run_boots,
#   data = data,
#   family = "poisson",
#   num_boots = 10,
#   set_seed = TRUE,
#   seed = 2025,
#   L = 49,
#   model_formula = as.formula(paste0('Y~', 'age_cat + sex_bin + race_cat'))
#   # .options = furrr_options(seed = TRUE)
# )
# cis = map2(
#   .x = res,
#   .y = betaHat,
#   .f = get_cis,
#   smooth_for_ci = TRUE,
#   L = 49
#   # .options = furrr_options(seed = TRUE)
# )
#
# plan(sequential)


plot = FALSE
if(plot){
  get_plot_df = function(sim_res) {
    num_var = nrow(sim_res$betaHat)
    map_dfr(
      .x = 1:num_var,
      .f = function(r) {
        data.frame(
          s = 1:length(sim_res$betaHat[r, ]),
          beta = sim_res$betaHat[r, ],
          lower = sim_res$betaHat[r, ] - 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
          upper = sim_res$betaHat[r, ] + 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
          lower.joint = sim_res$betaHat[r, ] - sim_res$qn[r] *
            sqrt(diag(sim_res$betaHat.var[, , r])),
          upper.joint = sim_res$betaHat[r, ] + sim_res$qn[r] * sqrt(diag(sim_res$betaHat.var[, , r]))
        ) %>%
          mutate(name = sim_res$betaHat %>% rownames() %>% .[r])
      }
    )
  }
  plt_df_boot = map(cis, .f = get_plot_df)


  plt_df_boot %>%
    bind_rows(.id = "model") %>%
    ggplot(aes(x = s, y = beta, group = model)) +
    geom_line(aes(color = model))  +
    facet_wrap(.~name, scales = "free_y") +
    geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = model), alpha = .3) +
    geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = model), alpha = .1) +
    labs(x = "L", y = "Estimate") +
    theme_classic() +
    geom_hline(aes(yintercept = 0))

}
