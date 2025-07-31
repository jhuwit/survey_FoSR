library(tidyverse)
library(tidyfun)
library(patchwork)
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "utils.R"))
force = FALSE
I = 10e5
L = 50
scenario = 1
family = "gaussian"
seed =4574
num_strata = 30
strata_sigma = 0.05
psu_sigma = 0.025
sd_beta = 0.125 # strata scaling factor
snr_b = 1 # signal noise ratio for random to fixed effects
snr_eps = 1
generate_superpopulation = function(I = 10e6, # size of superpopulation
                                    L = 50, # length of functional domain
                                    scenario = 1, # shape of function
                                    family = "gaussian",
                                    seed = 4574,
                                    num_strata = 30, # num total strata
                                    strata_sigma = 0.05,
                                    psu_sigma = 0.025,
                                    sd_beta = 0.125, # strata scaling factor
                                    snr_b = 1, # signal noise ratio for random to fixed effects
                                    snr_eps = 1 # signal to noise for gaussian
){
  stopifnot("scenario must be either 1 or 2" = scenario %in% c(1, 2))
  stopifnot("family must be either 'gaussian', 'poisson', or 'binomial'" = family %in% c("gaussian", "poisson", "binomial"))


  X_des <- cbind(1, rnorm(I, 0, 2))

  ## simulate true beta
  grid <- seq(0, 1, length = L)
  beta_fixed <- matrix(NA, 2, L)
  if (scenario == 1){
    beta_fixed[1, ] <- -0.15 - 0.1 * sin(2 * pi * grid) - 0.1 * cos(2 * pi * grid)
    beta_fixed[2, ] <- dnorm(grid, 0.6, 0.15) / 20
    # beta_fixed[2, ] <- dnorm(grid, 0.6, 0.15)
  } else if(scenario == 2) {
    beta_fixed[1,] <- 0.53 + 0.06*sin(3*grid*pi) - 0.03*cos(6.5*grid*pi)
    beta_fixed[2,] <- dnorm(grid, 0.2, .1)/60 + dnorm(grid, 0.35, .1)/200 -
      dnorm(grid, 0.65, .06)/250 + dnorm(grid, 1, .07)/60
  }
  rownames(beta_fixed) <- c("Intercept", "x")

  ## assign individuals to PSU and strata
  set.seed(seed)
  dirichlet_probs <- gtools::rdirichlet(1, rep(4, num_strata))
  set.seed(seed)
  stratum_assignments <- sample(1:num_strata, I, replace = TRUE, prob = dirichlet_probs)
  psu_assignments <- rep(NA, I)

  for (s in 1:num_strata) {
    set.seed(seed + s)
    num_in_strata = sum(stratum_assignments == s)
    num_psu = round(runif(1, 75, 125), 0)
    set.seed(seed + s)
    dps = gtools::rdirichlet(1, rep(10, num_psu))
    set.seed(seed + s)
    psu_in_stratum = sample(1:num_psu,
                            num_in_strata,
                            replace = TRUE,
                            prob = dps)
    psu_assignments[stratum_assignments == s] <- paste0(s, "_", psu_in_stratum)
  }

  ## create psu and strata-specific random effects
  nbasis <- 5
  basis <- fda::create.bspline.basis(c(0, 1), nbasis)
  Phi <- fda::eval.basis(grid, basis)

  set.seed(seed)
  strata_scores <- matrix(rnorm(num_strata * nbasis, 0, strata_sigma), num_strata, nbasis)
  strata_random_effects <- strata_scores %*% t(Phi)

  total_psu = length(unique(psu_assignments))

  set.seed(seed)
  psu_scores <- matrix(
    rnorm(total_psu * nbasis, 0, psu_sigma),
    total_psu,
    nbasis
  )
  psu_random_effects <- psu_scores %*% t(Phi)


  strata_effects_indiv <- strata_random_effects[stratum_assignments, ]
  psu_effects_indiv <- psu_random_effects[as.numeric(factor(psu_assignments)), ]
  random_effects <- strata_effects_indiv + psu_effects_indiv

  rm(strata_effects_indiv, psu_effects_indiv)

  ## add stratum-specific slope modifications
  set.seed(seed)
  stratum_scaling <- rnorm(num_strata, mean = 1, sd = sd_beta)

  beta1_by_stratum <- matrix(rep(stratum_scaling, each = L), nrow = num_strata) *
    matrix(rep(beta_fixed[2, ], times = num_strata),
           nrow = num_strata,
           byrow = TRUE)

  # assign to individuals
  beta1_by_indiv <- beta1_by_stratum[stratum_assignments, ]

  # adjust random effect based on signal to noise parameters

  fixef_signal <- matrix(rep(beta_fixed[1, ], I), nrow = I, byrow = TRUE) +
    X_des[, 2] * matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)

  # include stratum-specific slope variation in the random effects
  slope_re <- (stratum_scaling[stratum_assignments] - 1) *
    matrix(rep(beta_fixed[2, ], I), nrow = I, byrow = TRUE)
  ranef <- slope_re + random_effects
  ranef <- sd(fixef_signal) / sd(ranef) / snr_b * ranef
  rm(random_effects)


  # generate outcomes
  set.seed(seed)
  if (family == "gaussian") {
    sd_lp = sd(as.vector(t(fixef_signal)) + as.vector(t(ranef)))
    sigma = sd_lp / snr_eps

    Y_obs =  matrix(
      rnorm(
        n = I * L,
        mean = as.vector(t(fixef_signal)) + as.vector(t(ranef)),
        sd = sigma
      ),
      nrow = I,
      ncol = L,
      byrow = TRUE
    )
  } else if(family == "binomial") {
    p_true = plogis(as.vector(fixef_signal) + as.vector(t(ranef)))
    Y_obs <- matrix(
      rbinom(
        n = I * L,
        size = 1,
        prob = as.vector(t(p_true))
      ),
      nrow = I,
      ncol = L,
      byrow = TRUE
    )
  } else if (family == "poisson"){
    lam_true = exp(as.vector(fixef_signal) + as.vector(t(ranef)))
    Y_obs <- matrix(
      rpois(n = I * L,
            lambda = as.vector(t(lam_true))),
      nrow = I,
      ncol = L,
      byrow = TRUE
    )
  }

  return(list(Y_obs = Y_obs,
              X_des = X_des,
              stratum_assignments = stratum_assignments,
              psu_assignments = psu_assignments,
              dirichlet_probs = dirichlet_probs,
              beta_true = beta_fixed))
}

sample_from_population_wor <- function(X_des, # design matrix
                                       Y_obs, # y matrix
                                       I_n = 500, # subjects in each psu-strata combination
                                       num_strata = 30,  # total strata
                                       stratum_assignments, # assignment to ea strata
                                       num_selected_psu = 2,
                                       dirichlet_probs, # stratum probabilities
                                       psu_assignments, # assignment to psu (w/in strata)
                                       L = 50, # length of fnl domain
                                       seed = 1,
                                       inf_level = .5,
                                       compression = 3
){
  I = nrow(X_des)
  X1 = X_des[, 2]
  # select strata (use all for now)
  selected_strata = 1:num_strata
  p_strata_design <- dirichlet_probs[selected_strata]

  final_sample <- c() # to store final sample ids
  p1 <- rep(NA, I) # stage 1 selection probability
  p2 <- rep(NA, I) # stage 2 selection probability
  p_overall <- rep(NA, I) # overall selection probability
  psus <- c() # to store PSUs

  # within each strata, select PSU with replacement using PPS
  for (strata in 1:num_strata) {
    # strata = 1
    # individuals in the strata
    inds_in_stratum <- which(stratum_assignments == strata)

    # Get PSU sizes in this stratum
    psu_sizes <- table(psu_assignments[inds_in_stratum])
    psu_ids <- names(psu_sizes)

    # Sample PSUs WITH replacement using PPS
    set.seed(strata + seed)
    selected_psus <- sample(psu_ids,
                            size = num_selected_psu,
                            replace = FALSE,
                            prob = psu_sizes)

    psu_probs <- psu_sizes / sum(psu_sizes)  # PPS
    # probability of selection is 1 - (p(not selected both times))
    # psu_prob_selected <- 1 - (1 - psu_probs) ^ num_selected_psu
    # names(psu_prob_selected) <- psu_ids

    psu_prob_selected = map_dbl(.x = match(selected_psus, psu_ids),
                                .f = get_p_i,
                                psu_probs)

    names(psu_prob_selected) <- selected_psus
    # within each selected PSU select individuals based on X1
    for (psu in selected_psus) {
      inds_in_psu <- which(psu_assignments == psu &
                             stratum_assignments == strata)

      if (inf_level == 0) {
        inclusion_probs = rep(1 / length(inds_in_psu), length(inds_in_psu)) # uniform sampling
      } else {

        #
        # incl_score <- X1[inds_in_psu]  # Or your combination like: -0.6 * score1 + ...
        #
        # # Step 1: Center and scale
        # score_raw <- scale(incl_score)
        #
        # # Step 2: Compress range to avoid extreme logits
        # score_compressed <- score_raw * inf_level
        # score_compressed <- pmax(pmin(score_compressed, 3), -3)

        # Step 3: Apply plogis
        # inclusion_probs <- plogis(score_compressed)

        incl_score = rowMeans(Y_obs[inds_in_psu,]) * inf_level
        score_compressed <- pmax(pmin(incl_score, compression), -1 * compression)

        inclusion_probs = plogis(score_compressed)

      }
      # Informative sampling based on X1 or PC scores

      inclusion_probs <- inclusion_probs / sum(inclusion_probs) * I_n
      inclusion_probs[inclusion_probs > 1] <- 1

      #
      # inclusion_probs <- plogis(incl_score)
      # inclusion_probs <- inclusion_probs / sum(inclusion_probs)
      # inclusion_probs <- inclusion_probs * I_n
      # inclusion_probs[inclusion_probs > 1] <- 1
      set.seed(strata + seed + which(selected_psus == psu)) # ensure reproducibility
      sampled_units <- inds_in_psu[rbinom(length(inds_in_psu), 1, inclusion_probs) == 1]

      final_sample <- c(final_sample, sampled_units)
      psus <- c(psus, rep(psu, length(sampled_units)))
      p_psu <- psu_prob_selected[which(names(psu_prob_selected) == psu)]

      p1[inds_in_psu] <- p_psu
      p2[inds_in_psu] <- inclusion_probs
      p_overall[inds_in_psu] <- p_psu * inclusion_probs
    }
  }


  survey_weights <- 1 / p_overall


  dat.sim <- data.frame(
    ID = final_sample,
    X = X1[final_sample],
    strata = stratum_assignments[final_sample],
    psu = sub(".*\\_", "", psus),
    weight = survey_weights[final_sample],
    p_stage1 = p1[final_sample],
    p_stage2 = p2[final_sample]
  )


  Y_sample <- data.frame(Y_obs[final_sample, ])
  colnames(Y_sample) <- paste0("Y", 1:L)

  data =  cbind(dat.sim, Y_sample)
  return(data)
}

strata_sigma = 0.05
psu_sigma = sqrt(strata_sigma ^ 2 * 0.5)
scen = 1
snr_b = 1
snr_eps = 1
n = 500
get_p_i = function(i, probs) probs[i] * (1 + sum((probs[-i]) / (1-probs[-i])))

ifold = get_fold()



plot_results = function(scen, snr_b, snr_eps, sd_beta= 0.125,  inf_level = 0.5, n = 500, seed = 1,
                        compression = 3){
  lst = generate_superpopulation(
    scenario = scen,
    family = "gaussian",
    I = 10e6,
    L = 50,
    psu_sigma = psu_sigma,
    snr_b = snr_b,
    snr_eps = snr_eps,
    sd_beta = sd_beta
  )

  data = sample_from_population_wor(
    X_des = lst$X_des,
    Y_obs = lst$Y_obs,
    L = 50,
    I_n = n,
    num_strata = 30,
    stratum_assignments = lst$stratum_assignments,
    psu_assignments = lst$psu_assignments,
    dirichlet_probs = lst$dirichlet_probs,
    seed = seed,
    inf_level = inf_level,
    compression = compression
  )

  beta_true = lst$beta_true
  # rm(lst)
  betaTilde = map(
    .x = c("weighted", "unweighted"),
    .f = get_betatilde,
    data = data,
    family = "gaussian"
  )

  betaHat = map(.x = betaTilde, get_betahat, L = 50)
  # get bias
  bias1 = map_dbl(.x = betaHat,
                  .f = \(x) sum(abs(x[2,] - beta_true[2,])^2))
  bias0 = map_dbl(.x = betaHat,
                  .f = \(x) sum(abs(x[1,] - beta_true[1,])^2))
  bias1; bias0

  quants = quantile(data$weight, c(1/3, 2/3))

  Y_mat = data %>%
    select(starts_with("Y")) %>%
    as.matrix()


  mean_df =
    data %>%
    mutate(Y = I(Y_mat)) %>%
    mutate(quant = cut(weight, breaks = c(-Inf, quants, Inf), labels = c("low", "med", "high"))) %>%
    mutate(Y = matrix(Y, ncol = 50),
           Y_tf = tfd(Y, length = 50),
           Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
    group_by(quant) %>%
    summarize(Y_mean_strat = mean(Y_tf)) %>%
    ungroup() %>%
    mutate(Y_mean_strat = tf_smooth(Y_mean_strat, method = "lowess"))

  overall_df =
    data %>%
    mutate(Y = I(Y_mat)) %>%
    mutate(Y = matrix(Y, ncol = 50),
           Y_tf = tfd(Y, length = 50),
           Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
    summarize(Y_mean_strat = mean(Y_tf)) %>%
    mutate(Y_mean_strat = tf_smooth(Y_mean_strat, method = "lowess"))

  p1 = mean_df %>%
    ggplot() +
    geom_spaghetti(aes(y = Y_mean_strat, color = quant), linewidth = 1.1, alpha = 1) +
    scale_color_discrete(name = "Weight tertile") +
    labs(x = "Functional domain", y = "Y") +
    geom_spaghetti(data = overall_df, aes(y= Y_mean_strat), color = "black", linetype = "dashed") +
    theme_light()

  btt = beta_true %>%
    t() %>%
    as.data.frame() %>%
    magrittr::set_colnames(c("(Intercept)", "X")) %>%
    mutate(l = row_number()) %>%
    pivot_longer(cols = -l) %>%
    mutate(type = "truth")

  # beta_t = betaTilde %>%
  #   map(~ as.data.frame(t(.x))) %>%
  #   bind_rows(.id = "type") %>%
  #   group_by(type) %>%
  #   mutate(l = row_number()) %>%
  #   pivot_longer(cols = -c(type, l)) %>%
  #   mutate(type = if_else(type == 1, "weighted", "unweighted"))

  p2 =  betaHat %>%
    map(~ as.data.frame(t(.x))) %>%
    bind_rows(.id = "type") %>%
    group_by(type) %>%
    mutate(l = row_number()) %>%
    pivot_longer(cols = -c(type, l)) %>%
    mutate(type = if_else(type == 1, "weighted", "unweighted")) %>%
    ggplot(aes(x = l, y = value, color = type)) +
    labs(x = "Functional domain", y = "Coefficient") +
    geom_line() +
    facet_wrap(.~name, scales = "free_y") +
    geom_line(data = btt, aes(x = l, y = value, color = type), linetype = "dashed") +
    theme_light()  +
    scale_color_brewer(palette = "Dark2", name = "")

  rm(data)
  rm(Y_mat)
  rm(beta_true)
  rm(mean_df)
  gc()
  print(p1 + p2 + plot_layout(ncol = 1) +
          plot_annotation(title = paste0("Inf level: ", inf_level, ", n = ", n, ", snr_b = ", snr_b, ", snr_eps = ", snr_eps, ", sd_beta = ", sd_beta, ", scenario = ", scen),
                          subtitle = paste0("unwtd bias = ", round(bias0[2], 3), ", ", round(bias1[2], 3), " wtd bias = ", round(bias0[1], 4), ", ", round(bias1[1], 4))))
}

if(!dir.exists(here::here("simulation_debug"))){
  dir.create(here::here("simulation_debug"))
}

settings = expand_grid(scen = 1:2,
                       snr_b = c(.5, 1, 1.5, 5),
                       snr_eps = c(.5, 1, 1.5, 5),
                       sd_beta = .125,
                       inf_level = c(1, 5, 10),
                       compression = c(3, 4),
                       n = c(250, 500),
                       seed = c(1, 2)) %>%
  mutate(rn = row_number(),
         fold = ceiling(rn / 100))

settings = settings %>%
  filter(fold == ifold)

outfile = here::here("simulation_debug", paste0("plotv2_", sprintf("%02d", ifold), ".pdf"))
if(!file.exists(outfile) || force){
  pdf(outfile)

  for(row in 1:nrow(settings)){
    temp = settings[row,]
    plot_results(scen = temp$scen,
                 snr_b = temp$snr_b,
                 snr_eps = temp$snr_eps,
                 sd_beta = temp$sd_beta,
                 inf_level = temp$inf_level,
                 n = temp$n,
                 compression = temp$compression,
                 seed = temp$seed)
  }

  dev.off()
}
