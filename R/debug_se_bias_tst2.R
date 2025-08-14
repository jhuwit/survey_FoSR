library(future)
library(furrr)
library(tidyverse)
library(here)
library(devtools)
library(fastFMM)
library(haven)
library(dplyr)
library(survey)
library(progress)
library(lme4)
library(paletteer)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(tidyfun)
library(mvtnorm)
library(refund)
library(patchwork)
library(svrep)
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "00_data_gen_function_ff.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "create_settings.R"))

force = FALSE
force_plot = FALSE
nsim = 50
B = 500
ncores = parallelly::availableCores() - 1
fit_types = c('weighted', 'unweighted')
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'weighted', 'unweighted')
options(survey.lonely.psu = "adjust")

get_se = function(mod_output){
  mean_pw_se = apply(mod_output$betaHat.var, 3, function(mat) sqrt(diag(mat)))
  mean_pw_se %>%
    as_tibble() %>%
    magrittr::set_colnames(c("Intercept", "X")) %>%
    mutate(l = row_number()) %>%
    pivot_longer(cols = -l, names_to = "var")
}

## plot std error as function of domain

ifold = 25

temp = settings %>%
  filter(fold == ifold)

outfile = here::here("debug", "bias_se", paste0("fold_tst", sprintf("%03d", ifold), ".pdf"))
outfile2 = here::here("debug", "bias_se_dfs", paste0("fold_se_tst2", sprintf("%03d", ifold), ".rds"))
outfile3 = here::here("debug", "bias_se_dfs", paste0("fold_bias_tst2", sprintf("%03d", ifold), ".rds"))

# create partial dir to store ea. iteration

if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)
if (!dir.exists(dirname(outfile2))) dir.create(dirname(outfile2), recursive = TRUE)

biases = list()
ses = list()

if(!file.exists(outfile2) || force) {
  temp = settings[ifold, ] # set settings for this simulation

  lst = generate_superpopulation(
    scenario = temp$scenario,
    family = temp$family,
    I = 10e6,
    L = temp$len,
    snr_b = temp$snr_b,
    snr_eps = temp$snr_eps,
    strata_sigma = temp$strata_sigma,
    strata_scale = temp$strata_scale,
    seed = 111
  )

  for (iter in 1:nsim) {
    set.seed(iter)
    data = sample_from_population_wor(
      X_des = lst$X_des,
      Y_obs = lst$Y_obs,
      L = temp$len,
      I_n = temp$In,
      num_strata = 30,
      stratum_assignments = lst$stratum_assignments,
      psu_assignments = lst$psu_assignments,
      dirichlet_probs = lst$dirichlet_probs,
      seed = iter,
      inf_level = 20,
      compression = 2,
    )
    beta_true = lst$beta_true
    betaTilde = map(
      .x = fit_types,
      .f = get_betatilde,
      data = data,
      family = temp$family
    )

    betaHat = map(.x = betaTilde, get_betahat, L = temp$len)


    bias_df = map(
      betaHat,
      \(x) (x - beta_true) %>% as_tibble() %>% mutate(var = c("Intercept", "X")) %>% pivot_longer(cols = -var)
    ) %>%
      bind_rows(.id = "type") %>%
      mutate(type = if_else(type == "1", "weighted", "unweighted")) %>%
      mutate(l = as.numeric(name))

    biases[[iter]] <- bias_df

    res = future_map(
      .x = boot_types,
      .f = run_boots_fast,
      betaHat = betaHat[[1]],
      # don't actually use betahat in calculation so we can just use the weighted one
      data = data,
      family = temp$family,
      num_boots = B,
      seed = 2025,
      L = temp$len,
      samp_stages = c("PPSWOR", "Poisson"),
      .options = furrr_options(seed = TRUE)
    )

    fit_list = list(betaHat[[1]], betaHat[[1]], betaHat[[1]], betaHat[[2]])

    cis = future_map2(
      .x = res,
      .y = fit_list,
      .f = get_cis,
      L = temp$len,
      smooth_for_ci = TRUE,
      .options = furrr_options(seed = TRUE)
    )

    se = map(.x = cis, .f = get_se)
    se_df =
      se %>%
      bind_rows(.id = "boot_type") %>%
      mutate(
        boot_type =
          case_when(
            boot_type == "1" ~ "RWYB",
            boot_type == "2" ~ "BRR",
            boot_type == "3" ~ "Weighted",
            .default = "Unweighted"
          )
      )

    ses[[iter]] <- se_df
    rm(data, cis, res, bias_df, se_df)
    print(iter)
  }

  res_df =
    bind_rows(biases, .id = "iter")

  ses_df = bind_rows(ses, .id = "iter")

  write_rds(res_df, outfile3)
  write_rds(ses_df, outfile2)
  # plt = res_df %>%
  #   filter(var == "X") %>%
  #   ggplot(aes(x = l, y = value, color = type)) +
  #   geom_line() +
  #   facet_wrap(iter~., nrow = 10, ncol = 10) +
  #   geom_hline(aes(yintercept = 0))

  if(!file.exists(outfile) || force_plot){
    bias_plt =
      res_df %>%
      filter(var == "X") %>%
      ggplot(aes(x = l, y = value, color = type, fill = type)) +
      geom_smooth() +
      geom_hline(aes(yintercept = 0), linetype = "dashed") +
      labs(x = "Functional Domain", y = "Mean Bias") +
      scale_color_brewer(palette = "Dark2") +
      scale_fill_brewer(palette = "Dark2")

    se_plt =
      ses_df %>%
      filter(var == "X") %>%
      ggplot(aes(x = l, y = value, color = boot_type, fill = boot_type)) +
      geom_smooth() +
      facet_wrap(.~boot_type, nrow = 1) +
      labs(x = "Functional Domain", y = "Mean SE")

    se_plt_ov =
      ses_df %>%
      filter(var == "X") %>%
      ggplot(aes(x = l, y = value, color = boot_type, fill = boot_type)) +
      geom_smooth() +
      labs(x = "Functional Domain", y = "Mean SE")

    pdf(outfile)
    print(bias_plt / (se_plt + se_plt_ov + plot_layout(guides = "collect", axes = "collect")))
    dev.off()
  }


}

