## investigate informativeness
# library(future)
# library(furrr)
library(dplyr)
library(purrr)
library(readr)
# library(tidyverse)
library(here)
# library(devtools)
# library(fastFMM)
library(survey)
library(progress)
library(lme4)
# library(paletteer)
library(mgcv)
# library(ggplot2)
# library(gridExtra)
# library(tidyfun)
library(mvtnorm)
library(refund)
library(svrep)
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "00_data_gen_function_ff.R"))
source(here::here("R", "00_data_gen_function_fast.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "create_survey_settings.R"))


lst = generate_superpopulation(
  family = "gaussian",
  I = 10e5,
  L = 100,
  snr_b = 0.5,
  snr_eps = 1,
  strata_sigma = 0.125,
  strata_scale = 0.05,
  seed = 111
)

get_ratio = function(df, p1 = 0.25, p2 = 0.75){
  x1 =
    df |>
    filter(ptile <= p1) |>
    summarize(p = mean(p_overall)) |>
    pull(p)

  x2 =
    df |>
    filter(ptile >= p2) |>
    summarize(p = mean(p_overall)) |>
    pull(p)

  return(x2/x1)
}

get_stats = function(inf, iter = 1) {
  fold = settings |>
    filter(family == "gaussian", strata_scale > 0, strata_sigma > 0, snr_b == 0.5, snr_eps == 1, len == 100, In == 100, inf_level == inf) |>
    pull(fold)

  temp = settings[fold,]

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
    inf_level = temp$inf_level,
    compression = 2,
    family = temp$family
  )


  Y_mat =
    data |>
    select(starts_with("Y")) |>
    as.matrix()

  ybar = rowMeans(Y_mat)

  data_means =
    data |>
    select(-starts_with("Y")) |>
    mutate(ybar = ybar,
           p_overall = p_stage1 * p_stage2)


  c_pear = cor(data_means$p_overall, data_means$ybar, method = "pearson")
  c_spear = cor(data_means$p_overall, data_means$ybar, method = "spearman")

  data_means =
    data_means |>
    mutate(ptile = percent_rank(ybar))

  ratio_75 = get_ratio(data_means, 0.25, 0.75)
  ratio_9 = get_ratio(data_means, 0.1, 0.9)

  return(tibble(c_pear = c_pear,
                c_spear = c_spear,
                ratio_75 = ratio_75,
                ratio_9 = ratio_9,
                iter = iter,
                inf_level = inf))

}


comb = tidyr::expand_grid(iter = seq_len(10),
                          inf = c(0, 10, 15))


res = map2_dfr(.x = comb$inf, .y = comb$iter, .f = get_stats)

write_rds(res, here::here("results", "inf_level_exploration.rds"))
#
#
# res |>
#   pivot_longer(cols = contains("c_")) |>
#   ggplot(aes(x = name, y = value, fill = factor(inf_level))) +
#   geom_boxplot() +
#   scale_x_discrete(labels = c("Pearson", "Spearman"), name = "Correlation type") +
#   scale_y_continuous(name = "Correlation between functional outcome and selection probability") +
#   scale_fill_discrete(name = "Informativeness", labels = c("None", "Medium", "High"))
#
# res |>
#   pivot_longer(cols = contains("ratio")) |>
#   ggplot(aes(x = name, y = value, fill = factor(inf_level))) +
#   geom_boxplot() +
#   scale_x_discrete(labels = c("Ratio 75th:25th", "Ratio 90th:10th")) +
#   scale_fill_discrete(name = "Informativeness", labels = c("None", "Medium", "High"))



