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
library(svrep)
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "00_data_gen_function_twostage.R"))
source(here::here("R", "utils.R"))


options(survey.lonely.psu = "adjust")


settings1 = expand_grid(snr_b = c(0.25, 0.5, 1, 2.5, 5),
                        snr_eps = c(0.25, 0.5, 1, 2.5, 5),
                        psu_re = 0.5,
                        In = c(100, 250, 500),
                        len = c(50, 100),
                        scen = c(1, 2),
                        family = "gaussian")

settings2 = expand_grid(snr_b = c(0.25, 0.5, 1, 2.5, 5),
                        psu_re = 0.5,
                        In = c(100, 250, 500),
                        len = c(50, 100),
                        scen = c(1, 2),
                        family = c("poisson", "binomial"))



settings = bind_rows(settings1, settings2) %>%
  mutate(fold = row_number())
res_list = list()
for(ifold in 1:nrow(settings)){
  print(ifold)
  x = try({
    temp = settings[ifold,] # set settings for this simulation

  strata_sigma = 0.05
  psu_sigma = sqrt(strata_sigma ^ 2 * temp$psu_re)

  lst = generate_superpopulation(scenario = temp$scen,
                                 family = temp$family,
                                 I = 10e6,
                                 L = temp$len,
                                 psu_sigma = psu_sigma,
                                 snr_b = temp$snr_b,
                                 snr_eps = temp$snr_eps)



  # sim_res <- list() ## store simulation results
  # for (iter in 1:nsim) {
  #   set.seed(iter)
  #   x = try({
  data <- sample_from_population_wor(
    X_des = lst$X_des,
    Y_obs = lst$Y_obs,
    L = temp$len,
    ## dimension of the functional domain
    I_n = temp$In,
    num_strata = 30,
    stratum_assignments = lst$stratum_assignments,
    psu_assignments = lst$psu_assignments,
    dirichlet_probs = lst$dirichlet_probs,
    seed = 1
  ) %>%
    mutate(strata_psu = paste0(strata, "_", psu))

  data = data %>%
    group_by(strata_psu) %>%
    mutate(id = row_number()) %>%
    ungroup()
  ids = data$strata_psu
  visits = data$id

  Y_mat =
    data %>%
    select(starts_with("Y")) %>%
    as.matrix()

  mfpca_res =
    refund::mfpca.face(Y = Y_mat,
                       id = ids, # strata psu id
                       visit = visits) # "visit" within strata/psu
  eigenvalues_between <- mfpca_res$evalues[[1]]  # Between-cluster variance
  eigenvalues_within <- mfpca_res$evalues[[2]]

  cumulative_var_between <- cumsum(eigenvalues_between) / sum(eigenvalues_between)
  cumulative_var_within <- cumsum(eigenvalues_within) / sum(eigenvalues_within)

  npc_between <- which(cumulative_var_between >= 0.99)[1]
  npc_within <- which(cumulative_var_within >= 0.99)[1]

  total_var_bw = sum(eigenvalues_between[1:npc_between]); total_var_bw
  total_var_win = sum(eigenvalues_within[1:npc_within]); total_var_win


  PVE_between <- total_var_bw / (total_var_bw + total_var_win)
  PVE_within <- total_var_win / (total_var_bw + total_var_win)
  res = tibble(fold = ifold,
         pve_b = PVE_between,
         pve_w = PVE_within,
         tv_bw = total_var_bw,
         tv_win = total_var_win)
  rm(lst, data, Y_mat, mfpca_res, ids, visits,
     eigenvalues_between, eigenvalues_within, cumulative_var_between, cumulative_var_within,
     npc_between, npc_within, total_var_bw, total_var_win, PVE_between, PVE_within)
  res
  })
  res_list[[ifold]] = x
  rm(x)
}

write_rds(res_list, here::here("results", "mfpca_simulation.rds"))

# sim = read_rds(here::here("results", "mfpca_simulation.rds")) %>%
#   keep(is_tibble) %>%
#   bind_rows()
#
# sim = sim %>%
#   left_join(settings %>% mutate(fold = row_number()), by = "fold")
#
# sim %>%
#   filter(len == 50, family == "gaussian", psu_re == 0.5) %>%
#   ggplot(aes(x = factor(snr_b), y = pve_b, color = factor(snr_eps))) +
#   geom_point() +
#   facet_grid(In~scen)
#
# sim %>%
#   filter(len == 50, family != "gaussian", psu_re == 0.5) %>%
#   ggplot(aes(x = factor(snr_b), y = pve_b, color = factor(family))) +
#   geom_point() +
#   facet_grid(In~scen)
