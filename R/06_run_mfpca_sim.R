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
force = FALSE
force_iter = FALSE
source(here::here("R_cp", "01_sim_functions.R"))
source(here::here("R_cp", "00_data_gen_function_ff.R"))
source(here::here("R", "utils.R"))
source(here::here("R_cp", "create_survey_settings_final.R"))


ifold = get_fold()
options(survey.lonely.psu = "adjust")
nsim = 50


out_file = here::here("results", "mfpca_sim", paste0("fold_", sprintf("%03d", ifold), ".rds"))
if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file))

temp = settings[ifold, ] # set settings for this simulation

if (!file.exists(out_file) & temp$len != 1440) {

  partial_dir = here::here("results", "mfpca_sim", paste0("fold_", sprintf("%03d", ifold), "_partials"))
  if (!dir.exists(partial_dir)) dir.create(partial_dir)


  print(ifold)
  temp = settings[ifold, ] # set settings for this simulation

  lst = generate_superpopulation(
    family = temp$family,
    I = 10e6,
    L = temp$len,
    snr_b = temp$snr_b,
    snr_eps = temp$snr_eps,
    strata_sigma = temp$strata_sigma,
    strata_scale = temp$strata_scale,
    seed = 111
  )


  for(iter in 1:nsim){
    completed_iters = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$") %>%
      str_extract("\\d+") %>%
      as.integer()

    if (iter %in% completed_iters && !force_iter) {
      message("Skipping completed iteration: ", iter)
      next
    }

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
      refund::mfpca.face(Y = Y_mat, id = ids, # strata psu id
                         visit = visits) # "visit" within strata/psu
    eigenvalues_between <- mfpca_res$evalues[[1]]  # Between-cluster variance
    eigenvalues_within <- mfpca_res$evalues[[2]]

    cumulative_var_between <- cumsum(eigenvalues_between) / sum(eigenvalues_between)
    cumulative_var_within <- cumsum(eigenvalues_within) / sum(eigenvalues_within)

    npc_between <- which(cumulative_var_between >= 0.99)[1]
    npc_within <- which(cumulative_var_within >= 0.99)[1]

    total_var_bw = sum(eigenvalues_between[1:npc_between])
    total_var_win = sum(eigenvalues_within[1:npc_within])


    PVE_between <- total_var_bw / (total_var_bw + total_var_win)
    PVE_within <- total_var_win / (total_var_bw + total_var_win)
    res = tibble(
      fold = ifold,
      pve_b = PVE_between,
      pve_w = PVE_within,
      tv_bw = total_var_bw,
      tv_win = total_var_win
    )
    write_rds(res, file = file.path(partial_dir, sprintf("iter_%03d.rds", iter)))
    print(iter)
    rm(data, Y_mat, res)
  }


  partial_files = list.files(partial_dir, pattern = "^iter_\\d+\\.rds$", full.names = TRUE)

  sim_res = map(partial_files, \(x) read_rds(x) %>% mutate(id = sub(".*iter_(.+)\\.rds.*", "\\1", basename(x)))) %>%
    keep(., is.data.frame) %>%
    list_rbind()

  write_rds(sim_res, file = out_file)
  unlink(partial_dir, recursive = TRUE)

}



