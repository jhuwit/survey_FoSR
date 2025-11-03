# get correlation between strata
library(tidyverse)

if (!file.exists(here::here("data", "steps_covariates.rds"))) {
  source(here::here("R", "process_steps_data.R"))
}
steps_df = read_rds(here::here("data", "processed_steps_multilevel.rds"))
covars = read_rds(here::here("data", "covariates_accel_mortality_df.rds"))

covars = covars %>%
  select(SEQN, strata = masked_variance_pseudo_stratum, psu = masked_variance_pseudo_psu) %>%
  mutate(SEQN = as.numeric(SEQN))

steps_data_persub =
  steps_df %>%
  group_by(SEQN) %>%
  summarize(across(starts_with("min"), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")

steps_data_persub =
  steps_data_persub %>%
  left_join(covars, by = "SEQN") %>%
  mutate(strata_psu = paste0(strata, "_", psu)) %>%
  arrange(strata_psu) %>%
  group_by(strata_psu) %>%
  mutate(id = row_number()) %>%
  ungroup()

Y_mat = steps_data_persub %>%
  select(contains("min")) %>%
  as.matrix()

ids = steps_data_persub$strata_psu # treat the strata/psu as level 1 (analagous to subject in typical application)
visits = steps_data_persub$id # treat the person within as level 2 (analagous to visit in typical application)

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

PVE_between
PVE_within

# > PVE_between
# [1] 0.03088735
# > PVE_within
# [1] 0.9691126


## do the same thing for MIMS

mims_df = read_csv(here::here("data", "nhanes_1440_PAXMTSM.csv.xz"))


mims_df_filt =
  mims_df %>%
  right_join(steps_df %>% select(PAXDAYM, SEQN), by = c("PAXDAYM", "SEQN")) # to get rid of "bad" days

mims_data_persub =
  mims_df_filt %>%
  group_by(SEQN) %>%
  summarize(across(starts_with("min"), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")

mims_data_persub =
  mims_data_persub %>%
  left_join(covars) %>%
  mutate(strata_psu = paste0(strata, "_", psu)) %>%
  arrange(strata_psu) %>%
  group_by(strata_psu) %>%
  mutate(id = row_number()) %>%
  ungroup()



Y_mat = mims_data_persub %>%
  select(contains("min")) %>%
  as.matrix()

ids = mims_data_persub$strata_psu # treat the strata/psu as level 1 (analagous to subject in typical application)
visits = mims_data_persub$id # treat the person within as level 2 (analagous to visit in typical application)

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

PVE_between
PVE_within
#
# > PVE_between
# [1] 0.0231469
# > PVE_within
# [1] 0.9768531
