### process nhanes data to get minute-level steps per subject

library(tidyverse)
force = FALSE

if (!file.exists(here::here("data", "processed_steps_multilevel.rds")) || force) {
  steps = read_csv(here::here("data", "nhanes_1440_scsslsteps.csv.xz"))
  wear = read_csv(here::here("data", "nhanes_1440_PAXPREDM.csv.xz"))
  flags = read_csv(here::here("data", "nhanes_1440_PAXFLGSM.csv.xz"))
  mims = read_csv(here::here("data", "nhanes_1440_PAXMTSM.csv.xz"))

  all.equal(wear %>% select(SEQN, PAXDAYM), flags %>% select(SEQN, PAXDAYM))
  all.equal(mims %>% select(SEQN, PAXDAYM), flags %>% select(SEQN, PAXDAYM))
  all.equal(steps %>% select(SEQN, PAXDAYM), flags %>% select(SEQN, PAXDAYM))


  mims_mat = mims %>% select(starts_with("min")) %>% as.matrix()
  wear_mat = wear %>% select(starts_with("min")) %>% as.matrix()
  flag_mat = flags %>% select(starts_with("min")) %>% as.matrix()
  ## replace entires in wear mat with NA if flag is TRUE
  wear_mat_filt = ifelse(flag_mat, NA, wear_mat)
  ## replace entries in mims mat with NA if flag is TRUE
  mims_mat_filt = ifelse(flag_mat, NA, mims_mat)

  # create logical matrix where wear is 1, 2, or 4
  wear_mat_lgl = matrix(wear_mat_filt %in% c(1, 2, 4),
                        nrow = nrow(wear_mat_filt), ncol = ncol(wear_mat_filt))
  # create logical matrix where true if MIMS > 0
  mims_mat_lgl = matrix(mims_mat_filt > 0,
                        nrow = nrow(mims_mat_filt), ncol = ncol(mims_mat_filt))

  wake_mat_lgl = matrix(wear_mat_filt == 1,
                        nrow = nrow(wear_mat_filt), ncol = ncol(wear_mat_filt))

  # use rowsums to quickly get total number of minute per day with wear
  wear_summary_vec = rowSums(wear_mat_lgl, na.rm = TRUE)
  # use rowsums to get total number of minutes with nonzero MIMS
  mims_summary_vec = rowSums(mims_mat_lgl, na.rm = TRUE)
  # same for wake
  wake_summary_vec = rowSums(wake_mat_lgl, na.rm = TRUE)

  # get wear summary
  wear_summary =
    wear %>%
    select(SEQN, PAXDAYM) %>%
    bind_cols(non_flag_wear = wear_summary_vec) %>%
    bind_cols(non_zero_MIMS = mims_summary_vec) %>%
    bind_cols(wake_wear = wake_summary_vec) %>%
    mutate(include = non_flag_wear >= 1368 & wake_wear >= 420 & non_zero_MIMS >= 420)


  include_days =
    wear_summary %>%
    group_by(SEQN) %>%
    mutate(days_per_sub = sum(include)) %>%
    ungroup() %>%
    mutate(include_final = days_per_sub >= 3 & include) %>%
    filter(include_final) %>%
    select(SEQN, PAXDAYM)

  steps_small =
    steps %>%
    right_join(include_days, by = c("SEQN", "PAXDAYM"))

  wear_small =
    wear %>%
    right_join(include_days, by = c("SEQN", "PAXDAYM"))

  flags_small =
    flags %>%
    right_join(include_days, by = c("SEQN", "PAXDAYM"))

  # get rid of steps minutes that are nonwear OR flagged
  wear_mat = wear_small %>% select(starts_with("min")) %>% as.matrix()
  flag_mat = flags_small %>% select(starts_with("min")) %>% as.matrix()
  steps_mat = steps_small %>% select(starts_with("min")) %>% as.matrix()

  wear_mat_lgl = matrix(wear_mat == 3,
                        nrow = nrow(wear_mat), ncol = ncol(wear_mat))
  ## replace entires in wear mat with NA if flag is TRUE
  steps_mat_filt = ifelse(wear_mat_lgl, NA, steps_mat)
  ## replace entries in mims mat with NA if flag is TRUE
  steps_mat_filt2 = ifelse(flag_mat, NA, steps_mat_filt)

  steps_final =
    steps_small %>%
    select(-starts_with("min")) %>%
    bind_cols(steps_mat_filt2)


  write_rds(steps_final, here::here("data", "processed_steps_multilevel.rds"))
} else {
  steps_final = read_rds(here::here("data", "processed_steps_multilevel.rds"))
}



steps_persub =
  steps_final %>%
  group_by(SEQN) %>%
  summarize(across(starts_with("min"), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")

# write_rds(steps_persub, here::here("data", "processed_steps_sl.rds"))
# steps_persub = read_rds(here::here("data", "processed_steps_sl.rds"))
covars = read_rds(here::here("data", "covariates_accel_mortality_df.rds"))


covars_filt = covars %>%
  filter(SEQN %in% steps_persub$SEQN) %>%
  select(SEQN, data_release_cycle, gender, age_in_years_at_screening, full_sample_2_year_mec_exam_weight,
         race_hispanic_origin, masked_variance_pseudo_psu, masked_variance_pseudo_stratum,
         annual_household_income, cat_alcohol, cat_bmi, cat_smoke, bin_mobilityproblem,
         cat_education, mortstat) %>%
  left_join(steps_persub %>% mutate(SEQN = as.character(SEQN)), by = "SEQN")

write_rds(covars_filt, here::here("data", "steps_covariates.rds"))

#
# covars_filt %>%
#   slice(1:100) %>%
#   select(SEQN, starts_with("min")) %>%
#   pivot_longer(cols = starts_with("min")) %>%
#   mutate(s = as.numeric(sub(".*min\\_", "", name))) %>%
#   ggplot(aes(x = s, y = value)) +
#   geom_smooth()
