
library(dplyr)
library(survey)

df = read_rds(here::here("sim_data", "covariates_accel_mortality_df.rds")) %>%
  mutate(SEQN = as.numeric(SEQN)) %>%
  rename(psu = masked_variance_pseudo_psu,
         strata = masked_variance_pseudo_stratum)

## functional data
fnl_data = read_csv(here::here("sim_data", "nhanes_1440_PAXMTSM.csv.xz"))
inclusion = read_csv(here::here("sim_data", "inclusion_summary.csv.gz"))

days_keep =
  inclusion %>%
  filter(include) %>%
  select(SEQN, PAXDAYM)

fnl_data =
  fnl_data %>%
  right_join(days_keep, by = c("SEQN", "PAXDAYM"))

pa_df =
  fnl_data %>%
  group_by(SEQN) %>%
  summarize(across(contains("min"), ~mean(.x, na.rm = TRUE)),
            .groups = "drop")

df_small =
  df %>%
  select(SEQN, data_release_cycle, psu, strata,
         gender, age_in_years_at_screening, full_sample_2_year_interview_weight,
         full_sample_2_year_mec_exam_weight,
         cat_bmi)

pa_df =
  pa_df %>%
  left_join(df_small)  %>%
  select(SEQN, data_release_cycle, psu, strata,
         gender, age_in_years_at_screening, full_sample_2_year_interview_weight,
         full_sample_2_year_mec_exam_weight,
         cat_bmi, starts_with("min"))


write_rds(pa_df, here::here("sim_data", "pa_df_persub.rds"), compress = "gz")

