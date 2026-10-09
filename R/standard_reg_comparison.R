## work on real world data result
library(tidyverse)
library(patchwork)
library(svyfosr)
library(survey)
## application with precision weighting
pa_df = read_rds(here::here("data", "mims_covariates.rds")) %>%
  mutate(SEQN = as.character(SEQN))


pa_df_covars =
  pa_df %>%
  mutate(race_cat =
           case_when(race_hispanic_origin == "Non-Hispanic White" ~ "NH White",
                     race_hispanic_origin == "Non-Hispanic Black" ~ "NH Black",
                     race_hispanic_origin %in% c("Other Hispanic", "Mexican American") ~ "Hispanic",
                     .default = "Other"),
         race_cat = factor(race_cat, levels = c("NH White", "NH Black", "Hispanic", "Other"))) |>
  filter(age_in_years_at_screening >= 18) %>% # just keep adults
  mutate(gender = factor(gender, levels = c("Male", "Female")),
         weight = full_sample_2_year_mec_exam_weight / 2,
         age_cat = cut(age_in_years_at_screening, breaks=c(18, 30, 50, 65, 80), include.lowest = TRUE, right = FALSE)) |>
  select(-starts_with("min"))


pa_df_long = pa_df |>
  pivot_longer(cols = starts_with("min"),
               names_transform = ~as.numeric(sub(".*min\\_", "", .x)))

# compare fnl reg with taking the mean
pa_summary =
  pa_df_long |>
  group_by(SEQN) |>
  summarize(mean_pa = mean(value),
            .groups = "drop")

pa_df_mean =
  pa_df_covars |>
  left_join(pa_summary, by = "SEQN")


nhanes_des = svydesign(id = ~psu,  # Primary Sampling Units (PSU)
                       strata  = ~strata, # Stratification used in the survey
                       weights = ~weight,   # Survey weights
                       nest    = TRUE,      # Whether PSUs are nested within strata
                       data    = pa_df_mean)

mean_model = survey::svyglm(mean_pa ~ age_cat + gender + race_cat, data = pa_df_mean, design = nhanes_des,
                            family = gaussian())

summary(mean_model)

broom::tidy(mean_model) |>
  mutate(p.value = format.pval(p.value, digits = 2),
         across(c(estimate, std.error, statistic), ~round(.x, 3))) |>
  kableExtra::kbl(format = "latex", booktabs = TRUE)
