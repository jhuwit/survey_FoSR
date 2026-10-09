## work on real world data result
library(tidyverse)
library(patchwork)
library(svyfosr)
## application with precision weighting
pa_df = read_rds(here::here("data", "mims_covariates.rds")) %>%
  mutate(SEQN = as.character(SEQN))


pa_df

## try binning and compare

pa_df_long = pa_df |>
  pivot_longer(cols = starts_with("min"),
               names_transform = ~as.numeric(sub(".*min\\_", "", .x)))

# summarize in 5 min bins
bins_05 =
  pa_df_long |>
  mutate(bin5 = ceiling(name / 5)) |>
  group_by(SEQN, bin5) |>
  summarize(value = mean(value), .groups = "drop") |>
  pivot_wider(names_from = bin5, values_from = value)

bins_10 =
  pa_df_long |>
  mutate(bin5 = ceiling(name / 10)) |>
  group_by(SEQN, bin5) |>
  summarize(value = mean(value), .groups = "drop") |>
  pivot_wider(names_from = bin5, values_from = value)

bins_30 =
  pa_df_long |>
  mutate(bin5 = ceiling(name / 30)) |>
  group_by(SEQN, bin5) |>
  summarize(value = mean(value), .groups = "drop") |>
  pivot_wider(names_from = bin5, values_from = value)


pa_df_covars =
  pa_df %>%
  mutate(race_cat =
           case_when(race_hispanic_origin == "Non-Hispanic White" ~ "NH White",
                     race_hispanic_origin == "Non-Hispanic Black" ~ "NH Black",
                     race_hispanic_origin %in% c("Other Hispanic", "Mexican American") ~ "Hispanic",
                     .default = "Other"),
         race_cat = factor(race_cat, levels = c("NH White", "NH Black", "Hispanic", "Other"))) |>
  filter(age_in_years_at_screening >= 18) %>% # just keep adults
  mutate(sex_bin = if_else(gender == "Female", 1, 0),
         weight = full_sample_2_year_mec_exam_weight / 2,
         age_cat = cut(age_in_years_at_screening, breaks=c(18, 30, 50, 65, 80), include.lowest = TRUE, right = FALSE)) |>
  select(-starts_with("min"))

pa_full = as.matrix(pa_df %>% select(starts_with("min")))


dat_full = cbind(pa_df_covars, Y = pa_full)
# all.equal(pa_df_covars$SEQN, bins_05$SEQN)
dat_05 = cbind(pa_df_covars, Y = as.matrix(bins_05 |> select(-SEQN)))
dat_10 = cbind(pa_df_covars, Y = as.matrix(bins_10 |> select(-SEQN)))
dat_30 = cbind(pa_df_covars, Y = as.matrix(bins_30 |> select(-SEQN)))

tictoc::tic()
fit1 = svyfosr::svyfui(Y ~ age_cat + sex_bin + race_cat,
                      data = dat_full,
                      weights = weight,
                      family = gaussian(),
                      boot_type = "BRR",
                      num_boots = 500,
                      parallel = TRUE,
                      seed = 2213,
                      nknots_min = 15)
tictoc::toc()
# 165.841 sec elapsed
tictoc::tic()
fit2 = svyfosr::svyfui(Y ~ age_cat + sex_bin + race_cat,
                       data = dat_05,
                       weights = weight,
                       family = gaussian(),
                       boot_type = "BRR",
                       num_boots = 500,
                       parallel = TRUE,
                       seed = 2213,
                       nknots_min = 15)
tictoc::toc()
# 162.473 sec elapsed
tictoc::tic()
fit3 = svyfosr::svyfui(Y ~ age_cat + sex_bin + race_cat,
                       data = dat_10,
                       weights = weight,
                       family = gaussian(),
                       boot_type = "BRR",
                       num_boots = 500,
                       parallel = TRUE,
                       seed = 2213,
                       nknots_min = 15)
tictoc::toc()
# 162.907 sec elapsed
tictoc::tic()
fit4 = svyfosr::svyfui(Y ~ age_cat + sex_bin + race_cat,
                       data = dat_30,
                       weights = weight,
                       family = gaussian(),
                       boot_type = "BRR",
                       num_boots = 500,
                       parallel = TRUE,
                       seed = 2213,
                       nknots_min = 15)
tictoc::toc()
# 157.376 sec elapsed

## application with precision weighting
step_df = read_rds(here::here("data", "steps_covariates.rds")) %>%
  mutate(SEQN = as.character(SEQN))



## try binning and compare

step_df_long = step_df |>
  pivot_longer(cols = starts_with("min"),
               names_transform = ~as.numeric(sub(".*min\\_", "", .x)))

# summarize in 5 min bins
bins_05 =
  step_df_long |>
  mutate(bin5 = ceiling(name / 5)) |>
  group_by(SEQN, bin5) |>
  summarize(value = mean(value), .groups = "drop") |>
  pivot_wider(names_from = bin5, values_from = value)

bins_10 =
  step_df_long |>
  mutate(bin5 = ceiling(name / 10)) |>
  group_by(SEQN, bin5) |>
  summarize(value = mean(value), .groups = "drop") |>
  pivot_wider(names_from = bin5, values_from = value)

bins_30 =
  step_df_long |>
  mutate(bin5 = ceiling(name / 30)) |>
  group_by(SEQN, bin5) |>
  summarize(value = mean(value), .groups = "drop") |>
  pivot_wider(names_from = bin5, values_from = value)


step_df_covars =
  step_df %>%
  rename(strata = masked_variance_pseudo_stratum,
         psu = masked_variance_pseudo_psu) |>
  mutate(race_cat =
           case_when(race_hispanic_origin == "Non-Hispanic White" ~ "NH White",
                     race_hispanic_origin == "Non-Hispanic Black" ~ "NH Black",
                     race_hispanic_origin %in% c("Other Hispanic", "Mexican American") ~ "Hispanic",
                     .default = "Other"),
         race_cat = factor(race_cat, levels = c("NH White", "NH Black", "Hispanic", "Other"))) |>
  filter(age_in_years_at_screening >= 18) %>% # just keep adults
  mutate(sex_bin = if_else(gender == "Female", 1, 0),
         weight = full_sample_2_year_mec_exam_weight / 2,
         age_cat = cut(age_in_years_at_screening, breaks=c(18, 30, 50, 65, 80), include.lowest = TRUE, right = FALSE)) |>
  select(-starts_with("min"))

step_full = as.matrix(step_df %>% select(starts_with("min"))) |> round(digits = 0)


dat_full = cbind(step_df_covars, Y = step_full)
dat_05 = cbind(step_df_covars, Y = as.matrix(bins_05 |> select(-SEQN)) |> round(digits = 0))
dat_10 = cbind(step_df_covars, Y = as.matrix(bins_10 |> select(-SEQN)) |> round(digits = 0))
dat_30 = cbind(step_df_covars, Y = as.matrix(bins_30 |> select(-SEQN)) |> round(digits = 0))

tictoc::tic()
fit1 = svyfosr::svyfui(Y ~ age_cat + sex_bin + race_cat,
                       data = dat_full,
                       weights = weight,
                       family = poisson(),
                       boot_type = "BRR",
                       num_boots = 500,
                       parallel = TRUE,
                       seed = 2213,
                       nknots_min = 15)
tictoc::toc()
# 3000.727 sec elapsed
tictoc::tic()
fit2 = svyfosr::svyfui(Y ~ age_cat + sex_bin + race_cat,
                       data = dat_05,
                       weights = weight,
                       family = poisson(),
                       boot_type = "BRR",
                       num_boots = 500,
                       parallel = TRUE,
                       seed = 2213,
                       nknots_min = 15)
tictoc::toc()
# 3059.453 sec elapsed
tictoc::tic()
fit3 = svyfosr::svyfui(Y ~ age_cat + sex_bin + race_cat,
                       data = dat_10,
                       weights = weight,
                       family = poisson(),
                       boot_type = "BRR",
                       num_boots = 500,
                       parallel = TRUE,
                       seed = 2213,
                       nknots_min = 15)
tictoc::toc()
# 3136.109 sec elapsed
tictoc::tic()
fit4 = svyfosr::svyfui(Y ~ age_cat + sex_bin + race_cat,
                       data = dat_30,
                       weights = weight,
                       family = poisson(),
                       boot_type = "BRR",
                       num_boots = 500,
                       parallel = TRUE,
                       seed = 2213,
                       nknots_min = 15)
tictoc::toc()
# 3157.032 sec elapsed



# compare fnl reg with taking the mean
pa_summary =
  pa_df_long |>
  group_by(SEQN) |>
  summarize(mean_pa = mean(value),
            .groups = "drop")

pa_df_mean =
  pa_df_covars |>
  left_join(pa_summary, by = "SEQN")

library(survey)
nhanes_des = svydesign(id = ~psu,  # Primary Sampling Units (PSU)
                          strata  = ~strata, # Stratification used in the survey
                          weights = ~weight,   # Survey weights
                          nest    = TRUE,      # Whether PSUs are nested within strata
                          data    = pa_df_mean)

mean_model = survey::svyglm(mean_pa ~ age_cat + sex_bin + race_cat, data = pa_df_mean, design = nhanes_des,
                            family = gaussian())

summary(mean_model)
# sex bin: female 1 male 0
