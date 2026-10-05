## work on real world data result
library(tidyverse)
library(patchwork)
library(svyfosr)
## application with precision weighting
pa_df = read_rds(here::here("data", "mims_covariates.rds")) %>%
  mutate(SEQN = as.character(SEQN))


n_days = read_rds(here::here("data", "valid_days_per_subj.rds"))

pa_df =
  pa_df |>
  left_join(n_days, by = "SEQN")

pa_df =
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
         age_cat = cut(age_in_years_at_screening, breaks=c(18, 30, 50, 65, 80), include.lowest = TRUE, right = FALSE),
         pweight = weight * (n_days / 7))

pa_df |>
  select(SEQN, weight, pweight) |>
  pivot_longer(cols = -SEQN) |>
  ggplot(aes(x = value, fill = name, color = name)) +
  geom_density(alpha = .5)

Y = as.matrix(pa_df %>% select(starts_with("min")))

dat_full =
  pa_df %>%
  select(-starts_with("min"))

dat_full = cbind(dat_full, Y)

fit = svyfosr::svyfui(Y ~ age_cat + sex_bin + race_cat,
                      data = dat_full,
                      weights = weight,
                      family = gaussian(),
                      boot_type = "BRR",
                      num_boots = 500,
                      parallel = TRUE,
                      seed = 2213,
                      nknots_min = 15)

fit_pwt = svyfosr::svyfui(Y ~ age_cat + sex_bin + race_cat,
                      data = dat_full,
                      weights = pweight,
                      family = gaussian(),
                      boot_type = "BRR",
                      num_boots = 500,
                      parallel = TRUE,
                      seed = 2213,
                      nknots_min = 15)

fit_df = fit$tidy_df
fitpwt_df = fit_pwt$tidy_df

plot_df = fit_df |> mutate(type = "Survey weighted") |>
  bind_rows(fitpwt_df |> mutate(type = "Survey and precision weighted")) |>
  mutate(var_name = factor(var_name,
                       levels = c("(Intercept)", "age_cat[30,50)",   "age_cat[50,65)",   "age_cat[65,80]",  "sex_bin", "race_catNH Black",
                                  "race_catHispanic", "race_catOther"),
                       labels = c("Intercept", "Age 30-49", "Age 50-64", "Age 65+", "Sex: Female", "Non-Hispanic Black", "Hispanic", "Other race")))

col1 = "#009E73FF"
col2 = "#E69F00FF"
col2 = "#D55E00"
col1 = "#0072B2"
plot_df |>
  ggplot(aes(x = l, y = beta_hat, color = type)) +
  geom_line() +
  facet_wrap(.~var_name, scales = "free_y") +
  scale_x_continuous(breaks= c(0, seq(240, 1440, 6 * 60)), labels = c("", "06:00", "12:00",  "18:00", "24:00")) +
  labs(x = "Hour of Day", y = "Coefficient Estimate") +
  theme(strip.text = element_text(size = 14)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 14),
                   title = element_text(size = 18)) +
  theme_sub_axis_x(text = element_text(size = 14),
                   # text = element_text(size = 11, angle = 10),
                   title = element_text(size = 18)) +
  theme_sub_axis_y(text = element_text(size = 14),
                   title = element_text(size = 18)) +
  scale_color_manual(values = c(col1, col2), name = "Weighting type")

p =
  plot_df |>
  ggplot(aes(x = l, y = beta_hat)) +
  geom_ribbon(aes(ymin = lower_joint, ymax = upper_joint, fill = type), color = NA, alpha = .4) +
  geom_ribbon(aes(ymin = lower_pw, ymax = upper_pw, fill = type), color = NA, alpha = .6) +
  geom_line(aes(color = type), linewidth = 1.1) +
  facet_wrap(.~var_name, scales = "free_y") +
  scale_x_continuous(breaks= c(0, seq(240, 1440, 6 * 60)), labels = c("", "06:00", "12:00",  "18:00", "24:00")) +
  labs(x = "Hour of Day", y = "Coefficient Estimate") +
  theme(strip.text = element_text(size = 14)) +
  theme_sub_legend(position = c(0.8, 0.15),
                   text = element_text(size = 14),
                   title = element_text(size = 18)) +
  theme_sub_axis_x(text = element_text(size = 14),
                   # text = element_text(size = 11, angle = 10),
                   title = element_text(size = 18)) +
  theme_sub_axis_y(text = element_text(size = 14),
                   title = element_text(size = 18)) +
  scale_color_manual(values = c(col1, col2), name = "Weighting type") +
  scale_fill_manual(values = c(col1, col2), name = "Weighting type")



png(here::here("manuscript", "figures", "prec_wtd.png"),
    width = 12, height = 8, res = 400, units = "in")
p
dev.off()

