## work on real world data result
library(tidyverse)
library(patchwork)
#### ---- run models if not run ---- ###
if (!file.exists(here::here("results", "mims_application_1440.rds"))) {
  pa_df = read_rds(here::here("data", "mims_covariates.rds")) %>%
    mutate(SEQN = as.character(SEQN))

  ## sample size

  pa_df %>%
    count(gender)


  source(here::here("R", "01_sim_functions.R"))
  n_boot = 500

  L = 1440
  pa_df =
    pa_df %>%
    mutate(race_cat =
             case_when(race_hispanic_origin == "Non-Hispanic White" ~ "NH White",
                       race_hispanic_origin == "Non-Hispanic Black" ~ "NH Black",
                       race_hispanic_origin %in% c("Other Hispanic", "Mexican American") ~ "Hispanic",
                       .default = "Other"),
            race_cat = factor(race_cat, levels = c("NH White", "NH Black", "Hispanic", "Other")))



  fit_types = c('weighted', 'unweighted')
  boot_types = c('BRR', 'weighted', 'unweighted')


  data =
    pa_df %>%
    filter(age_in_years_at_screening >= 18) %>% # just keep adults
    mutate(sex_bin = if_else(gender == "Female", 1, 0),
           weight = full_sample_2_year_mec_exam_weight / 2,
           age_cat = cut(age_in_years_at_screening, breaks=c(18, 30, 50, 65, 80), include.lowest = TRUE, right = FALSE))
  # rename_with(.cols = starts_with("min"), ~str_replace(., "min_", "Y"))

  data %>% count(gender)
  data %>% count(age_cat)
  data %>% count(race_cat)
  data %>%
    group_by(gender) %>%
    summarize(mean_weight = mean(weight), sd_weight = sd(weight))

  data %>%
    group_by(age_cat) %>%
    summarize(mean_weight = mean(weight), sd_weight = sd(weight))

  data %>%
    group_by(race_cat) %>%
    summarize(mean_weight = mean(weight), sd_weight = sd(weight))

  Y = as.matrix(data %>% select(starts_with("min")))
  colnames(Y) = paste0("Y", 1:L)

  # betaTilde = map(.x = fit_types, .f = get_betatilde, data = data, family = "poisson")

  dat_full =
    data %>%
    select(-starts_with("min"))
  dat_full = cbind(dat_full, Y)


  betaTilde = map(.x = fit_types, .f = get_betatilde, data = dat_full, family = "gaussian",
                  model_formula = as.formula(paste0('Y~', 'age_cat + sex_bin + race_cat')))
  betaHat = map(.x = betaTilde, .f = get_betahat, L = L, nknots_min = 15)




  data = data %>%
    mutate(ID = SEQN)

  # Set the maximum global size to 1.13 GiB
  options(future.globals.maxSize = 3* 1024^3)  # Equivalent to 1.13 GiB

  tictoc::tic()
  res = map(
    .x = boot_types,
    .f = run_boots,
    data = dat_full %>% mutate(ID = SEQN),
    weights = dat_full$weight,
    family = "gaussian",
    num_boots = 500,
    parallel = TRUE,
    n_cores = 6,
    seed = 2025,
    model_formula = as.formula(paste0('Y~', 'age_cat + sex_bin + race_cat'))
  )

  fit_list = list(betaHat[[1]], betaHat[[1]], betaHat[[2]])

  cis = map2(
    .x = res,
    .y = fit_list,
    .f = get_cis_par,
    L = L,
    nknots_min = 15,
    nknots_min_fpca = 15,
    parallel = TRUE,
    n_cores = 7,
    smooth_for_ci = TRUE
  )
  elapsed_time = tictoc::toc(quiet = TRUE)

  write_rds(cis, here::here("results", "mims_application_1440.rds"))
} else  cis = read_rds(here::here("results", "mims_application_1440.rds"))

### ---------- make plots ------------ ###

get_plot_df = function(sim_res) {
  num_var = nrow(sim_res$betaHat)
  map_dfr(
    .x = 1:num_var,
    .f = function(r) {
      data.frame(
        s = 1:length(sim_res$betaHat[r, ]),
        beta = sim_res$betaHat[r, ],
        lower = sim_res$betaHat[r, ] - 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
        upper = sim_res$betaHat[r, ] + 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
        lower.joint = sim_res$betaHat[r, ] - sim_res$qn[r] *
          sqrt(diag(sim_res$betaHat.var[, , r])),
        upper.joint = sim_res$betaHat[r, ] + sim_res$qn[r] * sqrt(diag(sim_res$betaHat.var[, , r]))
      ) %>%
        mutate(name = sim_res$betaHat %>% rownames() %>% .[r])
    }
  )
}
plt_df_boot = map(cis, .f = get_plot_df)


col1 = "#009E73FF"
col2 = "#E69F00FF"



extra_rows =
  tibble(name = c("(Intercept)", "age_cat[30,50)", "age_cat[50,65)", "age_cat[65,80]", "sex_bin", "race_catNH Black", "race_catHispanic",
                  "race_catOther"),
         min =  c(0, rep(-6, 7)),
         max = c(0, rep(4, 7))) %>%
  mutate(name = factor(name,
                       levels = c("(Intercept)", "age_cat[30,50)",   "age_cat[50,65)",   "age_cat[65,80]",  "sex_bin", "race_catNH Black",
                                  "race_catHispanic", "race_catOther"),
                       labels = c("Intercept", "Age 30-49", "Age 50-64", "Age 65+", "Sex: Female", "Non-Hispanic Black", "Hispanic", "Other race")))

extra_rows =
  extra_rows %>% mutate(model = "BRR") %>%
  bind_rows(extra_rows %>% mutate(model = "Unweighted")) %>%
  bind_rows(extra_rows %>% mutate(model = "Weighted")) %>%
  mutate(modelf = factor(model, levels = c("Unweighted", "Weighted", "BRR")))

df =
  plt_df_boot %>%
  bind_rows(.id = "model") %>%
  as_tibble() %>%
  filter(model %in% c("1", "2", "3")) %>%
  mutate(modelf = recode_factor(model,
                "1" = "BRR",
                "2" = "Weighted",
                "3" = "Unweighted")) %>%
  mutate(name = factor(name,
                       levels = c("(Intercept)", "age_cat[30,50)",   "age_cat[50,65)",   "age_cat[65,80]",  "sex_bin", "race_catNH Black",
                                  "race_catHispanic", "race_catOther"),
                       labels = c("Intercept", "Age 30-49", "Age 50-64", "Age 65+", "Sex: Female", "Non-Hispanic Black", "Hispanic", "Other race")))
### calculate significance

## all three methods
t0 = df %>%
  mutate(sig = if_else((lower < 0 & upper < 0) | (lower > 0 & upper > 0), 1, 0)) %>%
  select(modelf, s, beta, sig, name) %>%
  select(-beta) %>%
  pivot_wider(names_from = modelf, values_from = sig) %>%
  # filter(name != "Intercept") %>%
  mutate(sum = BRR + Weighted + Unweighted) %>%
  rowwise() %>%
  mutate(agree = (sum == 3 | sum == 0)) %>%
  ungroup() %>%
  group_by(name) %>%
  summarize(pct_agree = sum(agree) / 1440 * 100)

## BRR and weighted
t1 = df %>%
  mutate(sig = if_else((lower < 0 & upper < 0) | (lower > 0 & upper > 0), 1, 0)) %>%
  select(modelf, s, beta, sig, name) %>%
  select(-beta) %>%
  pivot_wider(names_from = modelf, values_from = sig) %>%
  # filter(name != "Intercept") %>%
  mutate(sum = BRR + Weighted ) %>%
  rowwise() %>%
  mutate(agree = (sum == 2 | sum == 0)) %>%
  ungroup() %>%
  group_by(name) %>%
  summarize(pct_agree = sum(agree) / 1440 * 100)

# weighted and unweighted
t2 = df %>%
  mutate(sig = if_else((lower < 0 & upper < 0) | (lower > 0 & upper > 0), 1, 0)) %>%
  select(modelf, s, beta, sig, name) %>%
  select(-beta) %>%
  pivot_wider(names_from = modelf, values_from = sig) %>%
  # filter(name != "Intercept") %>%
  mutate(sum = Weighted + Unweighted ) %>%
  rowwise() %>%
  mutate(agree = (sum == 2 | sum == 0)) %>%
  ungroup() %>%
  group_by(name) %>%
  summarize(pct_agree = sum(agree) / 1440 * 100)


# BRR and unweighted
t3 = df %>%
  mutate(sig = if_else((lower < 0 & upper < 0) | (lower > 0 & upper > 0), 1, 0)) %>%
  select(modelf, s, beta, sig, name) %>%
  select(-beta) %>%
  pivot_wider(names_from = modelf, values_from = sig) %>%
  # filter(name != "Intercept") %>%
  mutate(sum = BRR + Unweighted ) %>%
  rowwise() %>%
  mutate(agree = (sum == 2 | sum == 0)) %>%
  ungroup() %>%
  group_by(name) %>%
  summarize(pct_agree = sum(agree) / 1440 * 100)

## agreement table
t0 %>% rename(all = pct_agree) %>%
  left_join(t1 %>% rename(bw = pct_agree)) %>%
  left_join(t2 %>% rename(wu = pct_agree)) %>%
  left_join(t3 %>% rename(bu = pct_agree)) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 1, nsmall = 0))) %>%
  kableExtra::kable(format = "latex", booktabs = TRUE)

# check CI widths - means
df %>%
  mutate(width = upper - lower) %>%
  group_by(modelf, name) %>%
  summarize(mw = mean(width)) %>%
  pivot_wider(names_from = modelf, values_from = mw) %>%
  mutate(pct = (BRR - Unweighted) / Unweighted * 100)
df %>%
  mutate(width = upper.joint - lower.joint) %>%
  group_by(modelf, name) %>%
  summarize(mw = mean(width)) %>%
  pivot_wider(names_from = modelf, values_from = mw) %>%
  mutate(pct = (BRR - Unweighted) / Unweighted * 100)

## create DF that is difference in CI widths
diff_df =
  df %>%
  mutate(width = upper - lower,
         width_joint = upper.joint - lower.joint) %>%
  select(modelf, s, name, width, width_joint, beta) %>%
  group_by(s, name) %>%
  mutate(beta_avg = beta[modelf == "Weighted"]) %>%
  ungroup() %>%
  select(-beta) %>%
  pivot_longer(cols = starts_with("width"), names_to = "type", values_to = "w") %>%
  pivot_wider(names_from = modelf, values_from = w) %>%
  mutate(brr_minus_unweighted = (BRR - Unweighted) / Unweighted * 100,
         brr_minus_weighted = (BRR - Weighted) / Weighted * 100,
         weighted_minus_unweighted = (Weighted - Unweighted) / Unweighted * 100) %>%
  pivot_longer(cols = c(brr_minus_unweighted, brr_minus_weighted, weighted_minus_unweighted),
               names_to = "diff_type",
               values_to = "diff_pct")

diff_df =
  diff_df %>%
  filter(type == "width")

diff_df %>%
  group_by(diff_type) %>%
  summarize(across(diff_pct, list(min = min, max = max, mean = mean, median = median), .names = "{col}_{fn}"))
range_pct  = range(diff_df$diff_pct, na.rm = TRUE)
range_beta  = range(diff_df$beta_avg, na.rm = TRUE)

b = diff(range_beta) / diff(range_pct)
a = range_beta[1] - b * range_pct[1]

mimscolor1 = "#0072B2FF"
mimscolor = "black"
pw = diff_df %>%
  # filter(name != "Intercept") %>%
  ggplot(aes(x = s)) +
  geom_line(aes(y = diff_pct, color = diff_type), linewidth = 1.1) +
  geom_line(aes(y = (beta_avg - a) / b), color = mimscolor, linetype = "3212", linewidth = .6) +  # rescale beta for plotting
  geom_hline(aes(yintercept = 0), color = "#666666", linetype = "dashed") +
  facet_wrap(name ~ ., ncol = 4) +
  scale_y_continuous(
    name = "% Difference in CI width",
    breaks = seq(-30, 120, 30),
    sec.axis = sec_axis(~ a + b * ., name = expression("Estimated " * hat(beta) * " from weighted method")) # reverse transform
  ) +
  scale_color_manual(values = c(col1, col2, mimscolor1), labels = c("BRR minus Unweighted", "BRR minus Weighted", "Weighted minus Unweighted"), name = "") +
  scale_x_continuous(breaks= c(0, seq(240, 1440, 6 * 60)), labels = c("", "06:00", "12:00",  "18:00", "24:00")) +
  labs(x = "Hour of Day") +
  theme(strip.text = element_text(size = 14),
        axis.title.y.right = element_text(color = mimscolor)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 14),
                   title = element_text(size = 18)) +
  theme_sub_axis_x(text = element_text(size = 14),
                  # text = element_text(size = 11, angle = 10),
                   title = element_text(size = 18)) +
  theme_sub_axis_y(text = element_text(size = 14),
                   title = element_text(size = 18))
pw = diff_df %>%
  filter(diff_type != "weighted_minus_unweighted") %>%
  # filter(name != "Intercept") %>%
  ggplot(aes(x = s)) +
  geom_line(aes(y = diff_pct, color = diff_type), linewidth = 1.1) +
  geom_line(aes(y = (beta_avg - a) / b), color = mimscolor, linetype = "3212", linewidth = .6) +  # rescale beta for plotting
  geom_hline(aes(yintercept = 0), color = "#666666", linetype = "dashed") +
  facet_wrap(name ~ ., ncol = 4) +
  scale_y_continuous(
    name = "% Difference in CI width",
    breaks = seq(-30, 120, 30),
    sec.axis = sec_axis(~ a + b * ., name = expression("Estimated " * hat(beta) * " from weighted method")) # reverse transform
  ) +
  scale_color_manual(values = c(col1, col2), labels = c("BRR minus Unweighted", "BRR minus Weighted", "Weighted minus Unweighted"), name = "") +
  scale_x_continuous(breaks= c(0, seq(240, 1440, 6 * 60)), labels = c("", "06:00", "12:00",  "18:00", "24:00")) +
  labs(x = "Hour of Day") +
  theme(strip.text = element_text(size = 14),
        axis.title.y.right = element_text(color = mimscolor)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 14),
                   title = element_text(size = 18)) +
  theme_sub_axis_x(text = element_text(size = 14),
                   # text = element_text(size = 11, angle = 10),
                   title = element_text(size = 18)) +
  theme_sub_axis_y(text = element_text(size = 14),
                   title = element_text(size = 18))

png(here::here("manuscript", "figures", "mims_application_width.png"),
    width = 12, height = 8, res = 400, units = "in")
pw
dev.off()

### create DF that is difference in beta coefs

diff_df_beta =
  df %>%
  select(modelf, s, name,  beta) %>%
  pivot_wider(names_from = modelf, values_from = beta) %>%
  mutate(weighted_minus_unweighted_pct = abs((Weighted - Unweighted)) / abs(Unweighted) * 100,
        weighted_minus_unweighted = Weighted - Unweighted)


diff_df_beta %>%
  group_by(name) %>%
  summarize(across(weighted_minus_unweighted, list(min = min, max = max, mean = mean, median = median), .names = "{col}_{fn}"))

range_pct  = range(diff_df_beta$weighted_minus_unweighted, na.rm = TRUE)
range_beta  = range(diff_df_beta$Weighted, na.rm = TRUE)

b = diff(range_beta) / diff(range_pct)
a = range_beta[1] - b * range_pct[1]

mimscolor1 = "#0072B2FF"
mimscolor = "black"
pw = diff_df %>%
  filter(diff_type != "weighted_minus_unweighted") %>%
  # filter(name != "Intercept") %>%
  ggplot(aes(x = s)) +
  geom_line(aes(y = diff_pct, color = diff_type), linewidth = 1.1) +
  geom_line(aes(y = (beta_avg - a) / b), color = mimscolor, linetype = "3212", linewidth = .6) +  # rescale beta for plotting
  geom_hline(aes(yintercept = 0), color = "#666666", linetype = "dashed") +
  facet_wrap(name ~ ., ncol = 4) +
  scale_y_continuous(
    name = "% Difference in CI width",
    breaks = seq(-30, 120, 30),
    sec.axis = sec_axis(~ a + b * ., name = expression("Estimated " * hat(beta) * " across methods")) # reverse transform
  ) +
  scale_color_manual(values = c(col1, col2), labels = c("BRR minus Unweighted", "BRR minus Weighted", "Weighted minus Unweighted"), name = "") +
  scale_x_continuous(breaks= c(0, seq(240, 1440, 6 * 60)), labels = c("", "06:00", "12:00",  "18:00", "24:00")) +
  labs(x = "Hour of Day") +
  theme(strip.text = element_text(size = 14),
        axis.title.y.right = element_text(color = mimscolor)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 14),
                   title = element_text(size = 18)) +
  theme_sub_axis_x(text = element_text(size = 14),
                   # text = element_text(size = 11, angle = 10),
                   title = element_text(size = 18)) +
  theme_sub_axis_y(text = element_text(size = 14),
                   title = element_text(size = 18))


df =
  plt_df_boot %>%
  bind_rows(.id = "model") %>%
  as_tibble() %>%
  filter(model %in% c("1", "3")) %>%
  mutate(modelf = recode_factor(model,
                                "1" = "BRR",
                                # "2" = "Weighted",
                                "3" = "Unweighted")) %>%
  mutate(name = factor(name,
                       levels = c("(Intercept)", "age_cat[30,50)",   "age_cat[50,65)",   "age_cat[65,80]",  "sex_bin", "race_catNH Black",
                                  "race_catHispanic", "race_catOther"),
                       labels = c("Intercept", "Age 30-49", "Age 50-64", "Age 65+", "Sex: Female", "Non-Hispanic Black", "Hispanic", "Other race")))  %>%
  mutate(modelf = fct_rev(modelf))

extra_rows =
  tibble(name = c("(Intercept)", "age_cat[30,50)", "age_cat[50,65)", "age_cat[65,80]", "sex_bin", "race_catNH Black", "race_catHispanic",
                  "race_catOther"),
         min =  c(0, rep(-6, 7)),
         max = c(0, rep(4, 7))) %>%
  mutate(name = factor(name,
                       levels = c("(Intercept)", "age_cat[30,50)",   "age_cat[50,65)",   "age_cat[65,80]",  "sex_bin", "race_catNH Black",
                                  "race_catHispanic", "race_catOther"),
                       labels = c("Intercept", "Age 30-49", "Age 50-64", "Age 65+", "Sex: Female", "Non-Hispanic Black", "Hispanic", "Other race")))

extra_rows =
  extra_rows %>% mutate(modelf = "BRR") %>%
  bind_rows(extra_rows %>% mutate(modelf = "Unweighted"))  %>%
  mutate(modelf = factor(modelf, levels = c("Unweighted", "BRR")))

plt_final=
  df %>%
  ggplot(aes(x = s, y = beta, group = modelf)) +
  facet_grid(name ~  modelf, scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = modelf), alpha = .5) +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = modelf), alpha = .3) +
  geom_line(aes(color = modelf), linewidth = 1.1)  +
  geom_point(data = extra_rows, aes(x = 60, y = min), alpha = 0) +
  geom_point(data = extra_rows, aes(x = 60, y = max), alpha = 0) +
  theme_grey() +
  geom_hline(
    data = tibble(name = factor(setdiff(levels(df$name), "Intercept")), levels = c("Age 30-49", "Age 50-64", "Age 65+", "Sex: Female", "Non-Hispanic Black", "Hispanic", "Other race")),
    aes(yintercept = 0),
    linetype = "dashed"
  ) +
  scale_y_continuous(breaks=seq(-100, 100, 3)) +
  scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
  labs(x = "Hour of Day", y = "Monitor Independent Movement Summary (MIMS)", fill = "Model", color = "Model") +
  theme(strip.text = element_text(size = 12),
        palette.color.discrete = c(col2, col1),
        palette.fill.discrete = c(col2, col1)) +
  theme_sub_legend(position = "none",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(size = 11, angle = 10),
                   title = element_text(size = 13)) +
  theme_sub_axis_y(text = element_text(size = 11),
                   title = element_text(size = 13))

df1 =
  df %>%
  filter(grepl("Int", name) | grepl("Age", name)) %>%
  mutate(name= factor(name,levels =  c("Intercept", "Age 30-49", "Age 50-64", "Age 65+")))

df2 =
  df %>%
  filter(!(grepl("Int", name) | grepl("Age", name))) %>%
  mutate(name= factor(name,levels =  c("Sex: Female", "Non-Hispanic Black", "Hispanic", "Other race")))

plt_final1 =
  df1 %>%
  ggplot(aes(x = s, y = beta, group = modelf)) +
  facet_grid(name ~  modelf, scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = modelf), alpha = .5) +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = modelf), alpha = .3) +
  geom_line(aes(color = modelf), linewidth = 1.1)  +
  geom_point(data = extra_rows %>% filter(grepl("Int", name) | grepl("Age", name)), aes(x = 60, y = min), alpha = 0) +
  geom_point(data = extra_rows %>% filter(grepl("Int", name) | grepl("Age", name)), aes(x = 60, y = max), alpha = 0) +
  theme_grey() +
  geom_hline(
    data = tibble(name = factor(setdiff(levels(df1$name), "Intercept")), levels = c("Age 30-49", "Age 50-64", "Age 65+")),
    aes(yintercept = 0),
    linetype = "dashed"
  ) +
  scale_y_continuous(breaks=seq(-100, 100, 2)) +
  scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
  labs(x = "Hour of Day", y = "Monitor Independent Movement Summary (MIMS)", fill = "Model", color = "Model") +
  theme(strip.text = element_text(size = 12),
        palette.color.discrete = c(col2, col1),
        palette.fill.discrete = c(col2, col1)) +
  theme_sub_legend(position = "none",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(size = 9, angle = 20),
                   title = element_text(size = 13)) +
  theme_sub_axis_y(text = element_text(size = 11),
                   title = element_text(size = 13))

plt_final1 =
  df1 %>%
  ggplot(aes(x = s, y = beta, group = modelf)) +
  facet_grid(name ~  ., scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = modelf), alpha = .5) +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = modelf), alpha = .3) +
  geom_line(aes(color = modelf), linewidth = 1.1)  +
  geom_point(data = extra_rows %>% filter(grepl("Int", name) | grepl("Age", name)), aes(x = 60, y = min), alpha = 0) +
  geom_point(data = extra_rows %>% filter(grepl("Int", name) | grepl("Age", name)), aes(x = 60, y = max), alpha = 0) +
  theme_grey() +
  geom_hline(
    data = tibble(name = factor(setdiff(levels(df1$name), "Intercept")), levels = c("Age 30-49", "Age 50-64", "Age 65+")),
    aes(yintercept = 0),
    linetype = "dashed"
  ) +
  scale_y_continuous(breaks=seq(-100, 100, 2)) +
  scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
  labs(x = "Hour of Day", y = "Monitor Independent Movement Summary (MIMS)", fill = "Model", color = "Model") +
  theme(strip.text = element_text(size = 12),
        palette.color.discrete = c(col2, col1),
        palette.fill.discrete = c(col2, col1)) +
  theme_sub_legend(position = "none",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(size = 9, angle = 20),
                   title = element_text(size = 13)) +
  theme_sub_axis_y(text = element_text(size = 11),
                   title = element_text(size = 13))

plt_final2 =
  df2 %>%
  ggplot(aes(x = s, y = beta, group = modelf)) +
  facet_grid(name ~  modelf, scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = modelf), alpha = .5) +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = modelf), alpha = .3) +
  geom_line(aes(color = modelf), linewidth = 1.1)  +
  geom_point(data = extra_rows %>% filter(!(grepl("Int", name) | grepl("Age", name))), aes(x = 60, y = min), alpha = 0) +
  geom_point(data = extra_rows %>% filter(!(grepl("Int", name) | grepl("Age", name))), aes(x = 60, y = max), alpha = 0) +
  theme_grey() +
  geom_hline(
    aes(yintercept = 0),
    linetype = "dashed"
  ) +
  scale_y_continuous(breaks=seq(-100, 100, 2)) +
  scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
  labs(x = "Hour of Day", y = "", fill = "Model", color = "Model") +
  theme(strip.text = element_text(size = 12),
        palette.color.discrete = c(col2, col1),
        palette.fill.discrete = c(col2, col1)) +
  theme_sub_legend(position = "none",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(size = 9, angle = 20),
                   title = element_text(size = 13)) +
  theme_sub_axis_y(text = element_text(size = 11),
                   title = element_text(size = 13))




png(here::here("manuscript", "figures", "mims_application_wide.png"),
    width = 12, height = 8, res = 400, units = "in")
plt_final1 + plt_final2
dev.off()


disagree_df =
  df1 %>%
  select(lower.joint, upper.joint, name, modelf, s) %>%
  mutate(sig = case_when(
    lower.joint > 0 & upper.joint > 0 ~ 1,
    upper.joint < 0 & lower.joint < 0 ~ 1,
    .default = 0
  )) %>%
  select(name, modelf, s, sig) %>%
  pivot_wider(names_from = modelf, values_from = sig) %>%
  mutate(brr_sig = BRR == 1 & Unweighted == 0,
         # unwt_sig = Unweighted == 1 & BRR == 0,
         disagree = BRR != Unweighted) %>%
  filter(disagree) %>%
  mutate(modelf = "BRR", beta = 1) # for plotting

disagree_df2 =
  df2 %>%
  select(lower.joint, upper.joint, name, modelf, s) %>%
  mutate(sig = case_when(
    lower.joint > 0 & upper.joint > 0 ~ 1,
    upper.joint < 0 & lower.joint < 0 ~ 1,
    .default = 0
  )) %>%
  select(name, modelf, s, sig) %>%
  pivot_wider(names_from = modelf, values_from = sig) %>%
  mutate(brr_sig = BRR == 1 & Unweighted == 0,
         # unwt_sig = Unweighted == 1 & BRR == 0,
         disagree = BRR != Unweighted) %>%
  filter(disagree) %>%
  mutate(modelf = "BRR", beta = 1) # for plotting

plt_final1_ov =
  ggplot(df1, aes(x = s, y = beta, group = modelf)) +
  facet_grid(name ~  ., scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = modelf), alpha = .5) +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = modelf), alpha = .3) +
  geom_line(aes(color = modelf), linewidth = 1.1)  +
  geom_point(data = extra_rows %>% filter(grepl("Int", name) | grepl("Age", name)), aes(x = 60, y = min), alpha = 0) +
  geom_point(data = extra_rows %>% filter(grepl("Int", name) | grepl("Age", name)), aes(x = 60, y = max), alpha = 0) +
  theme_grey() +
  geom_hline(
    data = tibble(name = factor(setdiff(levels(df1$name), "Intercept")), levels = c("Age 30-49", "Age 50-64", "Age 65+")),
    aes(yintercept = 0),
    linetype = "dashed"
  ) +
  scale_y_continuous(breaks=seq(-100, 100, 2)) +
  scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
  labs(x = "Hour of Day", y = "Monitor Independent Movement Summary (MIMS)", fill = "Model", color = "Model") +
  theme(strip.text = element_text(size = 12),
        palette.color.discrete = c(col2, col1),
        palette.fill.discrete = c(col2, col1)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(size = 9, angle = 20),
                   title = element_text(size = 13)) +
  theme_sub_axis_y(text = element_text(size = 11),
                   title = element_text(size = 13))

plt_final1_ov_v2 =
  ggplot(df1, aes(x = s, y = beta, group = modelf)) +
    facet_grid(name ~  ., scales = "free_y") +
    geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = modelf), alpha = .5) +
    geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = modelf), alpha = .3) +
    geom_line(aes(color = modelf), linewidth = 1.1)  +
    geom_point(data = extra_rows %>% filter(grepl("Int", name) | grepl("Age", name)), aes(x = 60, y = min), alpha = 0) +
    geom_point(data = extra_rows %>% filter(grepl("Int", name) | grepl("Age", name)), aes(x = 60, y = max), alpha = 0) +
    theme_grey() +
    geom_hline(
      data = tibble(name = factor(setdiff(levels(df1$name), "Intercept")), levels = c("Age 30-49", "Age 50-64", "Age 65+")),
      aes(yintercept = 0),
      linetype = "dashed"
    ) +
    scale_y_continuous(breaks=seq(-100, 100, 2)) +
    scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
    labs(x = "Hour of Day", y = "Monitor Independent Movement Summary (MIMS)", fill = "Model", color = "Model") +
    theme(strip.text = element_text(size = 12),
          palette.color.discrete = c(col2, col1),
          palette.fill.discrete = c(col2, col1)) +
    theme_sub_legend(position = "bottom",
                     text = element_text(size = 12),
                     title = element_text(size = 13)) +
    theme_sub_axis_x(text = element_text(size = 11),
                     title = element_text(size = 13)) +
    theme_sub_axis_y(text = element_text(size = 11),
                     title = element_text(size = 13)) +
    geom_rect(data = disagree_df %>% filter(brr_sig),
              aes(x = s, y = -6, height = 1, width = 1), fill = "#E51932") +
    geom_rect(data = disagree_df %>% filter(!brr_sig),
              aes(x = s, y = -6, height = 1, width = 1), fill = "#19B2FF")
    # RED: BRR signif, unwtd not signif
    # BLUE: unwtd signif, BRR not

pfov = df2 %>%
  ggplot(aes(x = s, y = beta, group = modelf)) +
  facet_grid(name ~ ., scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = modelf), alpha = .5) +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = modelf), alpha = .3) +
  geom_line(aes(color = modelf), linewidth = 1.1)  +
  geom_point(data = extra_rows %>% filter(!(grepl("Int", name) | grepl("Age", name))), aes(x = 60, y = min), alpha = 0) +
  geom_point(data = extra_rows %>% filter(!(grepl("Int", name) | grepl("Age", name))), aes(x = 60, y = max), alpha = 0) +
  theme_grey() +
  geom_hline(
    aes(yintercept = 0),
    linetype = "dashed"
  ) +
  scale_y_continuous(breaks=seq(-100, 100, 2)) +
  scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
  labs(x = "Hour of Day", y = "", fill = "Model", color = "Model") +
  theme(strip.text = element_text(size = 12),
        palette.color.discrete = c(col2, col1),
        palette.fill.discrete = c(col2, col1)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(size = 9, angle = 20),
                   title = element_text(size = 13)) +
  theme_sub_axis_y(text = element_text(size = 11),
                   title = element_text(size = 13))

pfov_v2 = df2 %>%
  ggplot(aes(x = s, y = beta, group = modelf)) +
  facet_grid(name ~ ., scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = modelf), alpha = .5) +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = modelf), alpha = .3) +
  geom_line(aes(color = modelf), linewidth = 1.1)  +
  geom_point(data = extra_rows %>% filter(!(grepl("Int", name) | grepl("Age", name))), aes(x = 60, y = min), alpha = 0) +
  geom_point(data = extra_rows %>% filter(!(grepl("Int", name) | grepl("Age", name))), aes(x = 60, y = max), alpha = 0) +
  theme_grey() +
  geom_hline(
    aes(yintercept = 0),
    linetype = "dashed"
  ) +
  scale_y_continuous(breaks=seq(-100, 100, 2)) +
  scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
  labs(x = "Hour of Day", y = "", fill = "Model", color = "Model") +
  theme(strip.text = element_text(size = 12),
        palette.color.discrete = c(col2, col1),
        palette.fill.discrete = c(col2, col1)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(size = 11),
                   title = element_text(size = 13)) +
  theme_sub_axis_y(text = element_text(size = 11),
                   title = element_text(size = 13)) +
  geom_rect(data = disagree_df2 %>% filter(brr_sig),
            aes(x = s, y = -6, height = 1, width = 1), fill = "#E51932") +
  geom_rect(data = disagree_df2 %>% filter(!brr_sig),
            aes(x = s, y = -6, height = 1, width = 1), fill = "#19B2FF")


png(here::here("manuscript", "figures", "mims_application_wide_overlap.png"),
    width = 12, height = 8, res = 400, units = "in")
plt_final1_ov + pfov + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()

png(here::here("manuscript", "figures", "mims_application_wide_highlight.png"),
    width = 12, height = 8, res = 400, units = "in")
plt_final1_ov_v2 + pfov_v2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()

svg(here::here("manuscript", "figures", "mims_application_wide_highlight.svg"),
    width = 12, height = 8)
plt_final1_ov_v2 + pfov_v2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()



########### old plots, don't use


p1 = plt_df_boot %>%
  bind_rows(.id = "model") %>%
  filter(model %in% c("1", "3")) %>%
  filter(grepl("age", name)) %>%
  mutate(model = factor(model, labels = c("BRR", "Unweighted")),
         model = factor(model, levels = c("Unweighted", "BRR"))) %>%
  mutate(name = factor(name, labels = paste0("Age ", c("30-49", "50-64", "65+")))) %>%
  # mutate(across(c(beta, lower, upper, lower.joint, upper.joint), ~exp(.x))) %>%
  mutate(model = factor(model, labels = c("BRR", "Unweighted")),
         model = factor(model, levels = c("Unweighted", "BRR")))%>%
  ggplot(aes(x = s, y = beta, group = model)) +
  facet_grid(.~name, scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = model), alpha = .3) +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = model), alpha = .1) +
  geom_line(aes(color = model, linetype = model), linewidth = 1.1) +
  theme_classic() +
  geom_hline(aes(yintercept = 0)) +
  scale_y_continuous(breaks=seq(-100, 100, 2)) +
  scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
  labs(x = "Hour of Day", y = "Monitor Independent Movement Summary (MIMS)", fill = "Model", color = "Model", linetype = "Model") +
  theme(strip.text = element_text(size = 12),
        palette.color.discrete = c(col2, col1),
        palette.fill.discrete = c(col2, col1)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(size = 11, angle = 10),
                   title = element_text(size = 13)) +
  theme_sub_axis_y(text = element_text(size = 11),
                   title = element_text(size = 13))


p2 = plt_df_boot %>%
  bind_rows(.id = "model") %>%
  filter(model %in% c("1", "3")) %>%
  filter(grepl("race", name)) %>%
  mutate(name = factor(name, labels = c("Hispanic", "Non-Hispanic Black", "Other"))) %>%
  # mutate(across(c(beta, lower, upper, lower.joint, upper.joint), ~exp(.x))) %>%
  mutate(model = factor(model, labels = c("BRR", "Unweighted")),
         model = factor(model, levels = c("Unweighted", "BRR"))) %>%
  ggplot(aes(x = s, y = beta, group = model)) +
  # geom_line(aes(color = model), linewidth = 1.1)  +
  facet_grid(.~name, scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = model), alpha = .3) +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = model), alpha = .1) +
  geom_line(aes(color = model, linetype = model), linewidth = 1.1)  +
  theme_classic() +
  geom_hline(aes(yintercept = 0)) +
  scale_y_continuous(breaks=seq(-100, 100, 2)) +
  scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
  labs(x = "Hour of Day", y = "Monitor Independent Movement Summary (MIMS)", fill = "Model", color = "Model", linetype = "Model") +
  theme(strip.text = element_text(size = 12),
        palette.color.discrete = c(col2, col1),
        palette.fill.discrete = c(col2, col1)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(size = 11, angle = 10),
                   title = element_text(size = 13)) +
  theme_sub_axis_y(text = element_text(size = 11),
                   title = element_text(size = 13))

p3 = plt_df_boot %>%
  bind_rows(.id = "model") %>%
  filter(model %in% c("1", "3")) %>%
  filter(grepl("sex", name)) %>%
  mutate(name = factor(name, labels = c("Sex: Female"))) %>%
  mutate(model = factor(model, labels = c("BRR", "Unweighted")),
         model = factor(model, levels = c("Unweighted", "BRR"))) %>%
  ggplot(aes(x = s, y = beta, group = model)) +
  facet_wrap(.~name, scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = model), alpha = .3) +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = model), alpha = .1) +
  geom_line(aes(color = model, linetype = model), linewidth = 1.1)  +
  theme_classic() +
  geom_hline(aes(yintercept = 0)) +
  scale_y_continuous(breaks=seq(-100, 100, .5)) +
  scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
  labs(x = "Hour of Day", y = "Monitor Independent Movement Summary (MIMS)", fill = "Model", color = "Model", linetype = "Model") +
  theme(strip.text = element_text(size = 12),
        palette.color.discrete = c(col2, col1),
        palette.fill.discrete = c(col2, col1)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(size = 11, angle = 10),
                   title = element_text(size = 13)) +
  theme_sub_axis_y(text = element_text(size = 11),
                   title = element_text(size = 13))

p4 = plt_df_boot %>%
  bind_rows(.id = "model") %>%
  filter(model %in% c("1", "3")) %>%
  # mutate(across(c(beta, lower, upper, lower.joint, upper.joint), ~exp(.x))) %>%
  mutate(model = factor(model, labels = c("BRR", "Unweighted")),
         model = factor(model, levels = c("Unweighted", "BRR"))) %>%
  filter(grepl("cept", name)) %>%
  mutate(name = factor(name, labels = c("Intercept"))) %>%
  # mutate(across(c(beta, lower, upper, lower.joint, upper.joint), ~exp(.x))) %>%
  ggplot(aes(x = s, y = beta, group = model)) +
  facet_wrap(.~name, scales = "free_y") +
  geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = model), alpha = .3) +
  geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = model), alpha = .1) +
  geom_line(aes(color = model, linetype = model), linewidth = 1.1)  +
  theme_classic() +
  scale_y_continuous(breaks=seq(-100, 100, 2)) +
  scale_x_continuous(breaks= c(60, seq(240, 1440, 4 * 60)), labels = c("01:00", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
  labs(x = "Hour of Day", y = "Monitor Independent Movement Summary (MIMS)", fill = "Model", color = "Model", linetype = "Model") +
  theme(strip.text = element_text(size = 12),
        palette.color.discrete = c(col2, col1),
        palette.fill.discrete = c(col2, col1)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(size = 11, angle = 10),
                   title = element_text(size = 13)) +
  theme_sub_axis_y(text = element_text(size = 11),
                   title = element_text(size = 13))

png(here::here("manuscript", "figures", "mims_application_all.png"),
    width = 12, height = 10, res = 300, units = "in")
(p4 | p3) / p1 / p2  + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
dev.off()
