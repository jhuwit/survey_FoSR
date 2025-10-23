## work on real world data result
library(tidyverse)
pa_df = read_rds(here::here("data", "mims_covariates.rds")) %>%
  mutate(SEQN = as.character(SEQN))
source(here::here("R_cp", "01_sim_functions.R"))
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


cis = read_rds(here::here("results", "mims_application_1440.rds"))

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

# check CI widths
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




########### old plots


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
library(patchwork)

png(here::here("manuscript", "figures", "mims_application_all.png"),
    width = 12, height = 10, res = 300, units = "in")
(p4 | p3) / p1 / p2  + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
dev.off()
