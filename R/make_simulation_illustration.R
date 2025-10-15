source(here::here("R", "00_data_gen_function_ff.R"))


lst = generate_superpopulation(I = 10e6, L = 50,
                         scenario = 1,
                         family = "gaussian",
                         seed = 4565,
                         num_strata = 30,
                         strata_sigma = 0.05,
                         strata_scale = 0.125,
                         snr_b = 1, snr_eps = 1)

lst2 = generate_superpopulation(I = 10e6, L = 50,
                               scenario = 1,
                               family = "gaussian",
                               seed = 4565,
                               num_strata = 30,
                               strata_sigma = 0,
                               strata_scale = 0,
                               snr_b = 1, snr_eps = 1)


psu_df = expand_grid(strata = 1:10,
                     psu = 1:5) %>%
  mutate(x = paste0(strata, "_", psu))

lst$Y_obs %>%
  as_tibble() %>%
  bind_cols(psu = lst$psu_assignments,
            strata = lst$stratum_assignments) %>%
  filter(strata %in% c(1:10),
         psu %in% psu_df$x) %>%
  group_by(strata, psu) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "ind",
               names_transform = ~as.integer(sub(".*V", "", .x))) %>%
  mutate(strata = paste0("Stratum ", strata),
         strata = factor(strata, levels = paste0("Stratum ", seq(1:10)))) %>%
  ggplot(aes(x = ind, y = value, color = factor(psu))) +
  facet_wrap(~strata, nrow = 2) +
  geom_smooth(se = FALSE, linewidth = 0.7) +
  labs(x = "Functional Domain", y = "Value", title = "Mean smoothed outcomes by stratum for 10 strata") +
  theme(legend.position = "none") +
  scale_color_viridis_d(option = "A")

lst2$Y_obs %>%
  as_tibble() %>%
  bind_cols(psu = lst2$psu_assignments,
            strata = lst2$stratum_assignments) %>%
  filter(strata %in% c(1:10),
         psu %in% psu_df$x) %>%
  group_by(strata, psu) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "ind",
               names_transform = ~as.integer(sub(".*V", "", .x))) %>%
  mutate(strata = paste0("Stratum ", strata),
         strata = factor(strata, levels = paste0("Stratum ", seq(1:10)))) %>%
  ggplot(aes(x = ind, y = value, color = factor(psu))) +
  facet_wrap(~strata, nrow = 2) +
  geom_smooth(se = FALSE, linewidth = 0.7) +
  labs(x = "Functional Domain", y = "Value", title = "Mean smoothed outcomes by stratum for 10 strata") +
  theme(legend.position = "none") +
  scale_color_viridis_d(option = "A")

lst2$Y_obs %>%
  as_tibble() %>%
  bind_cols(psu = lst2$psu_assignments,
            strata = lst2$stratum_assignments) %>%
  filter(strata %in% c(1:10),
         psu %in% psu_df$x) %>%
  group_by(strata, psu) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "ind",
               names_transform = ~as.integer(sub(".*V", "", .x))) %>%
  mutate(strata = paste0("Stratum ", strata),
         strata = factor(strata, levels = paste0("Stratum ", seq(1:10)))) %>%
  ggplot(aes(x = ind, y = value, color = factor(psu))) +
  facet_wrap(~strata, nrow = 2) +
  geom_smooth(se = FALSE, linewidth = 0.7) +
  labs(x = "Functional Domain", y = "Value", title = "Mean smoothed outcomes by stratum for 10 strata") +
  theme(legend.position = "none") +
  scale_color_viridis_d(option = "A")

no_re =
  lst2$Y_obs %>%
  as_tibble() %>%
  bind_cols(psu = lst2$psu_assignments,
            strata = lst2$stratum_assignments) %>%
  filter(strata %in% c(1:10),
         psu %in% psu_df$x) %>%
  group_by(strata, psu) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "ind",
               names_transform = ~as.integer(sub(".*V", "", .x))) %>%
  mutate(type = "No random effects")

ss = lst$Y_obs %>%
  as_tibble() %>%
  bind_cols(psu = lst$psu_assignments,
            strata = lst$stratum_assignments) %>%
  filter(strata %in% c(1:10),
         psu %in% psu_df$x) %>%
  group_by(strata, psu) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  pivot_longer(cols = starts_with("V"),
               names_to = "ind",
               names_transform = ~as.integer(sub(".*V", "", .x))) %>%
  mutate(type = "Stratum scaling and noise")




p1 = no_re %>% bind_rows(ss) %>%
  mutate(strata = paste0("Stratum ", strata),
         strata = factor(strata, levels = paste0("Stratum ", seq(1:10)))) %>%
  filter(strata %in% c("Stratum 1", "Stratum 2", "Stratum 3")) %>%
  ggplot(aes(x = ind, y = value, color = factor(psu))) +
  facet_grid(strata ~ type) +
  geom_smooth(se = FALSE, linewidth = 1) +
  labs(x = "Functional Domain", y = "Value", title = "Simulation",
       subtitle = "Mean smoothed outcomes by stratum and PSU") +
  theme_light(base_size = 14) +
  # theme(legend.position = "none") +
  scale_color_paletteer_d("colorBlindness::SteppedSequential5Steps", direction = -1,
                          name = "Stratum_PSU")


  # scale_color_viridis_d(option = "B", name = "stratum_psu")

pa_df = read_rds(here::here("data", "pa_df_persub.rds"))

strata_psu =
  pa_df %>% select(strata, psu) %>%
  distinct() %>%
  arrange(desc(strata)) %>%
  mutate(s2 = data.table::rleid(strata)) %>%
  filter(s2 <= 3)

pa_df_small =
  pa_df %>%
  filter(strata %in% strata_psu$strata,
         psu %in% strata_psu$psu) %>%
  pivot_longer(cols = starts_with("min"),
               names_to = "ind",
               names_transform = ~as.integer(sub(".*min\\_", "", .x)))

p2 = pa_df_small %>%
  mutate(strata = factor(strata, labels = c("1", "2", "3")),
          psu = paste0(strata, "_", psu),
         strata = paste0("Stratum ", strata)) %>%
  ggplot(aes(x = ind, y = value, color = factor(psu))) +
  facet_grid(strata ~ .) +
  geom_smooth(se = FALSE, linewidth = 1) +
  labs(x = "Functional Domain", y = "MIMS per minute", title = "NHANES",
  subtitle = "Mean smoothed outcomes by stratum and PSU") +
  theme_light(base_size = 14) +
  # theme(legend.position = "none") +
  scale_color_manual(values = c("#990F0FFF", "#CC5151FF", "#99540FFF", "#CC8E51FF", "#FF5500", "#FFCC65"),
                     name = "Stratum_PSU") +
  scale_x_continuous(limits = c(0,1440), breaks = seq(0, 1440, 240), labels = c("00:00", "04:00", "08:00",
                                                                                "12:00", "16:00", "20:00", "24:00"))

library(patchwork)

png(here::here("manuscript", "figures", "outcomes_strata_psu.png"), width = 14, height = 6, units = "in", res = 350)
p1 + p2 + plot_annotation(tag_levels = "A")
dev.off()

### plot outcomes by weight tertile

sample = sample_from_population_wor(
  X_des = lst$X_des,
  Y_obs = lst$Y_obs,
  L = 50,
  I_n = 100,
  num_strata = 30,
  stratum_assignments = lst$stratum_assignments,
  psu_assignments = lst$psu_assignments,
  dirichlet_probs = lst$dirichlet_probs,
  seed = iter,
  inf_level = 15,
  compression = 2,
  family = "gaussian"
)

sample_noninf = sample_from_population_wor(
  X_des = lst$X_des,
  Y_obs = lst$Y_obs,
  L = 50,
  I_n = 100,
  num_strata = 30,
  stratum_assignments = lst$stratum_assignments,
  psu_assignments = lst$psu_assignments,
  dirichlet_probs = lst$dirichlet_probs,
  seed = iter,
  inf_level = 0,
  compression = 2,
  family = "gaussian"
)


Y_mat = sample %>% select(starts_with("Y")) %>%
  as.matrix()

quants = quantile(sample$weight, c(1/3, 2/3))

mean_df =
  sample %>%
  mutate(Y = I(Y_mat)) %>%
  mutate(quant = cut(weight, breaks = c(-Inf, quants, Inf), labels = c("low", "med", "high"))) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50),
         Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(quant) %>%
  summarize(Y_mean_strat = mean(Y_tf)) %>%
  ungroup() %>%
  mutate(Y_mean_strat = tf_smooth(Y_mean_strat, method = "lowess"))

Y_mat_ni = sample_noninf %>% select(starts_with("Y")) %>%
  as.matrix()
quants = quantile(sample_noninf$weight, c(1/3, 2/3))

mean_df_ni =
  sample_noninf %>%
  mutate(Y = I(Y_mat_ni)) %>%
  mutate(quant = cut(weight, breaks = c(-Inf, quants, Inf), labels = c("low", "med", "high"))) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50),
         Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(quant) %>%
  summarize(Y_mean_strat = mean(Y_tf)) %>%
  ungroup() %>%
  mutate(Y_mean_strat = tf_smooth(Y_mean_strat, method = "lowess"))


p1 = mean_df %>%
  mutate(sample = "Sampling: informative") %>%
  bind_rows(mean_df_ni %>% mutate(sample ="Sampling: random")) %>%
  ggplot() +
  geom_spaghetti(aes(y = Y_mean_strat, color = quant), linewidth = 1.1, alpha = 1) +
  scale_color_manual(name = "Weight tertile",
                     values = c("#075AFFFF", "#FF9900FF", "#FF0000FF"),
                     labels = c("1", "2", "3")) +
  facet_grid(.~sample)+
  labs(x = "Functional domain", y = "Value",
       title = "Simulation",
       subtitle = "Mean outcomes by weight tertile") +
  theme_light(base_size = 14) +
  theme(legend.position = "bottom")



Y_mat_nh = pa_df %>% select(starts_with("min")) %>%
  as.matrix()
quants = quantile(pa_df$full_sample_2_year_mec_exam_weight/2, c(1/3, 2/3))

mean_df_nh =
  pa_df %>%
  select(-starts_with("min")) %>%
  mutate(weight = full_sample_2_year_mec_exam_weight/2) %>%
  mutate(Y = I(Y_mat_nh)) %>%
  mutate(quant = cut(weight, breaks = c(-Inf, quants, Inf), labels = c("low", "med", "high"))) %>%
  mutate(Y = matrix(Y, ncol = 1440),
         Y_tf = tfd(Y, length = 1440),
         Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(quant) %>%
  summarize(Y_mean_strat = mean(Y_tf)) %>%
  ungroup() %>%
  mutate(Y_mean_strat = tf_smooth(Y_mean_strat, method = "lowess"))

# paletteer_d("colorBlindness::Blue2OrangeRed14Steps")



p2 = mean_df_nh %>%
  ggplot() +
  geom_spaghetti(aes(y = Y_mean_strat, color = quant, linetype = quant), linewidth = 1.1, alpha = 1) +
  scale_color_manual(name = "Weight tertile",
                     values = c("#075AFFFF", "#FF9900FF", "#FF0000FF"),
                     labels = c("1", "2", "3")) +
  scale_linetype_manual(values = c(1, 2, 3),
                        labels = c("1", "2", "3"),
                        name = "Weight tertile") +
  labs(x = "Functional domain", y = "Value",
       title = "NHANES",
       subtitle = "Mean outcomes by weight tertile") +
  labs(x = "Functional Domain", y = "MIMS per minute") +
  theme_light(base_size = 14) +
  scale_x_continuous(limits = c(0,1440), breaks = seq(0, 1440, 240), labels = c("00:00", "04:00", "08:00",
                                                                                "12:00", "16:00", "20:00", "24:00")) +
  theme(legend.position = "bottom")

png(here::here("manuscript", "figures", "weight_strata_psu.png"), width = 10, height = 6, units = "in", res = 350)
p1 | p2 + plot_annotation(tag_levels = "A")
dev.off()

## get non-informative from no res

sample2 = sample_from_population_wor(
  X_des = lst2$X_des,
  Y_obs = lst2$Y_obs,
  L = 50,
  I_n = 100,
  num_strata = 30,
  stratum_assignments = lst2$stratum_assignments,
  psu_assignments = lst2$psu_assignments,
  dirichlet_probs = lst2$dirichlet_probs,
  seed = iter,
  inf_level = 15,
  compression = 2,
  family = "gaussian"
)

sample_noninf2 = sample_from_population_wor(
  X_des = lst2$X_des,
  Y_obs = lst2$Y_obs,
  L = 50,
  I_n = 100,
  num_strata = 30,
  stratum_assignments = lst2$stratum_assignments,
  psu_assignments = lst2$psu_assignments,
  dirichlet_probs = lst2$dirichlet_probs,
  seed = iter,
  inf_level = 0,
  compression = 2,
  family = "gaussian"
)


Y_mat = sample2 %>% select(starts_with("Y")) %>%
  as.matrix()

quants = quantile(sample2$weight, c(1/3, 2/3))

mean_df2 =
  sample2 %>%
  mutate(Y = I(Y_mat)) %>%
  mutate(quant = cut(weight, breaks = c(-Inf, quants, Inf), labels = c("low", "med", "high"))) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50),
         Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(quant) %>%
  summarize(Y_mean_strat = mean(Y_tf)) %>%
  ungroup() %>%
  mutate(Y_mean_strat = tf_smooth(Y_mean_strat, method = "lowess"))

Y_mat_ni = sample_noninf2 %>% select(starts_with("Y")) %>%
  as.matrix()
quants = quantile(sample_noninf2$weight, c(1/3, 2/3))

mean_df_ni2 =
  sample_noninf2 %>%
  mutate(Y = I(Y_mat_ni)) %>%
  mutate(quant = cut(weight, breaks = c(-Inf, quants, Inf), labels = c("low", "med", "high"))) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50),
         Y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(quant) %>%
  summarize(Y_mean_strat = mean(Y_tf)) %>%
  ungroup() %>%
  mutate(Y_mean_strat = tf_smooth(Y_mean_strat, method = "lowess"))


p1 = mean_df %>%
  mutate(sample = "Sampling: informative", re = "Strata scaling and noise") %>%
  bind_rows(mean_df_ni %>% mutate(sample ="Sampling: random", re = "Strata scaling and noise")) %>%
  bind_rows(mean_df_ni2 %>% mutate(sample ="Sampling: random", re = "No random effects")) %>%
  bind_rows(mean_df2 %>% mutate(sample ="Sampling: informative", re = "No random effects")) %>%
  ggplot() +
  geom_spaghetti(aes(y = Y_mean_strat, color = quant, linetype = quant), linewidth = 1.1, alpha = 1) +
  scale_color_manual(name = "Weight tertile",
                     values = c("#075AFFFF", "#FF9900FF", "#FF0000FF"),
                     labels = c("1", "2", "3")) +
  scale_linetype_manual(values = c(1, 2, 3),
                        labels = c("1", "2", "3"),
                        name = "Weight tertile") +
  facet_grid(sample ~ re)+
  labs(x = "Functional domain", y = "Value",
       title = "Simulation",
       subtitle = "Mean outcomes by weight tertile") +
  theme_light(base_size = 14) +
  theme(legend.position = "bottom")

png(here::here("manuscript", "figures", "weight_strata_psu.png"), width = 10, height = 6, units = "in", res = 350)
p1 + p2 + plot_annotation(tag_levels = "A")
dev.off()

