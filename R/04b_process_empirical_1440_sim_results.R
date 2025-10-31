# getremote survey_FoSR/results/simulations/empirical_sim_1440/*.rds results/simulations/empirical_sim_1440/

library(tidyverse)
library(patchwork)
nsim = 200
cols = c("#E69F00FF",  "#56B4E9FF", "#009E73FF", "#CC79A7FF")
L = 1440

if (!file.exists(here::here("results", "simulations", "all_empirical_sim_res.rds")) || force) {
  settings = expand_grid(iter = 1:nsim,
                        weight_type = c("uniform", "nh_weights", "mims_weights", "combo_weights")) %>%
    mutate(rn = row_number())

  files = list.files(here::here("results", "simulations", "empirical_sim_1440"), full.names = TRUE)
  res_df = map(.x = files,
               .f = \(x) read_rds(x) %>% mutate(fold = sub(".*fold\\_(.+).rds.*", "\\1", x))) %>%
    bind_rows() %>%
    mutate(fold = as.numeric(fold))

  folds = unique(res_df$fold)
  all.equal(folds, 1:800) # check


  res_df =
    res_df %>%
    select(-iter) %>%
    left_join(settings, by = c("fold" = "rn")) %>%
    rename(setting = weight_type)

  write_rds(res_df, here::here("results", "simulations", "all_empirical_sim_res.rds"), compress = "xz")

} else {
  res_df = read_rds(here::here("results", "simulations", "all_empirical_sim_res.rds"))
}


res_summ =
  res_df %>%
  group_by(boot_type, setting, var) %>%
  mutate(cover_pw = cover_pw / L) %>%
  summarize(across(c(MISE, mean_joint_se, mean_pw_se, cover_pw), mean),
            across(cover_joint, ~sum(.x) / n()),
            .groups = "drop")

p1 = res_df %>%
  mutate(cover_pw = cover_pw / L) %>%
  mutate(boot_type = factor(boot_type, levels = c("unweighted", "weighted", "BRR"), labels = c("Unweighted", "Weighted", "BRR")),
         setting = factor(setting, levels = c("uniform", "nh_weights", "mims_weights", "combo_weights"),
                          labels = c("Uniform", "NHANES", "MIMS", "NHANES&MIMS")),
         var = factor(var, labels = c("Intercept", "Sex Coefficient"))) %>%
  ggplot(aes(x = setting, y = cover_pw, color = boot_type, fill = boot_type)) +
  facet_grid(. ~ var) +
  geom_boxplot(
    box.colour = "black",
    median.linewidth = 1,
    median.color = "black",
    staplewidth = 0.5, # show staple
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               aes(shape = boot_type),
               size = 2.5,
               color = "black",
               position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c(8, 1, 2)) +
  theme(palette.color.discrete = cols[1:3],
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  scale_y_continuous(breaks = c(seq(0.2, 0.8, 0.2), 0.95, 1)) +
  labs(x = "Weighting Scheme", y = "Pointwise Coverage", color = "Inference Type", fill = "Inference Type", shape = "Inference Type") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  geom_hline(aes(yintercept = 0.95), color = cols[4], linetype = "dashed", linewidth = 1.2) +
  theme_sub_axis_x(text = element_text(angle = 30))


p2 = res_df %>%
  mutate(MISE = log(MISE)) %>%
  mutate(boot_type = factor(boot_type, levels = c("unweighted", "weighted", "BRR"), labels = c("Unweighted", "Weighted", "BRR")),
         setting = factor(setting, levels = c("uniform", "nh_weights", "mims_weights", "combo_weights"),
                          labels = c("Uniform", "NHANES", "MIMS", "NHANES&MIMS")),
         var = factor(var, labels = c("Intercept", "Sex Coefficient"))) %>%
  filter(boot_type != "BRR") %>%
  ggplot(aes(x = setting, y = MISE, color = boot_type, fill = boot_type)) +
  facet_grid(. ~ var) +
  geom_boxplot(
    box.colour = "black",
    median.linewidth = 1,
    median.color = "black",
    staplewidth = 0.5, # show staple
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               aes(shape = boot_type),
               size = 2.5,
               color = "black",
               position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c(8, 1)) +
  theme(palette.color.discrete = cols[1:3],
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = "Weighting Scheme", y = "log Mean Integrated Squared Error (MISE)", color = "Estimation Type", fill = "Estimation Type", shape = "Estimation Type") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_sub_axis_x(text = element_text(angle = 30))



p3 = res_df %>%
  mutate(boot_type = factor(boot_type, levels = c("unweighted", "weighted", "BRR"), labels = c("Unweighted", "Weighted", "BRR")),
         setting = factor(setting, levels = c("uniform", "nh_weights", "mims_weights", "combo_weights"),
                          labels = c("Uniform", "NHANES", "MIMS", "NHANES&MIMS")),
         var = factor(var, labels = c("Intercept", "Sex Coefficient"))) %>%
  ggplot(aes(x = setting, y = mean_pw_se, color = boot_type, fill = boot_type)) +
  facet_grid(. ~ var) +
  geom_boxplot(
    box.colour = "black",
    median.linewidth = 1,
    median.color = "black",
    staplewidth = 0.5, # show staple
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               aes(shape = boot_type),
               size = 2.5,
               color = "black",
               position = position_dodge(width = 0.75)) +
  scale_shape_manual(values = c(8, 1, 2)) +
  theme(palette.color.discrete = cols[1:3],
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_axis_x(text = element_text(angle = 30)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = "Weighting Scheme", y = "Mean Pointwise Standard Error", color = "Inference Type", fill = "Inference Type", shape = "Inference Type") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

(p1 + p2) / (p3 + plot_spacer()) + plot_layout(axis_titles = "collect")

p4 = res_df %>%
  mutate(boot_type = factor(boot_type, levels = c("unweighted", "weighted", "BRR"), labels = c("Unweighted", "Weighted", "BRR")),
         setting = factor(setting, levels = c("uniform", "nh_weights", "mims_weights", "combo_weights"),
                          labels = c("Uniform", "NHANES", "MIMS", "NHANES&MIMS")),
         var = factor(var, labels = c("Intercept", "Sex Coefficient"))) %>%
  group_by(var, boot_type, setting) %>%
  summarize(cover_joint = mean(cover_joint)) %>%
  ggplot(aes(x = setting, y = cover_joint, color = boot_type, shape = boot_type)) +
  facet_grid(. ~ var) +
  scale_shape_manual(values = c(8, 19, 17)) +
  geom_jitter(width = .2, size = 3) +
  theme(palette.color.discrete = cols[1:3],
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_sub_axis_x(text = element_text(angle = 30)) +
  labs(x = "Weighting Scheme", y = "Joint Coverage", color = "Inference Type", shape = "Inference Type") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  geom_hline(aes(yintercept = 0.95), color = cols[4], linetype = "dashed", linewidth = 1.2) +
  scale_y_continuous(breaks = c(seq(0, 0.8, 0.2), 0.95, 1))

png(here::here("manuscript", "figures", "empirical_sim_fig.png"), width = 12, height = 10, res = 350,
    units = "in")
(p1 + p2) / (p3 + p4) + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "A")
dev.off()
# svg 1300 x 850

##### ---- table -----
# for results

res_df %>%
  mutate(boot_type = factor(boot_type, levels = c("unweighted", "weighted", "BRR"), labels = c("Unweighted", "Weighted", "BRR")),
         setting = factor(setting, levels = c("uniform", "nh_weights", "mims_weights", "combo_weights"),
                          labels = c("Uniform", "NHANES", "MIMS", "NHANES&MIMS")),
         var = factor(var, labels = c("Intercept", "Sex Coefficient"))) %>%
  group_by(var, boot_type, setting) %>%
  summarize(cover_pw = mean(cover_pw / L),
            cover_joint  = mean(cover_joint)) %>%
  pivot_longer(cols = contains("cover")) %>%
  pivot_wider(names_from = boot_type, values_from = value) %>%
  arrange(setting)
