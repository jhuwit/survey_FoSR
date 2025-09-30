library(tidyverse)
library(paletteer)
library(patchwork)
files = list.files(here::here("results", "simulations", "empirical_sim"),
                   pattern = ".rds",
                   full.names = TRUE)
L = 96
# paletteer_d("ggthemes::colorblind")
cols = c("#E69F00FF",  "#56B4E9FF", "#009E73FF", "#CC79A7FF")

weight_types_key = tibble(setting = c("uniform", "nh_weights", "mims_weights", "combo_weights"),
                          fold = c("1", "2", "3", "4"))

fit_types = c('weighted', 'unweighted')
boot_types = c('BRR', 'weighted', 'unweighted')

res_df = map(.x = files,
                .f = read_rds) %>%
  bind_rows(.id = "type") %>%
  left_join(weight_types_key, by = c("type" = "fold"))


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
  ggplot(aes(x = setting, y = cover_pw, color = boot_type)) +
  facet_grid(. ~ var) +
  geom_boxplot(
    # whisker.linetype = "dashed",
    box.colour = "black",
    median.linewidth = 1,
    staplewidth = 0.5, # show staple
    # staple.colour = "grey50",
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 8,
               position = position_dodge(width = 0.75)) +
  theme(palette.color.discrete = cols[1:3],
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  scale_y_continuous(breaks = c(seq(0.2, 0.8, 0.2), 0.95, 1)) +
  labs(x = "Weighting Scheme", y = "Pointwise Coverage", color = "Inference Type") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  geom_hline(aes(yintercept = 0.95), color = cols[4], linetype = "dashed", linewidth = 1.2)


p2 = res_df %>%
  mutate(MISE = log(MISE)) %>%
  mutate(boot_type = factor(boot_type, levels = c("unweighted", "weighted", "BRR"), labels = c("Unweighted", "Weighted", "BRR")),
         setting = factor(setting, levels = c("uniform", "nh_weights", "mims_weights", "combo_weights"),
                          labels = c("Uniform", "NHANES", "MIMS", "NHANES&MIMS")),
         var = factor(var, labels = c("Intercept", "Sex Coefficient"))) %>%
  filter(boot_type != "BRR") %>%
  ggplot(aes(x = setting, y = MISE, color = boot_type)) +
  facet_grid(. ~ var) +
  geom_boxplot(
    # whisker.linetype = "dashed",
    box.colour = "black",
    median.linewidth = 1,
    staplewidth = 0.5, # show staple
    # staple.colour = "grey50",
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 8,
               position = position_dodge(width = 0.75)) +
  theme(palette.color.discrete = cols[1:3],
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = "Weighting Scheme", y = "log Mean Integrated Squared Error (MISE)", color = "Estimation Type") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))



p3 = res_df %>%
  mutate(boot_type = factor(boot_type, levels = c("unweighted", "weighted", "BRR"), labels = c("Unweighted", "Weighted", "BRR")),
         setting = factor(setting, levels = c("uniform", "nh_weights", "mims_weights", "combo_weights"),
                          labels = c("Uniform", "NHANES", "MIMS", "NHANES&MIMS")),
         var = factor(var, labels = c("Intercept", "Sex Coefficient"))) %>%
  ggplot(aes(x = setting, y = mean_pw_se, color = boot_type)) +
  facet_grid(. ~ var) +
  geom_boxplot(
    # whisker.linetype = "dashed",
    box.colour = "black",
    median.linewidth = 1,
    staplewidth = 0.5, # show staple
    # staple.colour = "grey50",
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 8,
               position = position_dodge(width = 0.75)) +
  theme(palette.color.discrete = cols[1:3],
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = "Weighting Scheme", y = "Mean Pointwise Standard Error", color = "Inference Type") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

(p1 + p2) / (p3 + plot_spacer()) + plot_layout(axis_titles = "collect")

p4 = res_df %>%
  mutate(boot_type = factor(boot_type, levels = c("unweighted", "weighted", "BRR"), labels = c("Unweighted", "Weighted", "BRR")),
         setting = factor(setting, levels = c("uniform", "nh_weights", "mims_weights", "combo_weights"),
                          labels = c("Uniform", "NHANES", "MIMS", "NHANES&MIMS")),
         var = factor(var, labels = c("Intercept", "Sex Coefficient"))) %>%
  group_by(var, boot_type, setting) %>%
  summarize(cover_joint = mean(cover_joint)) %>%
  ggplot(aes(x = setting, y = cover_joint, color = boot_type)) +
  facet_grid(. ~ var) +
  geom_jitter(width = .2, size = 3) +
  theme(palette.color.discrete = cols[1:3],
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  labs(x = "Weighting Scheme", y = "Joint Coverage", color = "Inference Type") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  geom_hline(aes(yintercept = 0.95), color = cols[4], linetype = "dashed", linewidth = 1.2) +
  scale_y_continuous(breaks = c(seq(0, 0.8, 0.2), 0.95, 1))

png(here::here("manuscript", "empirical_sim_fig.png"), width = 12, height = 10, res = 350,
    units = "in")
(p1 + p2) / (p3 + p4) + plot_layout(axis_titles = "collect") + plot_annotation(tag_levels = "A")
dev.off()
# svg 1300 x 850

## add table as 4th element
