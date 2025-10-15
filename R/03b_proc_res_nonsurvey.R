# getremote survey_FoSR/results/simulations/non_survey3/*.rds results/simulations/non_survey3/
# getremote survey_FoSR/results/simulations/non_survey2/*.rds results/simulations/non_survey/

library(tidyverse)
library(patchwork)
library(paletteer)
source(here::here("R", "create_nonsurvey_settings.R"))
force = FALSE
# paletteer_d("ggthemes::colorblind")


col1 = "#D55E00FF"
col2 = "#0072B2FF"

if(!file.exists(here::here("results", "simulations", "all_nonsurvey_sim_res.rds")) || force){
  sim_res_fui = list.files(here::here("results", "simulations", "non_survey3"),
                           full.names = TRUE, recursive = TRUE)
  sim_res_fui = sim_res_fui[!grepl("partial", sim_res_fui)]


  sim_res_famm = list.files(here::here("results", "simulations", "non_survey"),
                            full.names = TRUE, recursive = TRUE)


  process_one_file = function(file){
    x = read_rds(file)
    xfold = as.numeric(sub(".*fold\\_(.+).rds.*", "\\1", basename(file)))
    settings = sim_settings %>% filter(fold == xfold)

    x %>%
      mutate(scen = settings$scenario,
             n = settings$n,
             B = settings$num_boots,
             L = settings$L,
             family = settings$family,
             sigma = settings$sigma,
             fold = xfold) %>%
      select(scen, n, L, family, sigma, fold, everything())
  }

  all_res_famm = map_dfr(sim_res_famm,
                         process_one_file) %>%
    filter(type == "FAMM")

  all_res_fui = map_dfr(sim_res_fui,
                        process_one_file)

  all_res =
    bind_rows(all_res_famm, all_res_fui)
  write_rds(all_res, here::here("results", "simulations", "all_nonsurvey_sim_res.rds"), compress = "xz")
} else{
  all_res = read_rds(here::here("results", "simulations", "all_nonsurvey_sim_res.rds"))
}



all_res =
  all_res %>%
  mutate(snr_epslab = case_when(sigma == 0.5 ~ "SNR[epsilon]==0.5",
                                           sigma == 1 ~ "SNR[epsilon]==1",
                                           sigma == 5 ~ "SNR[epsilon]==5"),
         lMISE = log10(MISE))

## create plot fn
plot_mise = function(df, x_var, facet_x, facet_y, scales = "free_y", xlab, tl = NULL, lnr = 1) {
  if(!is.null(facet_x) & !is.null(facet_y)){
    df = df %>%
      rename(xvar := !!sym(x_var), facet_xf := !!sym(facet_x), facet_yf := !!sym(facet_y))
  } else if (is.null(facet_x)) {
    df = df %>%
      rename(xvar := !!sym(x_var), facet_yf := !!sym(facet_y))
  } else if (is.null(facet_y)) {
    df = df %>%
      rename(xvar := !!sym(x_var), facet_xf := !!sym(facet_x))
  }
  p = df %>%
    ggplot(aes(x = factor(xvar), y = lMISE, color = type, fill = type)) +
    geom_boxplot(
      box.colour = "black",
      median.color = "black",
      median.linewidth = .9,
      staplewidth = 0.5, # show staple
      outlier.size = 0.5,
      outlier.alpha = 0.5,
    ) +
    stat_summary(fun = mean,
                 geom = "point",
                 color = "black",
                 aes(shape = type),
                 size = 2,
                 position = position_dodge(width = 0.75)) +
    theme(palette.color.discrete = c(col1, col2),
          strip.text = element_text(size = 14)) +
    theme_sub_legend(position = "bottom",
                     text = element_text(size = 14),
                     title = element_text(size = 15)) +
    theme_sub_axis(text= element_text(size = 12),
                   title = element_text(size = 15)) +
    theme_sub_panel(grid.major.x = element_blank(),
                    grid.minor.x = element_blank()) +
    labs(x = xlab, y = expression(log[10]~"Mean Integrated Squared Error"), color = "Model", fill = "Model", shape = "Model") +
    guides(color = guide_legend(nrow = lnr, byrow = TRUE)) +
    scale_shape_manual(values = c(1, 8))

  if(is.null(facet_x)) {
    p = p +
      facet_grid(. ~ facet_yf,  scales = scales, labeller = custom_labeller)
  } else if (is.null(facet_y)) {
    p = p +
      facet_grid(facet_xf ~ ., scales = scales, labeller = custom_labeller)
  } else{
    p = p +
      facet_grid(facet_xf ~ facet_yf, scales = scales, labeller = custom_labeller)
  }
  if(!is.null(tl)) {
    p = p + labs(title = paste0(tl))
  }
  p
}


###############
# plots
###############

# ------- SNR plots ------

custom_labeller = labeller(
  facet_xf = label_parsed,   # parse snr_blab
  facet_yf  = label_parsed    # keep inf_lab plain
)


p1 = all_res %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(B == 500, family == "gaussian", var == "x", L == 100, n == 500) %>%
  plot_mise(x_var = "scen", facet_y = "snr_epslab", facet_x = NULL,
            scales = "fixed", xlab = "Scenario", tl = "Gaussian response", lnr = 1)

p2 = all_res %>%
  filter(family != "gaussian", var == "x", L == 100, n == 500) %>%
  mutate(scen = paste0("Scenario ", scen),
         family = if_else(family == "binomial", "Bernoulli", "Poisson")) %>%
  plot_mise(x_var = "scen", facet_y = "family", facet_x = NULL,
            scales = "free_y", xlab = "Scenario", tl = "Non-Gaussian response", lnr = 1)


# SNR final plot
png(here::here("manuscript", "figures", "snr_nosurvey.png"),
    width = 10, height = 6, res = 300, units = "in")
p1 / p2 + plot_layout(axis_titles = "collect", guides = "collect")  & theme(legend.position = "bottom")
dev.off()

# sample size final plot
custom_labeller = labeller(
  facet_xf = label_value,   # parse snr_blab
  facet_yf  = label_value    # keep inf_lab plain
)

p1 =
  all_res %>%
  mutate(n = paste0("n = ", n),
         n = factor(n, levels = c("n = 500", "n = 1000", "n = 5000", "n = 10000"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "gaussian", var == "x", sigma == 1, L == 100) %>%
  plot_mise(x_var = "scen", facet_y = "n", facet_x = NULL,
            scales = "free_y", xlab = "Scenario", lnr = 1)

p1a = all_res %>%
  mutate(L = paste0("L = ", L),
         L = factor(L, levels = c("L = 50", "L = 100", "L = 500"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "gaussian", var == "x", sigma == 1, n == 500) %>%
  plot_mise(x_var = "scen", facet_y = "L", facet_x = NULL,
            scales = "free_y", xlab = "Scenario", lnr = 1)

png(here::here("manuscript", "figures", "samplesize_nosurvey.png"),
    width = 10, height = 6, res = 300, units = "in")
p1 / p1a + plot_layout(guides = "collect", axis_titles =  "collect") & theme(legend.position = "bottom")
dev.off()

# -------- sample size plots -----


## for binomial
p1 = all_res %>%
  mutate(n = paste0("n = ", n),
         n = factor(n, levels = c("n = 500", "n = 1000", "n = 5000", "n = 10000"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "binomial", var == "x",  L == 100) %>%
  plot_mise(x_var = "scen", facet_y = "n", facet_x = NULL,
            scales = "free_y", xlab = "Scenario", lnr = 1)

p1a = all_res %>%
  mutate(L = paste0("L = ", L),
         L = factor(L, levels = c("L = 50", "L = 100", "L = 500"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "binomial", var == "x", n == 500) %>%
  plot_mise(x_var = "scen", facet_y = "L", facet_x = NULL,
            scales = "free_y", xlab = "Scenario", lnr = 1)

png(here::here("manuscript", "figures", "samplesize_nosurvey_bin.png"),
    width = 10, height = 6, res = 300, units = "in")
p1 / p1a + plot_layout(guides = "collect", axis_titles =  "collect") + plot_annotation(title = "Bernoulli response")& theme(legend.position = "bottom")
dev.off()

## for binomial
p1 = all_res %>%
  mutate(n = paste0("n = ", n),
         n = factor(n, levels = c("n = 500", "n = 1000", "n = 5000", "n = 10000"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "poisson", var == "x",  L == 100) %>%
  plot_mise(x_var = "scen", facet_y = "n", facet_x = NULL,
            scales = "free_y", xlab = "Scenario", lnr = 1)

p1a = all_res %>%
  mutate(L = paste0("L = ", L),
         L = factor(L, levels = c("L = 50", "L = 100", "L = 500"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "poisson", var == "x", n == 500) %>%
  plot_mise(x_var = "scen", facet_y = "L", facet_x = NULL,
            scales = "free_y", xlab = "Scenario", lnr = 1)

png(here::here("manuscript", "figures", "samplesize_nosurvey_poisson.png"),
    width = 10, height = 6, res = 300, units = "in")
p1 / p1a + plot_layout(guides = "collect", axis_titles =  "collect") + plot_annotation(title = "Poisson response")& theme(legend.position = "bottom")
dev.off()



# -------- standard error plots ------

p1 = all_res %>%
  mutate(scen = paste0("Gaussian scenario ", scen)) %>%
  filter(B == 500, family == "gaussian", var == "x", L == 100, n == 500) %>%
  ggplot(aes(x = factor(scen), y = mean_pw_se, color = type)) +
  geom_boxplot(
    box.colour = "black",
    median.linewidth = 1,
    staplewidth = 0.5, # show staple
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 8,
               position = position_dodge(width = 0.75)) +
  theme(palette.color.discrete = c(col1, col2),
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  facet_wrap(~snr_epslab, labeller = label_parsed, scales = "free_y") +
  labs(x = "Scenario",
       y = "Mean Pointwise Standard Error",
       color = "Model")
# x = "\u03c3 for gaussian response",
# y = bquote(Mean~Integrated~Squared~Error%.%10^4))

p2 = all_res %>%
  filter(family != "gaussian", var == "x", L == 100, n == 500) %>%
  mutate(scen = paste0("Scenario ", scen),
         family = if_else(family == "binomial", "Binomial", "Poisson")) %>%
  ggplot(aes(x = factor(scen), y = mean_pw_se, color = type)) +
  facet_wrap(~family, scales = "free_y") +
  geom_boxplot(
    box.colour = "black",
    median.linewidth = 1,
    staplewidth = 0.5, # show staple
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 8,
               position = position_dodge(width = 0.75)) +
  theme(palette.color.discrete = c(col1, col2),
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  labs(x = "Scenario",
       y = "Mean Pointwise Standard Error",
       color = "Model")
# ggpubr::stat_compare_means(aes(group = type), label = "p.format", method = "t.test", label.y.npc = 0.95, size = 5)


# SNR final plot
png(here::here("manuscript", "figures", "snr_nosurvey_se.png"),
    width = 8, height = 6, res = 300, units = "in")
p1 / p2 + plot_layout(axis_titles = "collect", guides = "collect")  & theme(legend.position = "bottom")
dev.off()
# sample size final plot
p1 = all_res %>%
  mutate(n = paste0("n = ", n),
         n = factor(n, levels = c("n = 500", "n = 1000", "n = 5000", "n = 10000"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "gaussian", var == "x", sigma == 1, L == 100) %>%
  ggplot(aes(x = factor(scen), y = mean_pw_se, color = type)) +
  facet_wrap(.~ n, scales = "free_y", nrow = 1) +
  geom_boxplot(
    box.colour = "black",
    median.linewidth = 1,
    staplewidth = 0.5, # show staple
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 8,
               position = position_dodge(width = 0.75)) +
  theme(palette.color.discrete = c(col1, col2),
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  labs(x = "Scenario",
       y = "Mean Pointwise Standard Error",
       color = "Model")
p1a = all_res %>%
  mutate(L = paste0("L = ", L),
         L = factor(L, levels = c("L = 50", "L = 100", "L = 500"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "gaussian", var == "x", sigma == 1, n == 500) %>%
  ggplot(aes(x = factor(scen), y = mean_pw_se, color = type)) +
  facet_wrap(. ~ L, scales = "free_y", nrow = 1) +
  geom_boxplot(
    box.colour = "black",
    median.linewidth = 1,
    staplewidth = 0.5, # show staple
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 8,
               position = position_dodge(width = 0.75)) +
  theme(palette.color.discrete = c(col1, col2),
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  labs(x = "Scenario",
       y = "Mean Pointwise Standard Error",
       color = "Model")



png(here::here("manuscript", "figures", "samplesize_nosurvey_se.png"),
    width = 8, height = 6, res = 300, units = "in")
p1 / p1a + plot_layout(guides = "collect", axis_titles =  "collect") & theme(legend.position = "bottom")
dev.off()

## for binomial
p1 = all_res %>%
  mutate(n = paste0("n = ", n),
         n = factor(n, levels = c("n = 500", "n = 1000", "n = 5000", "n = 10000"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "binomial", var == "x",  L == 100) %>%
  ggplot(aes(x = factor(scen), y = mean_pw_se, color = type)) +
  facet_wrap(.~ n, scales = "free_y", nrow = 1) +
  geom_boxplot(
    box.colour = "black",
    median.linewidth = 1,
    staplewidth = 0.5, # show staple
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 8,
               position = position_dodge(width = 0.75)) +
  theme(palette.color.discrete = c(col1, col2),
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  labs(x = "Scenario",
       y = "Mean Pointwise Standard Error",
       color = "Model")
p1a = all_res %>%
  mutate(L = paste0("L = ", L),
         L = factor(L, levels = c("L = 50", "L = 100", "L = 500"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "binomial", var == "x", n == 500) %>%
  ggplot(aes(x = factor(scen), y = mean_pw_se, color = type)) +
  facet_wrap(.~L, scales = "free_y", nrow = 1) +
  geom_boxplot(
    box.colour = "black",
    median.linewidth = 1,
    staplewidth = 0.5, # show staple
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 8,
               position = position_dodge(width = 0.75)) +
  theme(palette.color.discrete = c(col1, col2),
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  labs(x = "Scenario",
       y = "Mean Pointwise Standard Error",
       color = "Model")



png(here::here("manuscript", "figures", "samplesize_nosurvey_bin_se.png"),
    width = 8, height = 6, res = 300, units = "in")
p1 / p1a + plot_layout(guides = "collect", axis_titles =  "collect") & theme(legend.position = "bottom")
dev.off()

## for poisson
p1 = all_res %>%
  mutate(n = paste0("n = ", n),
         n = factor(n, levels = c("n = 500", "n = 1000", "n = 5000", "n = 10000"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "poisson", var == "x", L == 100) %>%
  ggplot(aes(x = factor(scen), y = mean_pw_se, color = type)) +
  facet_wrap(.~ n, scales = "free_y", nrow = 1) +
  geom_boxplot(
    box.colour = "black",
    median.linewidth = 1,
    staplewidth = 0.5, # show staple
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 8,
               position = position_dodge(width = 0.75)) +
  theme(palette.color.discrete = c(col1, col2),
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  labs(x = "Scenario",
       y = "Mean Pointwise Standard Error",
       color = "Model")
p1a = all_res %>%
  mutate(L = paste0("L = ", L),
         L = factor(L, levels = c("L = 50", "L = 100", "L = 500"))) %>%
  mutate(scen = paste0("Scenario ", scen)) %>%
  filter(family == "poisson", var == "x",  n == 500) %>%
  ggplot(aes(x = factor(scen), y = mean_pw_se, color = type)) +
  facet_wrap(.~ L, scales = "free_y", nrow = 1) +
  geom_boxplot(
    box.colour = "black",
    median.linewidth = 1,
    staplewidth = 0.5, # show staple
    outlier.size = 0.5,
    outlier.alpha = 0.5
  ) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 8,
               position = position_dodge(width = 0.75)) +
  theme(palette.color.discrete = c(col1, col2),
        strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  labs(x = "Scenario",
       y = "Mean Pointwise Standard Error",
       color = "Model")



png(here::here("manuscript", "figures", "samplesize_nosurvey_poisson_se.png"),
    width = 8, height = 6, res = 300, units = "in")
p1 / p1a + plot_layout(guides = "collect", axis_titles =  "collect") & theme(legend.position = "bottom")
dev.off()





# ==============
# tables
# ==============
# ------ coverage -----
# baseline setting: n = 500, L = 100, sigma = 1, scenario = 1, response = gaussian

colnames_to_row <- function(df) {
  out <- rbind(colnames(df), df)
  colnames(out) <- paste0("V", seq_len(ncol(df)))
  out
}

### scenario 1
t1 = all_res %>%
  filter(B == 500, family == "gaussian", sigma == 1, scen == 1, L == 100, var == "x") %>%
  group_by(type, n) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = n, values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  colnames_to_row()
  # kableExtra::kable(format = "latex")


t2 = all_res %>%
  filter(B == 500, family == "gaussian", sigma == 1, scen == 1, n == 500, var == "x") %>%
  group_by(type, L) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = L, values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  colnames_to_row()
  # kableExtra::kable(format = "latex")



t3 = all_res %>%
  filter(B == 500, scen == 1, n == 500, L == 100, var == "x") %>%
  filter(!(sigma != 1 & family == "gaussian")) %>%
  group_by(type, family) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = family, values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  select(type, metric, gaussian, everything()) %>%
  colnames_to_row()

  # kableExtra::kable(format = "latex")

t4 = all_res %>%
  filter(B == 500,scen == 1, n == 500, L == 100, family == "gaussian", var == "x") %>%
  group_by(type, sigma) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = sigma, values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  colnames_to_row()


tb = reduce(list(t1, t2, t4, t3), bind_rows) %>%
  mutate(V2 = case_when(grepl("pw", V2) ~ "Pointwise coverage",
                        grepl("joint", V2) ~ "Joint Coverage",
                        .default = V2))
tb %>%
  kable(format = "latex", booktabs = TRUE)

### scenario 2
t1 = all_res %>%
  filter(B == 500, family == "gaussian", sigma == 1, scen == 2, L == 100, var == "x") %>%
  group_by(type, n) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = n, values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  colnames_to_row()
# kableExtra::kable(format = "latex")


t2 = all_res %>%
  filter(B == 500, family == "gaussian", sigma == 1, scen == 2, n == 500, var == "x") %>%
  group_by(type, L) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = L, values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  colnames_to_row()
# kableExtra::kable(format = "latex")



t3 = all_res %>%
  filter(B == 500, scen == 2, n == 500, L == 100, var == "x") %>%
  filter(!(sigma != 1 & family == "gaussian")) %>%
  group_by(type, family) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = family, values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  select(type, metric, gaussian, everything()) %>%
  colnames_to_row()

# kableExtra::kable(format = "latex")

t4 = all_res %>%
  filter(B == 500,scen == 2, n == 500, L == 100, family == "gaussian", var == "x") %>%
  group_by(type, sigma) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = sigma, values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  colnames_to_row()


tb = reduce(list(t1, t2, t4, t3), bind_rows) %>%
  mutate(V2 = case_when(grepl("pw", V2) ~ "Pointwise coverage",
                        grepl("joint", V2) ~ "Joint Coverage",
                        .default = V2))
tb %>%
  kable(format = "latex", booktabs = TRUE)


## scenarios together
t1 = all_res %>%
  filter(B == 500, family == "gaussian", sigma == 1, L == 100, var == "x") %>%
  group_by(type, n, scen) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = c(scen, n), values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  colnames_to_row()
# kableExtra::kable(format = "latex")


t2 = all_res %>%
  filter(B == 500, family == "gaussian", sigma == 1, n == 500, var == "x") %>%
  group_by(type, L, scen) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = c(scen, L), values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  colnames_to_row()
# kableExtra::kable(format = "latex")



t3 = all_res %>%
  filter(B == 500, n == 500, L == 100, var == "x") %>%
  filter(!(sigma != 1 & family == "gaussian")) %>%
  group_by(type, family, scen) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = c(scen, family), values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  select(type, metric, contains("gaussian"), everything()) %>%
  colnames_to_row()

# kableExtra::kable(format = "latex")

t4 = all_res %>%
  filter(B == 500, n == 500, L == 100, family == "gaussian", var == "x") %>%
  group_by(type, sigma, scen) %>%
  summarize(mean_pw = mean(cover_pw),
            mean_joint = mean(cover_joint),
            .groups = "drop") %>%
  filter(type != "Bayes")  %>%
  pivot_longer(cols = starts_with("mean"),
               names_to = "metric",
               values_to = "value") %>%
  drop_na() %>%
  pivot_wider(names_from = c(scen, sigma), values_from = value) %>%
  mutate(across(where(is.numeric), ~sprintf("%.2f", signif(.x, 2)))) %>%
  colnames_to_row()


tb = reduce(list(t1, t2, t4, t3), bind_rows) %>%
  mutate(V2 = case_when(grepl("pw", V2) ~ "Pointwise coverage",
                        grepl("joint", V2) ~ "Joint Coverage",
                        .default = V2))
tb %>%
  kable(format = "latex", booktabs = TRUE)

# ------ timing results -----
## table
all_res %>%
  filter(var == "x") %>%
  group_by(scen, n, L, family, sigma, type) %>%
  summarize(across(time, c(median = median, mean = mean, min=min, max=max),
                   .names = "{col}_{fn}")) %>%
  mutate(across(starts_with("time"),  ~signif(.x, 2))) %>%
  mutate(mean_range = paste0(time_mean, " (", time_min, ",", time_max, ")")) %>%
  mutate(mean_range = paste0(round(time_mean,1), " (", round(time_min,1), "-", round(time_max,1), ")")) %>%
  select(-starts_with("time")) %>%
  pivot_wider(names_from = type, values_from = mean_range)


all_res %>%
  filter(sigma == 1, family == "gaussian", scen == 1, n == 500, var == "x") %>%
  ggplot(aes(x = factor(L), y = time, color = type)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = "Length of Functional Domain")

all_res %>%
  filter(sigma == 1, family == "gaussian", scen == 1, L == 100, var == "x") %>%
  ggplot(aes(x = factor(n), y = time, color = type)) +
  geom_boxplot() +
  scale_y_log10()

