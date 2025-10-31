library(tidyverse)
library(kableExtra)
library(patchwork)
source(here::here("R", "create_survey_settings_final.R"))

res = read_rds(here::here("results", "simulations", "all_survey_res.rds"))


all_res_prc =
  res %>%
  mutate(boot_type = if_else(boot_type == "RWYB", "Rao-Wu-Yue-Beaumont", boot_type)) %>%
  mutate(cover_pw = cover_pw / len,
         boot_type =  factor(boot_type, levels = c("unweighted", "weighted", "BRR", "Rao-Wu-Yue-Beaumont"),
                             labels = c("Unweighted", "Weighted", "BRR", "RWYB")),
         type = case_when(
           strata_scale == 0 & strata_sigma == 0 ~ "No random effects",
           strata_scale > 0 & strata_sigma > 0 ~ "Strata scaling and noise",
           strata_scale > 0 & strata_sigma == 0 ~ "Strata scaling only",
           strata_scale == 0 & strata_sigma > 0 ~ "Strata noise only"),
         snr_b = if_else(type == "No random effects", NA_real_, snr_b),
         inf_lab = case_when(inf_level == 0 ~ "Uniform sampling",
                             inf_level == 10 ~ "Medium informativeness",
                             inf_level == 15 ~ "High informativeness"),
         inf_lab = factor(inf_lab, levels = c("Uniform sampling", "Medium informativeness", "High informativeness")),
         snr_blab = case_when(snr_b == 0.5 ~ "SNR[b]==0.5",
                              snr_b == 1 ~ "SNR[b]==1",
                              snr_b == 5 ~ "SNR[b]==5"),
         l_lab = paste0("L = ", len),
         l_lab = factor(l_lab, levels = c("L = 50", "L = 100", "L = 1440")),
         n_lab = case_when(In == 100 ~ "I[n]==100",
                           In == 500 ~ "I[n]==500"),

         # factor(paste0("n = ", In)),
         lMISE = log10(MISE))

## start with gaussian only
all_res_g =
  all_res_prc %>%
  filter(family == "gaussian")

all_res_ng =
  all_res_prc %>%
  filter(family != "gaussian") %>%
  mutate(family = factor(family, levels = c("poisson", "binomial"), labels = c("Poisson", "Bernoulli")))



plot_mise = function(df, x_var, facet_x, facet_y, scales = "free_y", custom_labeller, xlab, tl = NULL, lnr = 1) {
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
    ggplot(aes(x = factor(xvar), y = lMISE, color = boot_type, fill = boot_type)) +
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
                 aes(shape = boot_type),
                 size = 2,
                 position = position_dodge(width = 0.75)) +
    theme(palette.color.discrete = cols[1:2],
          strip.text = element_text(size = 14)) +
    theme_sub_legend(position = "bottom",
                     text = element_text(size = 14),
                     title = element_text(size = 15)) +
    theme_sub_axis(text= element_text(size = 12),
                   title = element_text(size = 15)) +
    theme_sub_panel(grid.major.x = element_blank(),
                    grid.minor.x = element_blank()) +
    labs(x = xlab, y = expression(log[10]~"Mean Integrated Squared Error"), color = "Estimation Type", fill = "Estimation Type", shape = "Estimation Type") +
    guides(color = guide_legend(nrow = lnr, byrow = TRUE)) +
    scale_shape_manual(values = c(1, 8))

  if(is.null(facet_x)) {
    p = p +
      facet_grid(. ~ facet_yf, labeller = custom_labeller, scales = scales)
  } else if (is.null(facet_y)) {
    p = p +
      facet_grid(facet_xf ~ ., labeller = custom_labeller, scales = scales)
  } else{
    p = p +
      facet_grid(facet_xf ~ facet_yf, labeller = custom_labeller, scales = scales)
  }
  if(!is.null(tl)) {
    p = p + labs(title = paste0(tl))
  }
  p
}

plot_se = function(df, x_var, facet_x, facet_y, scales = "free_y", custom_labeller, xlab, tl = NULL, lnr = 1) {
  if(is.null(facet_x)) {
    p = df %>%
      rename(xvar := !!sym(x_var), facet_yf := !!sym(facet_y)) %>%
      ggplot(aes(x = factor(xvar), y = mean_pw_se, color = boot_type)) +
      facet_grid(. ~ facet_yf, labeller = custom_labeller, scales = scales) +
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
      theme(palette.color.discrete = cols[1:4],
            strip.text = element_text(size = 12)) +
      theme_sub_legend(position = "bottom",
                       text = element_text(size = 12),
                       title = element_text(size = 13)) +
      theme_sub_axis(text= element_text(size = 10),
                     title = element_text(size = 13)) +
      theme_sub_panel(grid.major.x = element_blank(),
                      grid.minor.x = element_blank()) +
      labs(x = xlab, y = "Mean Pointwise Standard Error", color = "Inference Type") +
      guides(color = guide_legend(nrow = lnr, byrow = TRUE))
  } else if (is.null(facet_y)) {
    p = df %>%
      rename(xvar := !!sym(x_var), facet_xf := !!sym(facet_x)) %>%
      ggplot(aes(x = factor(xvar), y = mean_pw_se, color = boot_type)) +
      facet_grid(facet_xf ~ ., labeller = custom_labeller, scales = scales) +
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
      theme(palette.color.discrete = cols[1:4],
            strip.text = element_text(size = 12)) +
      theme_sub_legend(position = "bottom",
                       text = element_text(size = 12),
                       title = element_text(size = 13)) +
      theme_sub_axis(text= element_text(size = 10),
                     title = element_text(size = 13)) +
      theme_sub_panel(grid.major.x = element_blank(),
                      grid.minor.x = element_blank()) +
      labs(x = xlab, y = "Mean Pointwise Standard Error", color = "Inference Type") +
      guides(color = guide_legend(nrow = lnr, byrow = TRUE))
  } else{
    p = df %>%
      rename(xvar := !!sym(x_var), facet_xf := !!sym(facet_x), facet_yf := !!sym(facet_y)) %>%
      ggplot(aes(x = factor(xvar), y = mean_pw_se, color = boot_type)) +
      facet_grid(facet_xf ~ facet_yf, labeller = custom_labeller, scales = scales) +
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
      theme(palette.color.discrete = cols[1:4],
            strip.text = element_text(size = 12)) +
      theme_sub_legend(position = "bottom",
                       text = element_text(size = 12),
                       title = element_text(size = 13)) +
      theme_sub_axis(text= element_text(size = 10),
                     title = element_text(size = 13)) +
      theme_sub_panel(grid.major.x = element_blank(),
                      grid.minor.x = element_blank()) +
      labs(x = xlab, y = "Mean Pointwise Standard Error",  color = "Inference Type") +
      guides(color = guide_legend(nrow = lnr, byrow = TRUE))
  }
  if(!is.null(tl)) {
    p = p + labs(title = paste0(tl))
  }
  p
}

# ------ MISE plots -------

#####  Gaussian plots
custom_labeller = labeller(
  facet_xf = label_parsed,   # parse snr_blab
  facet_yf  = label_value     # keep inf_lab plain
)
cols = c("#E69F00FF",  "#56B4E9FF", "#009E73FF", "#CC79A7FF")


p1 =
  all_res_g %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "Strata scaling and noise") %>%
  plot_mise(x_var = "snr_eps", facet_y = "inf_lab", facet_x = "snr_blab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]))

p2 = all_res_g %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "Strata scaling and noise") %>%
  plot_mise(x_var = "snr_eps", facet_y = "inf_lab", facet_x = "snr_blab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]))


png(here::here("manuscript", "figures", "mise_gaussian.png"), width = 10, height = 6, units = "in", res = 350)
p1
dev.off()

png(here::here("manuscript", "figures", "mise_gaussian_int.png"), width = 10, height = 6, units = "in", res = 350)
p2
dev.off()
## varying random effects
p1 = all_res_g %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "Strata noise only") %>%
  plot_mise(x_var = "snr_eps", facet_y = "inf_lab", facet_x = "snr_blab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]),
            tl = "Random effects structure: Strata/PSU noise only")



p2 = all_res_g %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "No random effects") %>%
  plot_mise(x_var = "snr_eps", facet_y = "inf_lab", facet_x = NULL,
            scales = "free", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]),
            tl = "Random effects structure: No random effects")


png(here::here("manuscript", "figures", "mise_gaussian_supp.png"), width = 12, height = 8, units = "in", res = 350)
p1 / p2 + plot_layout(heights = c(2, 1),  axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
dev.off()

p1 = all_res_g %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "Strata noise only") %>%
  plot_mise(x_var = "snr_eps", facet_y = "inf_lab", facet_x = "snr_blab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]),
            tl = "Random effects structure: Strata/PSU noise only")



p2 = all_res_g %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "No random effects") %>%
  plot_mise(x_var = "snr_eps", facet_y = "inf_lab", facet_x = NULL,
            scales = "free", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]),
            tl = "Random effects structure: No random effects")


png(here::here("manuscript", "figures", "mise_gaussian_supp_int.png"), width = 12, height = 8, units = "in", res = 350)
p1 / p2 + plot_layout(heights = c(2, 1),  axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
dev.off()

custom_labeller = labeller(
  facet_yf = label_value,   # parse snr_blab
  facet_xf  = label_parsed     # keep inf_lab plain
)

## varying sample size parameters
p1 = all_res_g %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, snr_eps == 1, type == "Strata scaling and noise") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = "Length of functional domain")


png(here::here("manuscript", "figures", "mise_gaussian_n.png"),
    width = 10, height = 6, units = "in", res = 350)
p1
dev.off()

p1 = all_res_g %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, snr_eps == 1, type == "Strata scaling and noise") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = "Length of functional domain")


png(here::here("manuscript", "figures", "mise_gaussian_n_int.png"),
    width = 10, height = 6, units = "in", res = 350)
p1
dev.off()

## varying random effects
p1 = all_res_g %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, snr_eps == 1, type == "Strata noise only") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = "Length of functional domain",
            tl = "Random effects structure: Strata/PSU noise only")

p2 = all_res_g %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_eps == 1, type == "No random effects") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = "Length of functional domain",
            tl = "Random effects structure: No random effects")
  # scale_y_continuous(limits = c(-8.66, -3.23))

png(here::here("manuscript", "figures", "mise_gaussian_n_supp.png"), width = 12, height = 8, units = "in", res = 350)
p1 / p2 + plot_layout(heights = c(2, 1.25),  axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
dev.off()

p1 = all_res_g %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, snr_eps == 1, type == "Strata noise only") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = "Length of functional domain",
            tl = "Random effects structure: Strata/PSU noise only")

p2 = all_res_g %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_eps == 1, type == "No random effects") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = "Length of functional domain",
            tl = "Random effects structure: No random effects")
  # scale_y_continuous(limits = c(-8.66, -3.23))

png(here::here("manuscript", "figures", "mise_gaussian_n_supp_int.png"), width = 12, height = 8, units = "in", res = 350)
p1 / p2 + plot_layout(heights = c(2, 1.25),  axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
dev.off()

#####  Non Gaussian plots

p1 = all_res_ng %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "Strata scaling and noise") %>%
  plot_mise(x_var = "snr_b", facet_y = "inf_lab", facet_x = "family",
            scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[b]))


p2 = all_res_ng %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, family == "Poisson", type == "Strata scaling and noise") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Poisson")

p3 = all_res_ng %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, family == "Bernoulli", type == "Strata scaling and noise") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Bernoulli")


png(here::here("manuscript", "figures", "mise_ng.png"),
    width = 16, height = 12, units = "in", res = 350)
p1 | (p2  / p3) + plot_layout(axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
dev.off()

p1 = all_res_ng %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "Strata scaling and noise") %>%
  plot_mise(x_var = "snr_b", facet_y = "inf_lab", facet_x = "family",
            scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[b]))


p2 = all_res_ng %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, family == "Poisson", type == "Strata scaling and noise") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Poisson")

p3 = all_res_ng %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, family == "Bernoulli", type == "Strata scaling and noise") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "fixed", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Bernoulli")


png(here::here("manuscript", "figures", "mise_ng_int.png"),
    width = 16, height = 12, units = "in", res = 350)
p1 | (p2  / p3) + plot_layout(axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
dev.off()

#### --- supp plots ---
p1 = all_res_ng %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "Strata noise only") %>%
  plot_mise(x_var = "snr_b", facet_y = "inf_lab", facet_x = "family",
            tl = "Random effects structure: Strata/PSU noise only",
            scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[b]))



p3 = all_res_ng %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "No random effects") %>%
  plot_mise(x_var = "snr_b", facet_y = "inf_lab", facet_x = "family",
            scales = "free", custom_labeller = custom_labeller, xlab = "", tl = "Random effects structure: No random effects") +
  theme_sub_axis_x(text = element_blank(),
                   ticks = element_blank())


png(here::here("manuscript", "figures", "mise_ng_supp.png"),
    width = 12, height = 8, units = "in", res = 350)
p1 / p3 + plot_layout(heights = c(1, 1), axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
dev.off()

p1 = all_res_ng %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "Strata noise only") %>%
  plot_mise(x_var = "snr_b", facet_y = "inf_lab", facet_x = "family",
            tl = "Random effects structure: Strata/PSU noise only",
            scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[b]))



p3 = all_res_ng %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), In == 100, len == 50, type == "No random effects") %>%
  plot_mise(x_var = "snr_b", facet_y = "inf_lab", facet_x = "family",
            scales = "free", custom_labeller = custom_labeller, xlab = "", tl = "Random effects structure: No random effects") +
  theme_sub_axis_x(text = element_blank(),
                   ticks = element_blank())


png(here::here("manuscript", "figures", "mise_ng_supp_int.png"),
    width = 12, height = 8, units = "in", res = 350)
p1 / p3 + plot_layout(heights = c(1, 1), axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
dev.off()

custom_labeller = labeller(
  facet_yf = label_value,   # parse snr_blab
  facet_xf  = label_parsed
)

## varying sample size parameters
p1 = all_res_ng %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, family == "Poisson", type == "Strata noise only") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Poisson: Strata/PSU noise only")

p2 = all_res_ng %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, family == "Bernoulli", type == "Strata noise only") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Bernoulli: Strata/PSU noise only")

p3 = all_res_ng %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), family == "Poisson", type == "No random effects") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Poisson: No random effects")

p4 = all_res_ng %>%
  filter(var == "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), family == "Bernoulli", type == "No random effects") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Bernoulli: No random effects")



png(here::here("manuscript", "figures", "mise_n_ng_supp.png"),
    width = 16, height = 8, units = "in", res = 350)
(p1 + p2 + plot_layout(axis_titles = "collect")) / (p3 + p4 + plot_layout(axis_titles = "collect")) + plot_layout(axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
# p1 / p2 + plot_layout(heights = c(1, 1), axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
dev.off()

p1 = all_res_ng %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, family == "Poisson", type == "Strata noise only") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Poisson: Strata/PSU noise only")

p2 = all_res_ng %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), snr_b == 0.5, family == "Bernoulli", type == "Strata noise only") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Bernoulli: Strata/PSU noise only")

p3 = all_res_ng %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), family == "Poisson", type == "No random effects") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Poisson: No random effects")

p4 = all_res_ng %>%
  filter(var != "x") %>%
  filter(boot_type %in% c("Weighted", "Unweighted"), family == "Bernoulli", type == "No random effects") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Bernoulli: No random effects")



png(here::here("manuscript", "figures", "mise_n_ng_supp_int.png"),
    width = 16, height = 8, units = "in", res = 350)
(p1 + p2 + plot_layout(axis_titles = "collect")) / (p3 + p4 + plot_layout(axis_titles = "collect")) + plot_layout(axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
# p1 / p2 + plot_layout(heights = c(1, 1), axis_titles = "collect", guides = "collect") & theme(legend.position = "bottom")
dev.off()
# ------- coverage table --------

# gaussian
colnames_to_row <- function(df) {
  r1 = colnames(df) %>% t() %>% as_tibble()
  colnames(df) = paste0("V", seq_len(ncol(df)))
  out = rbind(r1, df)
  # colnames(out) <- paste0("V", seq_len(ncol(df)))
  out
}



t1 =
  all_res_g %>%
  filter(type == "Strata scaling and noise", snr_b == 0.5, snr_eps == 1, In == 100, len == 50) %>%
  group_by(inf_lab, boot_type, var) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(var, inf_lab), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()



t2 =
  all_res_g %>%
  filter(inf_level == 10, (snr_b == 0.5 | is.na(snr_b)), snr_eps == 1, In == 100, len == 50) %>%
  group_by(type, boot_type, var) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(var, type), values_from = value) %>%
  janitor::clean_names() %>%
  dplyr::select(boot_type, name, contains("no_random_effects"), contains("strata_scaling_only"), everything()) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()


t3 =
  all_res_g %>%
  filter(inf_level == 10, type == "Strata scaling and noise", snr_eps == 1, In == 100, len == 50) %>%
  group_by(snr_blab, boot_type, var) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(var, snr_blab), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()


t4 = all_res_g %>%
  filter(inf_level == 10, type == "Strata scaling and noise", snr_b == 0.5, In == 100, len == 50) %>%
  group_by(snr_eps, boot_type, var) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(var, snr_eps), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()


t5 = all_res_g %>%
  filter(inf_level == 10, type == "Strata scaling and noise", snr_eps == 1, snr_b == 0.5, len == 50) %>%
  group_by(In, boot_type, var) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(var, In), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()


t6 = all_res_g %>%
  filter(inf_level == 10, type == "Strata scaling and noise", snr_eps == 1, snr_b == 0.5, In == 100) %>%
  group_by(l_lab, boot_type, var) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(var, l_lab), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()




library(kableExtra)
tb = reduce(list(t1, t2, t3, t4, t5, t6), bind_rows)
tb %>%
  mutate(V2 = case_when(V2 == "cover_pw" ~ "Pointwise", V2 == "cover_joint" ~ "Joint", TRUE ~ V2),
         V1 = case_when(
           V1 == "Unweighted" ~ "Unweighted",
           V1 == "Weighted" ~ "Weighted",
           .default = V1
         )) %>%
  filter(V2 %in% c("Pointwise", "name")) %>%
  select(-V2) %>%
  kable(format = "latex", booktabs = TRUE)

tb %>%
  mutate(V2 = case_when(V2 == "cover_pw" ~ "Pointwise", V2 == "cover_joint" ~ "Joint", TRUE ~ V2),
         V1 = case_when(
           V1 == "Unweighted" ~ "Unweighted",
           V1 == "Weighted" ~ "Weighted",
           .default = V1
         )) %>%
  filter(V2 %in% c("Joint", "name")) %>%
  select(-V2) %>%
  kable(format = "latex", booktabs = TRUE)
# kable_styling(latex_options = c("hold_position", "striped"))

# non gaussian

t1 =
  all_res_ng %>%
  filter(var == "x") %>%
  filter(type == "Strata scaling and noise", snr_b == 0.5, In == 100, len == 50) %>%
  group_by(inf_lab, boot_type, family) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from =  c(family, inf_lab), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  # ungroup() %>%
  colnames_to_row()


t2 =
  all_res_ng %>%
  filter(var == "x") %>%
  filter(inf_level == 10, (snr_b == 0.5 | is.na(snr_b)), In == 100, len == 50) %>%
  group_by(type, boot_type, family) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(family, type), values_from = value) %>%
  janitor::clean_names() %>%
  # select(-strata_scaling_only) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  ungroup() %>%
  colnames_to_row()

t3 =
  all_res_ng %>%
  filter(var == "x") %>%
  filter(inf_level == 10, type == "Strata scaling and noise", In == 100, len == 50) %>%
  group_by(snr_blab, boot_type, family) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(family, snr_blab), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  ungroup() %>%
  colnames_to_row()


t5 = all_res_ng %>%
  filter(var == "x") %>%
  filter(inf_level == 10, type == "Strata scaling and noise",  snr_b == 0.5, len == 50) %>%
  group_by(In, boot_type, family) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(family, In), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  ungroup() %>%
  colnames_to_row()

t6 = all_res_ng %>%
  filter(var == "x") %>%
  filter(inf_level == 10, type == "Strata scaling and noise",  snr_b == 0.5, In == 100) %>%
  group_by(l_lab, boot_type, family) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(family, l_lab), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  ungroup() %>%
  colnames_to_row()


tb = reduce(list(t1, t2, t3, t5, t6), bind_rows)
tb %>%
  mutate(V2 = case_when(V2 == "cover_pw" ~ "Pointwise", V2 == "cover_joint" ~ "Joint", TRUE ~ V2),
         V1 = case_when(
           V1 == "Unweighted" ~ "Unweighted",
           V1 == "Weighted" ~ "Weighted",
           .default = V1
         )) %>%
  filter(V2 %in% c("Pointwise", "name")) %>%
  select(-V2) %>%
  kable(format = "latex", booktabs = TRUE)

tb %>%
  mutate(V2 = case_when(V2 == "cover_pw" ~ "Pointwise", V2 == "cover_joint" ~ "Joint", TRUE ~ V2),
         V1 = case_when(
           V1 == "Unweighted" ~ "Unweighted",
           V1 == "Weighted" ~ "Weighted",
           .default = V1
         )) %>%
  filter(V2 %in% c("Joint", "name")) %>%
  select(-V2) %>%
  kable(format = "latex", booktabs = TRUE)


### probabilities for fnl intercept

t1 =
  all_res_ng %>%
  filter(var != "x") %>%
  filter(type == "Strata scaling and noise", snr_b == 0.5, In == 100, len == 50) %>%
  group_by(inf_lab, boot_type, family) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from =  c(family, inf_lab), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  # ungroup() %>%
  colnames_to_row()


t2 =
  all_res_ng %>%
  filter(var != "x") %>%
  filter(inf_level == 10, (snr_b == 0.5 | is.na(snr_b)), In == 100, len == 50) %>%
  group_by(type, boot_type, family) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(family, type), values_from = value) %>%
  janitor::clean_names() %>%
  # select(-strata_scaling_only) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  ungroup() %>%
  colnames_to_row()

t3 =
  all_res_ng %>%
  filter(var != "x") %>%
  filter(inf_level == 10, type == "Strata scaling and noise", In == 100, len == 50) %>%
  group_by(snr_blab, boot_type, family) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(family, snr_blab), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  ungroup() %>%
  colnames_to_row()


t5 = all_res_ng %>%
  filter(var != "x") %>%
  filter(inf_level == 10, type == "Strata scaling and noise",  snr_b == 0.5, len == 50) %>%
  group_by(In, boot_type, family) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(family, In), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  ungroup() %>%
  colnames_to_row()

t6 = all_res_ng %>%
  filter(var != "x") %>%
  filter(inf_level == 10, type == "Strata scaling and noise",  snr_b == 0.5, In == 100) %>%
  group_by(l_lab, boot_type, family) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = c(family, l_lab), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  ungroup() %>%
  colnames_to_row()


tb = reduce(list(t1, t2, t3, t5, t6), bind_rows)
tb %>%
  mutate(V2 = case_when(V2 == "cover_pw" ~ "Pointwise", V2 == "cover_joint" ~ "Joint", TRUE ~ V2),
         V1 = case_when(
           V1 == "Unweighted" ~ "Unweighted",
           V1 == "Weighted" ~ "Weighted",
           .default = V1
         )) %>%
  filter(V2 %in% c("Pointwise", "name")) %>%
  select(-V2) %>%
  kable(format = "latex", booktabs = TRUE)

tb %>%
  mutate(V2 = case_when(V2 == "cover_pw" ~ "Pointwise", V2 == "cover_joint" ~ "Joint", TRUE ~ V2),
         V1 = case_when(
           V1 == "Unweighted" ~ "Unweighted",
           V1 == "Weighted" ~ "Weighted",
           .default = V1
         )) %>%
  filter(V2 %in% c("Joint", "name")) %>%
  select(-V2) %>%
  kable(format = "latex", booktabs = TRUE)


### ----- mfpca results ----- ####

mfpca_res = read_rds(here::here("results", "mfpca_sim.rds"))
mfpca_res = mfpca_res %>%
  mutate(type = case_when(
    strata_scale == 0 & strata_sigma == 0 ~ "No random effects",
    strata_scale > 0 & strata_sigma > 0 ~ "Strata scaling and noise",
    strata_scale > 0 & strata_sigma == 0 ~ "Strata scaling only",
    strata_scale == 0 & strata_sigma > 0 ~ "Strata noise only"
  ))

p1 = mfpca_res %>%
  mutate(type = factor(type)) %>%
  ggplot(aes(x = pve_b * 100)) +
  geom_density(aes(fill = type)) +
  facet_grid(type ~ ., scales = "free_y") +
  scale_fill_manual(values = c("#E51932FF", "#19B2FFFF", "#654CFFFF")) +
  labs(x = "Percent variance explained by strata/PSU", y = "Density")  +
  scale_x_continuous(breaks=seq(0,100,10),limits=c(0,100)) +
  theme(strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "none",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank())
# paletteer_d("colorBlindness::PairedColor12Steps")


png(here::here("manuscript", "figures", "pve_mfpca.png"),
    width = 10, height = 8, units = "in", res = 350)
p1
dev.off()

# paletteer_d("colorBlindness::ModifiedSpectralScheme11Steps")
p1 = mfpca_res %>%
  mutate(type = factor(type),
         # type = forcats::fct_relevel(type, "Strata scaling only", after = 1L),
         snr_b = if_else(type == "No random effects", "Not applicable", paste0(snr_b)),
         snr_b = factor(snr_b, levels = c("Not applicable", "0.5", "1", "5"),
                        labels = c("Not applicable", "High", "Medium", "Low"))) %>%
  ggplot(aes(x = pve_b * 100)) +
  geom_density(aes(fill = snr_b, color = snr_b), alpha = 0.5) +
  facet_grid(type ~ ., scales = "free_y") +
  scale_color_manual(values = c("#3FA0FFFF", "#D82632FF", "#F76D5EFF", "#FFAD72FF"), name = "Relative importance of random effects") +
  scale_fill_manual(values = c("#3FA0FFFF", "#D82632FF", "#F76D5EFF", "#FFAD72FF"), name = "Relative importance of random effects") +
  labs(x = "Percent Variance Explained by Strata and PSU", y = "Density")  +
  scale_x_continuous(breaks=seq(0,100,10),limits=c(0,100)) +
  theme(strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_panel(grid.minor.y = element_blank(),
                  grid.minor.x = element_blank())

png(here::here("manuscript", "figures", "pve_mfpca2.png"),
    width = 10, height = 8, units = "in", res = 350)
p1
dev.off()


### find where unweighted or weighted did the worst
lo = mfpca_res %>% arrange(desc(pve_w)) %>% slice(1:10)
lo %>% slice(1) %>% t()
hi = mfpca_res %>% arrange(desc(pve_b)) %>% slice(1:10)
hi %>% slice(1) %>% t()
res_sort = res %>% arrange(desc(pve_w))


diffs_rwyb_x =
  all_res_g %>%
  filter(var == "x") %>%
  select(id, cover_pw, boot_type, fold) %>%
  pivot_wider(names_from = boot_type, values_from = cover_pw) %>%
  janitor::clean_names() %>%
  mutate(diff = rwyb - weighted) %>%
  group_by(fold) %>%
  summarize(across(diff, mean)) %>%
  left_join(settings, by = join_by(fold)) %>%
  left_join(mfpca_res %>% select(fold, pve_b, pve_w), by = join_by(fold)) %>%
  mutate(type = case_when(
    strata_scale == 0 & strata_sigma == 0 ~ "No random effects",
    strata_scale > 0 & strata_sigma > 0 ~ "Strata scaling and noise",
    strata_scale > 0 & strata_sigma == 0 ~ "Strata scaling only",
    strata_scale == 0 & strata_sigma > 0 ~ "Strata noise only"
  )) %>%
  mutate(type = factor(type),
         # type = forcats::fct_relevel(type, "Strata scaling only", after = 1L),
         snr_b = if_else(type == "No random effects", "Not applicable", paste0(snr_b)),
         snr_b = factor(snr_b, levels = c("Not applicable", "0.5", "1", "5"),
                        labels = c("Not applicable", "High", "Medium", "Low")),
         var = "Functional coefficient")
diffs_rwyb_int =
  all_res_g %>%
  filter(var != "x") %>%
  select(id, cover_pw, boot_type, fold) %>%
  pivot_wider(names_from = boot_type, values_from = cover_pw) %>%
  janitor::clean_names() %>%
  mutate(diff = rwyb - weighted) %>%
  group_by(fold) %>%
  summarize(across(diff, mean)) %>%
  left_join(settings, by = join_by(fold)) %>%
  left_join(mfpca_res %>% select(fold, pve_b, pve_w), by = join_by(fold)) %>%
  mutate(type = case_when(
    strata_scale == 0 & strata_sigma == 0 ~ "No random effects",
    strata_scale > 0 & strata_sigma > 0 ~ "Strata scaling and noise",
    strata_scale > 0 & strata_sigma == 0 ~ "Strata scaling only",
    strata_scale == 0 & strata_sigma > 0 ~ "Strata noise only"
  )) %>%
  mutate(type = factor(type),
         # type = forcats::fct_relevel(type, "Strata scaling only", after = 1L),
         snr_b = if_else(type == "No random effects", "Not applicable", paste0(snr_b)),
         snr_b = factor(snr_b, levels = c("Not applicable", "0.5", "1", "5"),
                        labels = c("Not applicable", "High", "Medium", "Low")),
         var = "Functional intercept")


p1 = diffs_rwyb_x %>% bind_rows(diffs_rwyb_int) %>%
  ggplot(aes(x = pve_b * 100, y = diff * 100))  +
  geom_point(aes(color = snr_b, shape = type), size = 3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_shape_discrete(name = "Random effects") +
  facet_wrap(.~var, scales = "free_y") +
  scale_color_manual(values = c("#3FA0FFFF", "#D82632FF", "#F76D5EFF", "#FFAD72FF"), name = "Importance of random effects") +
  geom_smooth(color = "darkgrey") +
  labs(x = "Percent Variance Explained by Strata and PSU", y = "Difference in Coverage: RWYB - Weighted") +
  scale_x_continuous(breaks=seq(0,100,10),limits=c(0,100)) +
  theme(strip.text = element_text(size = 12)) +
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 12),
                   title = element_text(size = 13)) +
  theme_sub_axis(text= element_text(size = 10),
                 title = element_text(size = 13)) +
  theme_sub_panel(grid.minor.y = element_blank(),
                  grid.minor.x = element_blank()) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"),
         shape = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"))
png(here::here("manuscript", "figures", "pve_mfpca_diff.png"),
    width = 10, height = 8, units = "in", res = 350)
p1
dev.off()



