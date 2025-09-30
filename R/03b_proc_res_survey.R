# get results from cluster
# getremote survey_FoSR/results/simulations/survey_sim2/* results/simulations/survey_sim2/

library(tidyverse)
library(kableExtra)
library(patchwork)
source(here::here("R", "create_survey_settings.R"))
cols = c("#E69F00FF",  "#56B4E9FF", "#009E73FF", "#CC79A7FF")


force = FALSE
# save the raw simulation results if they don't already exist
if (!file.exists(here::here("results", "simulations", "all_survey_sim_res.rds")) || force ){
  files = list.files(here::here("results", "simulations", "survey_sim2"), full.names = TRUE)
  files = files[!grepl("partial", files)]

  all_res =
    map_dfr(.x = files,
            .f = \(x) read_rds(x) %>%
              mutate(fold = sub(".*fold\\_(\\d+)\\.rds.*", "\\1", x),
                     fold = as.integer(fold))) %>%
    select(id, MISE, starts_with("mean"), starts_with("cover"), var, boot_type, n, n_boot, fold) %>%
    right_join(settings_new, by = join_by(fold))

  # check to make sure we have 200 iterations * 4 methods * 2 variables = 1600 iterations
  wrong_nrow =
    all_res %>%
    group_by(fold) %>%
    count() %>%
    filter(n != 1600)

  # make sure we have all 420 folds
  test_unique = all.equal(unique(all_res$fold), 1:420)

  if(nrow(wrong_nrow) > 0) stop("Some folds do not have 1600 iterations")
  if(!test_unique) stop("Not all folds are present")
  # if we made it here, everything is good
  write_rds(all_res, here::here("results", "simulations", "all_survey_sim_res.rds"))
} else {
  all_res = read_rds(here::here("results", "simulations", "all_survey_sim_res.rds"))
}

##### now process for easier plotting
all_res_prc =
  all_res %>%
  mutate(cover_pw = cover_pw / len,
         boot_type =  factor(boot_type, levels = c("unweighted", "weighted", "BRR", "Rao-Wu-Yue-Beaumont")),
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
         l_lab = factor(l_lab, levels = c("L = 50", "L = 100")),
         n_lab = factor(paste0("n = ", In)),
         lMISE = log(MISE))

## start with gaussian only
all_res_g =
  all_res_prc %>%
  filter(family == "gaussian" & var == "x")

all_res_ng =
  all_res_prc %>%
  filter(family != "gaussian" & var == "x")

plot_mise = function(df, x_var, facet_x, facet_y, scales = "free_y", custom_labeller, xlab, tl = NULL, lnr = 1) {
  if(is.null(facet_x)) {
    p = df %>%
      rename(xvar := !!sym(x_var), facet_yf := !!sym(facet_y)) %>%
      ggplot(aes(x = factor(xvar), y = lMISE, color = boot_type)) +
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
      theme(palette.color.discrete = cols[1:2],
            strip.text = element_text(size = 12)) +
      theme_sub_legend(position = "bottom",
                       text = element_text(size = 12),
                       title = element_text(size = 13)) +
      theme_sub_axis(text= element_text(size = 10),
                     title = element_text(size = 13)) +
      theme_sub_panel(grid.major.x = element_blank(),
                      grid.minor.x = element_blank()) +
      labs(x = xlab, y = "log MISE", color = "Estimation Type") +
      guides(color = guide_legend(nrow = lnr, byrow = TRUE))
  } else if (is.null(facet_y)) {
    p = df %>%
      rename(xvar := !!sym(x_var), facet_xf := !!sym(facet_x)) %>%
      ggplot(aes(x = factor(xvar), y = lMISE, color = boot_type)) +
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
      theme(palette.color.discrete = cols[1:2],
            strip.text = element_text(size = 12)) +
      theme_sub_legend(position = "bottom",
                       text = element_text(size = 12),
                       title = element_text(size = 13)) +
      theme_sub_axis(text= element_text(size = 10),
                     title = element_text(size = 13)) +
      theme_sub_panel(grid.major.x = element_blank(),
                      grid.minor.x = element_blank()) +
      labs(x = xlab, y = "log MISE", color = "Estimation Type") +
      guides(color = guide_legend(nrow = lnr, byrow = TRUE))
  } else{
    p = df %>%
      rename(xvar := !!sym(x_var), facet_xf := !!sym(facet_x), facet_yf := !!sym(facet_y)) %>%
      ggplot(aes(x = factor(xvar), y = lMISE, color = boot_type)) +
      facet_grid(facet_xf ~ facet_yf, labeller = custom_labeller, scales = scales) +
      geom_boxplot(
        box.colour = "black",
        median.linewidth = 1,
        staplewidth = 0.5, # show staple
        outlier.size = 0.5,
        outlier.alpha = 0.5,
      ) +
      stat_summary(fun = mean,
                   geom = "point",
                   shape = 8,
                   position = position_dodge(width = 0.75)) +
      theme(palette.color.discrete = cols[1:2],
            strip.text = element_text(size = 12)) +
      theme_sub_legend(position = "bottom",
                       text = element_text(size = 12),
                       title = element_text(size = 13)) +
      theme_sub_axis(text= element_text(size = 10),
                     title = element_text(size = 13)) +
      theme_sub_panel(grid.major.x = element_blank(),
                      grid.minor.x = element_blank()) +
      labs(x = xlab, y = "log MISE", color = "Estimation Type") +
      guides(color = guide_legend(nrow = lnr, byrow = TRUE))
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
  facet_yf = label_parsed,   # parse snr_blab
  facet_xf  = label_value     # keep inf_lab plain
)

## varying sampling level and SNR parameters
p1 = all_res_g %>%
  filter(boot_type %in% c("weighted", "unweighted"), In == 100, len == 50, type == "Strata scaling and noise") %>%
  plot_mise(x_var = "snr_eps", facet_x = "inf_lab", facet_y = "snr_blab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]))

p2 = all_res_g %>%
  filter(In == 100, len == 50, type == "Strata scaling and noise") %>%
  plot_se(x_var = "snr_eps", facet_x = "inf_lab", facet_y = "snr_blab", lnr = 2,
            scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]))

png(here::here("manuscript", "figures", "mise_se_sampling.png"),
    width = 10, height = 8, units = "in", res = 350)
(p1 | p2) & plot_annotation(tag_levels = "A")
dev.off()

## varying random effects
p1 = all_res_g %>%
  filter(boot_type %in% c("weighted", "unweighted"), In == 100, len == 50, type == "Strata scaling only") %>%
  plot_mise(x_var = "snr_eps", facet_x = "inf_lab", facet_y = "snr_blab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]),
            tl = "MISE: Strata scaling only")

p2 = all_res_g %>%
  filter(boot_type %in% c("weighted", "unweighted"), In == 100, len == 50, type == "No random effects") %>%
  plot_mise(x_var = "snr_eps", facet_x = "inf_lab", facet_y = NULL,
            scales = "free", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]),
            tl = "MISE: No random effects")

p3 =  all_res_g %>%
  filter(In == 100, len == 50, type == "Strata scaling only") %>%
  plot_se(x_var = "snr_eps", facet_x = "inf_lab", facet_y = "snr_blab", lnr = 2,
          scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]),
          tl = "Standard errors: Strata scaling only")

p4 = all_res_g %>%
  filter(In == 100, len == 50, type == "No random effects") %>%
  plot_se(x_var = "snr_eps", facet_x = "inf_lab", facet_y =  NULL, lnr = 2,
          scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[epsilon]),
          tl = "Standard errors: No random effects")

png(here::here("manuscript", "figures", "mise_se_sampling_supp.png"), width = 12, height = 16, units = "in", res = 350)
(p1 | p3) / (p2 | p4) & plot_annotation(tag_levels = "A")
dev.off()

custom_labeller = labeller(
  facet_yf = label_value,   # parse snr_blab
  facet_xf  = label_value     # keep inf_lab plain
)

## varying sample size parameters
p1 = all_res_g %>%
  filter(boot_type %in% c("weighted", "unweighted"), snr_b == 0.5, snr_eps == 1, type == "Strata scaling and noise") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain")

p2 = all_res_g %>%
  filter(snr_b == 0.5, snr_eps == 1, type == "Strata scaling and noise") %>%
  plot_se(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab", lnr = 2,
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain")

png(here::here("manuscript", "figures", "mise_se_sampsize.png"),
    width = 14, height = 6, units = "in", res = 350)
(p1 | p2) & plot_annotation(tag_levels = "A")
dev.off()

#####  Non Gaussian plots

p1 = all_res_ng %>%
  filter(boot_type %in% c("weighted", "unweighted"), In == 100, len == 50, type == "Strata scaling and noise") %>%
  plot_mise(x_var = "snr_b", facet_x = "inf_lab", facet_y = "family",
            scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[b]))

p2 = all_res_ng %>%
  filter(In == 100, len == 50, type == "Strata scaling and noise") %>%
  plot_se(x_var = "snr_b", facet_x = "inf_lab", facet_y = "family", lnr = 2,
            scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[b]))

png(here::here("manuscript", "figures", "mise_se_sampling_ng.png"),
    width = 10, height = 8, units = "in", res = 350)
(p1 | p2) & plot_annotation(tag_levels = "A")
dev.off()

p1 = all_res_ng %>%
  filter(boot_type %in% c("weighted", "unweighted"), In == 100, len == 50, type == "Strata scaling only") %>%
  plot_mise(x_var = "snr_b", facet_x = "inf_lab", facet_y = "family",
            scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[b]))

p2 = all_res_ng %>%
  filter(In == 100, len == 50, type == "Strata scaling only") %>%
  plot_se(x_var = "snr_b", facet_x = "inf_lab", facet_y = "family", lnr = 2,
          scales = "free_y", custom_labeller = custom_labeller, xlab = expression(SNR[b]))

png(here::here("manuscript", "figures", "mise_se_sampling_ng_supp1.png"),
    width = 10, height = 8, units = "in", res = 350)
(p1 | p2) & plot_annotation(tag_levels = "A", title = "Strata scaling only")
dev.off()

p3 = all_res_ng %>%
  filter(boot_type %in% c("weighted", "unweighted"), In == 100, len == 50, type == "No random effects") %>%
  plot_mise(x_var = "family", facet_x = "inf_lab", facet_y = NULL,
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Family")

p4 = all_res_ng %>%
  filter(In == 100, len == 50, type == "No random effects") %>%
  plot_se(x_var = "family", facet_x = "inf_lab", facet_y = NULL, lnr = 2,
          scales = "free_y", custom_labeller = custom_labeller, xlab = "Family")

png(here::here("manuscript", "figures", "mise_se_sampling_ng_supp2.png"),
    width = 10, height = 8, units = "in", res = 350)
(p3 | p4) & plot_annotation(tag_levels = "A", title = "No random effects")
dev.off()


## varying sample size parameters
p1 = all_res_ng %>%
  filter(boot_type %in% c("weighted", "unweighted"), snr_b == 0.5, family == "poisson", type == "Strata scaling and noise") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Poisson")

p2 = all_res_ng %>%
  filter(snr_b == 0.5, family == "poisson", type == "Strata scaling and noise") %>%
  plot_se(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab", lnr = 2,
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Poisson")


png(here::here("manuscript", "figures", "mise_se_sampsize_poisson.png"),
    width = 10, height = 5, units = "in", res = 350)
(p1 | p2) & plot_annotation(tag_levels = "A")
dev.off()

## varying sample size parameters
p1 = all_res_ng %>%
  filter(boot_type %in% c("weighted", "unweighted"), snr_b == 0.5, family == "binomial", type == "Strata scaling and noise") %>%
  plot_mise(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab",
            scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Binomial")

p2 = all_res_ng %>%
  filter(snr_b == 0.5, family == "binomial", type == "Strata scaling and noise") %>%
  plot_se(x_var = "l_lab", facet_x = "n_lab", facet_y = "inf_lab", lnr = 2,
          scales = "free_y", custom_labeller = custom_labeller, xlab = "Length of functional domain", tl = "Binomial")


png(here::here("manuscript", "figures", "mise_se_sampsize_bin.png"),
    width = 10, height = 5, units = "in", res = 350)
(p1 | p2) & plot_annotation(tag_levels = "A")
dev.off()


# ------- coverage table --------

# gaussian
colnames_to_row <- function(df) {
  out <- rbind(colnames(df), df)
  colnames(out) <- paste0("V", seq_len(ncol(df)))
  out
}

t1 =
  all_res_g %>%
  filter(type == "Strata scaling and noise", snr_b == 0.5, snr_eps == 1, In == 100, len == 50) %>%
  group_by(inf_lab, boot_type) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = inf_lab, values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()

t2 =
  all_res_g %>%
  filter(inf_level == 10, (snr_b == 0.5 | is.na(snr_b)), snr_eps == 1, In == 100, len == 50) %>%
  group_by(type, boot_type) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  janitor::clean_names() %>%
  dplyr::select(boot_type, name, no_random_effects, strata_scaling_only, everything()) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()

t3 =
  all_res_g %>%
  filter(inf_level == 10, type == "Strata scaling and noise", snr_eps == 1, In == 100, len == 50) %>%
  group_by(snr_blab, boot_type) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = snr_blab, values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()

t4 = all_res_g %>%
  filter(inf_level == 10, type == "Strata scaling and noise", snr_b == 0.5, In == 100, len == 50) %>%
  group_by(snr_eps, boot_type) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = snr_eps, values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()

t5 = all_res_g %>%
  filter(inf_level == 10, type == "Strata scaling and noise", snr_eps == 1, snr_b == 0.5, len == 50) %>%
  group_by(In, boot_type) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = In, values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()

t6 = all_res_g %>%
  filter(inf_level == 10, type == "Strata scaling and noise", snr_eps == 1, snr_b == 0.5, In == 100) %>%
  group_by(l_lab, boot_type) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from = l_lab, values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row()


tb = reduce(list(t1, t2, t3, t4, t5, t6), bind_rows)
tb %>%
  mutate(V2 = case_when(V2 == "cover_pw" ~ "Pointwise", V2 == "cover_joint" ~ "Joint", TRUE ~ V2),
         V1 = case_when(
           V1 == "unweighted" ~ "Unweighted",
           V1 == "weighted" ~ "Weighted",
           .default = V1
         )) %>%
  filter(V2 %in% c("Pointwise", "name")) %>%
  select(-V2) %>%
  kable(format = "latex", booktabs = TRUE)

tb %>%
  mutate(V2 = case_when(V2 == "cover_pw" ~ "Pointwise", V2 == "cover_joint" ~ "Joint", TRUE ~ V2),
         V1 = case_when(
           V1 == "unweighted" ~ "Unweighted",
           V1 == "weighted" ~ "Weighted",
           .default = V1
         )) %>%
  filter(V2 %in% c("Joint", "name")) %>%
  select(-V2) %>%
  kable(format = "latex", booktabs = TRUE)
# kable_styling(latex_options = c("hold_position", "striped"))

# non gaussian

t1 =
  all_res_ng %>%
  filter(type == "Strata scaling and noise", snr_eps == 1, snr_b == 0.5, In == 100, len == 50) %>%
  group_by(inf_lab, boot_type, family) %>%
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) %>%
  pivot_longer(cols = starts_with("cover")) %>%
  pivot_wider(names_from =  c(family, inf_lab), values_from = value) %>%
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  ungroup() %>%
  colnames_to_row()


t2 =
  all_res_ng %>%
  filter(inf_level == 10, (snr_b == 0.5 | is.na(snr_b)), snr_eps == 1, In == 100, len == 50) %>%
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
  filter(inf_level == 10, type == "Strata scaling and noise", snr_eps == 1, In == 100, len == 50) %>%
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
  filter(inf_level == 10, type == "Strata scaling and noise", snr_eps == 1, snr_b == 0.5, len == 50) %>%
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
  filter(inf_level == 10, type == "Strata scaling and noise", snr_eps == 1, snr_b == 0.5, In == 100) %>%
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
           V1 == "unweighted" ~ "Unweighted",
           V1 == "weighted" ~ "Weighted",
           .default = V1
         )) %>%
  filter(V2 %in% c("Pointwise", "name")) %>%
  select(-V2) %>%
  kable(format = "latex", booktabs = TRUE)

tb %>%
  mutate(V2 = case_when(V2 == "cover_pw" ~ "Pointwise", V2 == "cover_joint" ~ "Joint", TRUE ~ V2),
         V1 = case_when(
           V1 == "unweighted" ~ "Unweighted",
           V1 == "weighted" ~ "Weighted",
           .default = V1
         )) %>%
  filter(V2 %in% c("Joint", "name")) %>%
  select(-V2) %>%
  kable(format = "latex", booktabs = TRUE)

### ----- finding large diffs -----

# getremote survey_FoSR/results/mfpca_sim2/*.rds results/mfpca_sim2/
source(here::here("R", "create_survey_settings.R"))
files = list.files(here::here("results", "mfpca_sim2"), full.names = TRUE)


res = map_dfr(.x = files, .f = \(x) read_rds(x) %>%
                mutate(fold = as.numeric(sub(".*fold\\_(.+).rds.*", "\\1", basename(x)))))


res = res %>%
  left_join(settings_new, by = "fold")

res = res %>%
  mutate(type = case_when(
    strata_scale == 0 & strata_sigma == 0 ~ "No random effects",
    strata_scale > 0 & strata_sigma > 0 ~ "Strata scaling and noise",
    strata_scale > 0 & strata_sigma == 0 ~ "Strata scaling only",
    strata_scale == 0 & strata_sigma > 0 ~ "Strata noise only"
  ))

p1 = res %>%
  mutate(type = factor(type),
         type = forcats::fct_relevel(type, "Strata scaling only", after = 1L)) %>%
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
p1 = res %>%
  mutate(type = factor(type),
         type = forcats::fct_relevel(type, "Strata scaling only", after = 1L),
         snr_b = if_else(type == "No random effects", "Not applicable", paste0(snr_b)),
         snr_b = factor(snr_b, levels = c("Not applicable", "0.5", "1", "5"),
                        labels = c("Not applicable", "High", "Medium", "Low"))) %>%
  ggplot(aes(x = pve_b * 100)) +
  geom_density(aes(fill = snr_b, color = snr_b), alpha = 0.5) +
  facet_grid(type ~ ., scales = "free_y") +
  scale_color_manual(values = c("#3FA0FFFF", "#D82632FF", "#F76D5EFF", "#FFAD72FF"), name = "Relative importance of random effects") +
  scale_fill_manual(values = c("#3FA0FFFF", "#D82632FF", "#F76D5EFF", "#FFAD72FF"), name = "Relative importance of random effects") +
  labs(x = "Percent variance explained by strata/PSU", y = "Density")  +
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
lo = res %>% arrange(desc(pve_w)) %>% slice(1:10)
hi = res %>% arrange(desc(pve_b)) %>% slice(1:10)
res_sort = res %>% arrange(desc(pve_w))


diffs_rwyb =
  all_res_g %>%
  filter(var == "x") %>%
  select(id, cover_pw, boot_type, fold) %>%
  pivot_wider(names_from = boot_type, values_from = cover_pw) %>%
  janitor::clean_names() %>%
  mutate(diff = rao_wu_yue_beaumont - weighted) %>%
  group_by(fold) %>%
  summarize(across(diff, mean)) %>%
  left_join(settings_new, by = join_by(fold)) %>%
  left_join(res %>% select(fold, pve_b, pve_w), by = join_by(fold))

p1 = diffs_rwyb %>%
  mutate(type = case_when(
    strata_scale == 0 & strata_sigma == 0 ~ "No random effects",
    strata_scale > 0 & strata_sigma > 0 ~ "Strata scaling and noise",
    strata_scale > 0 & strata_sigma == 0 ~ "Strata scaling only",
    strata_scale == 0 & strata_sigma > 0 ~ "Strata noise only"
  )) %>%
  mutate(type = factor(type),
         type = forcats::fct_relevel(type, "Strata scaling only", after = 1L),
         snr_b = if_else(type == "No random effects", "Not applicable", paste0(snr_b)),
         snr_b = factor(snr_b, levels = c("Not applicable", "0.5", "1", "5"),
                        labels = c("Not applicable", "High", "Medium", "Low"))) %>%
  ggplot(aes(x = pve_b * 100, y = diff * 100))  +
  geom_point(aes(color = snr_b, shape = type), size = 3) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_shape_discrete(name = "Random effects") +
  scale_color_manual(values = c("#3FA0FFFF", "#D82632FF", "#F76D5EFF", "#FFAD72FF"), name = "Importance of random effects") +
  geom_smooth(color = "darkgrey") +
  labs(x = "% Variation Explained by Strata and PSU", y = "Difference in Coverage: RWYB - Weighted") +
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

