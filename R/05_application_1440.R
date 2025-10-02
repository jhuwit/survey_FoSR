## work on real world data result
library(tidyverse)
steps_df = read_rds(here::here("data", "steps_covariates.rds")) %>%
  mutate(SEQN = as.character(SEQN))
source(here::here("R", "01_sim_functions.R"))
n_boot = 500

L = 1440
steps_df =
  steps_df %>%
  mutate(race_cat =
           case_when(race_hispanic_origin == "Non-Hispanic White" ~ "NH White",
                     race_hispanic_origin == "Non-Hispanic Black" ~ "NH Black",
                     race_hispanic_origin %in% c("Other Hispanic", "Mexican American") ~ "Hispanic",
                     .default = "Other"),
          race_cat = factor(race_cat, levels = c("NH White", "NH Black", "Hispanic", "Other")))



fit_types = c('weighted', 'unweighted')
boot_types = c('BRR', 'weighted', 'unweighted')


data =
  steps_df %>%
  mutate(sex_bin = if_else(gender == "Female", 1, 0),
         weight = full_sample_2_year_mec_exam_weight / 2,
         age_cat = cut(age_in_years_at_screening, breaks=c(18, 30, 50, 65, 80), include.lowest = TRUE, right = FALSE)) %>%
  rename(psu = masked_variance_pseudo_psu,
         strata = masked_variance_pseudo_stratum)
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
betaHat = map(.x = betaTilde, .f = get_betahat, L = L)




data = data %>%
  mutate(ID = SEQN)

library(future)
library(furrr)
cores = parallelly::availableCores() - 1
plan(multisession, workers = cores)
tictoc::tic()
res = future_map(
  .x = boot_types,
  betaHat = betaHat[[1]],
  .f = run_boots_xfast,
  data = dat_full %>% mutate(ID = SEQN),
  family = "gaussian",
  num_boots = 500,
  seed = 2025,
  L = L,
  model_formula = as.formula(paste0('Y~', 'age_cat + sex_bin + race_cat')),
  .options = furrr_options(seed = TRUE)
)

fit_list = list(betaHat[[1]], betaHat[[1]], betaHat[[2]])

cis = future_map2(
  .x = res,
  .y = fit_list,
  .f = get_cis,
  L = L,
  smooth_for_ci = TRUE,
  .options = furrr_options(seed = TRUE)
)
elapsed_time = tictoc::toc(quiet = TRUE)
plan(sequential)
# 179.122 sec elapsed
# write_rds(cis, here::here("results", "application.rds"))
#
#
#
# # getremote survey_FoSR/results/application.rds results/
# cis = read_rds(here::here("results", "application.rds"))
# # res = map2(
# #   .x = boot_types,
# #   .y = betaHat,
# #   .f = run_boots,
# #   data = data,
# #   family = "poisson",
# #   num_boots = 10,
# #   set_seed = TRUE,
# #   seed = 2025,
# #   L = 49,
# #   model_formula = as.formula(paste0('Y~', 'age_cat + sex_bin + race_cat'))
# #   # .options = furrr_options(seed = TRUE)
# # )
# # cis = map2(
# #   .x = res,
# #   .y = betaHat,
# #   .f = get_cis,
# #   smooth_for_ci = TRUE,
# #   L = 49
# #   # .options = furrr_options(seed = TRUE)
# # )
# #
# # plan(sequential)
#
#
# plot = FALSE
# if(plot){
#   get_plot_df = function(sim_res) {
#     num_var = nrow(sim_res$betaHat)
#     map_dfr(
#       .x = 1:num_var,
#       .f = function(r) {
#         data.frame(
#           s = 1:length(sim_res$betaHat[r, ]),
#           beta = sim_res$betaHat[r, ],
#           lower = sim_res$betaHat[r, ] - 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
#           upper = sim_res$betaHat[r, ] + 2 * sqrt(diag(sim_res$betaHat.var[, , r])),
#           lower.joint = sim_res$betaHat[r, ] - sim_res$qn[r] *
#             sqrt(diag(sim_res$betaHat.var[, , r])),
#           upper.joint = sim_res$betaHat[r, ] + sim_res$qn[r] * sqrt(diag(sim_res$betaHat.var[, , r]))
#         ) %>%
#           mutate(name = sim_res$betaHat %>% rownames() %>% .[r])
#       }
#     )
#   }
#   plt_df_boot = map(cis, .f = get_plot_df)
#
#
#   plt_df_boot %>%
#     bind_rows(.id = "model") %>%
#     mutate(model = factor(model, labels = c("BRR", "Wtd", "Unwtd"))) %>%
#     ggplot(aes(x = s, y = beta, group = model)) +
#     geom_line(aes(color = model))  +
#     facet_wrap(.~name, scales = "free_y") +
#     geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = model), alpha = .3) +
#     geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = model), alpha = .1) +
#     labs(x = "L", y = "Estimate") +
#     theme_classic() +
#     geom_hline(aes(yintercept = 0))
#
#   plt_df_boot %>%
#     bind_rows(.id = "model") %>%
#     # mutate(across(c(beta, lower, upper, lower.joint, upper.joint), ~exp(.x))) %>%
#     mutate(model = factor(model, labels = c("BRR", "Wtd", "Unwtd"))) %>%
#     ggplot(aes(x = s, y = beta, group = model)) +
#     geom_line(aes(color = model))  +
#     facet_wrap(.~name, scales = "free_y") +
#     geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = model), alpha = .3) +
#     geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = model), alpha = .1) +
#     labs(x = "L", y = "Estimate") +
#     theme_classic() +
#     geom_hline(aes(yintercept = 1))
#
#   library(paletteer)
#
#   col1 = "#009E73FF"
#   col2 = "#E69F00FF"
#
#   p1 = plt_df_boot %>%
#     bind_rows(.id = "model") %>%
#     filter(model %in% c("1", "3")) %>%
#     filter(grepl("age", name)) %>%
#     mutate(model = factor(model, labels = c("BRR", "Unweighted")),
#            model = factor(model, levels = c("Unweighted", "BRR"))) %>%
#     mutate(name = factor(name, labels = paste0("Age ", c("30-49", "50-64", "65+")))) %>%
#     # mutate(across(c(beta, lower, upper, lower.joint, upper.joint), ~exp(.x))) %>%
#     mutate(model = factor(model, labels = c("BRR", "Unweighted")),
#            model = factor(model, levels = c("Unweighted", "BRR")))%>%
#     ggplot(aes(x = s, y = beta, group = model)) +
#     facet_grid(.~name, scales = "free_y") +
#     geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = model), alpha = .3) +
#     geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = model), alpha = .1) +
#     geom_line(aes(color = model), linewidth = 1.1) +
#     theme_classic() +
#     geom_hline(aes(yintercept = 0)) +
#     scale_y_continuous(breaks=seq(-100, 100, 15)) +
#     scale_x_continuous(breaks= seq(0, 96, 16), labels = c("", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
#     labs(x = "Hour of Day", y = "Steps", fill = "Model", color = "Model") +
#     theme(strip.text = element_text(size = 12),
#           palette.color.discrete = c(col2, col1),
#           palette.fill.discrete = c(col2, col1)) +
#     theme_sub_legend(position = "bottom",
#                      text = element_text(size = 12),
#                      title = element_text(size = 13)) +
#     theme_sub_axis(text= element_text(size = 12),
#                    title = element_text(size = 13))
#
#   png(here::here("manuscript", "figures", "steps_application3.png"),
#       width = 12, height = 8, res = 300, units = "in")
#   p1
#   dev.off()
#
#   p2 = plt_df_boot %>%
#     bind_rows(.id = "model") %>%
#     filter(model %in% c("1", "3")) %>%
#     filter(grepl("race", name)) %>%
#     mutate(name = factor(name, labels = c("Hispanic", "Non-Hispanic Black", "Other"))) %>%
#     # mutate(across(c(beta, lower, upper, lower.joint, upper.joint), ~exp(.x))) %>%
#     mutate(model = factor(model, labels = c("BRR", "Unweighted")),
#            model = factor(model, levels = c("Unweighted", "BRR"))) %>%
#     ggplot(aes(x = s, y = beta, group = model)) +
#     # geom_line(aes(color = model), linewidth = 1.1)  +
#     facet_grid(.~name, scales = "free_y") +
#     geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = model), alpha = .3) +
#     geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = model), alpha = .1) +
#     geom_line(aes(color = model), linewidth = 1.1)  +
#     theme_classic() +
#     geom_hline(aes(yintercept = 0)) +
#     scale_y_continuous(breaks=seq(-50, 50, 10)) +
#     scale_x_continuous(breaks= seq(0, 96, 16), labels = c("", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
#     labs(x = "Hour of Day", y = "Steps", fill = "Model", color = "Model") +
#     theme(strip.text = element_text(size = 12),
#           palette.color.discrete = c(col2, col1),
#           palette.fill.discrete = c(col2, col1)) +
#     theme_sub_legend(position = "bottom",
#                      text = element_text(size = 12),
#                      title = element_text(size = 13)) +
#     theme_sub_axis(text= element_text(size = 12),
#                    title = element_text(size = 13))
#   png(here::here("manuscript", "figures", "steps_application2.png"),
#       width = 12, height = 8, res = 300, units = "in")
#   p2
#   dev.off()
#   p3 = plt_df_boot %>%
#     bind_rows(.id = "model") %>%
#     filter(model %in% c("1", "3")) %>%
#     filter(grepl("sex", name)) %>%
#     mutate(name = factor(name, labels = c("Female sex"))) %>%
#     mutate(model = factor(model, labels = c("BRR", "Unweighted")),
#            model = factor(model, levels = c("Unweighted", "BRR"))) %>%
#     ggplot(aes(x = s, y = beta, group = model)) +
#     facet_wrap(.~name, scales = "free_y") +
#     geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = model), alpha = .3) +
#     geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = model), alpha = .1) +
#     geom_line(aes(color = model), linewidth = 1.1)  +
#     theme_classic() +
#     geom_hline(aes(yintercept = 0)) +
#     scale_y_continuous(breaks=seq(-50, 50, 10)) +
#     scale_x_continuous(breaks= seq(0, 96, 16), labels = c("", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
#     labs(x = "Hour of Day", y = "Steps", fill = "Model", color = "Model") +
#     theme(strip.text = element_text(size = 12),
#           palette.color.discrete = c(col2, col1),
#           palette.fill.discrete = c(col2, col1)) +
#     theme_sub_legend(position = "bottom",
#                      text = element_text(size = 12),
#                      title = element_text(size = 13)) +
#     theme_sub_axis(text= element_text(size = 12),
#                    title = element_text(size = 13))
#
#   p4 = plt_df_boot %>%
#     bind_rows(.id = "model") %>%
#     filter(model %in% c("1", "3")) %>%
#     # mutate(across(c(beta, lower, upper, lower.joint, upper.joint), ~exp(.x))) %>%
#     mutate(model = factor(model, labels = c("BRR", "Unweighted")),
#            model = factor(model, levels = c("Unweighted", "BRR"))) %>%
#     filter(grepl("cept", name)) %>%
#     mutate(name = factor(name, labels = c("Intercept"))) %>%
#     # mutate(across(c(beta, lower, upper, lower.joint, upper.joint), ~exp(.x))) %>%
#     ggplot(aes(x = s, y = beta, group = model)) +
#     facet_wrap(.~name, scales = "free_y") +
#     geom_ribbon(aes(x = s, ymin = lower, ymax = upper, fill = model), alpha = .3) +
#     geom_ribbon(aes(x = s, ymin = lower.joint, ymax = upper.joint, fill = model), alpha = .1) +
#     geom_line(aes(color = model), linewidth = 1.1)  +
#     theme_classic() +
#     scale_y_continuous(breaks=seq(-50, 200, 25)) +
#     scale_x_continuous(breaks= seq(0, 96, 16), labels = c("", "04:00", "08:00", "12:00", "16:00", "20:00", "24:00")) +
#     labs(x = "Hour of Day", y = "Steps", fill = "Model", color = "Model") +
#     theme(strip.text = element_text(size = 12),
#           palette.color.discrete = c(col2, col1),
#           palette.fill.discrete = c(col2, col1)) +
#     theme_sub_legend(position = "bottom",
#                      text = element_text(size = 12),
#                      title = element_text(size = 13)) +
#     theme_sub_axis(text= element_text(size = 12),
#                    title = element_text(size = 13))
#
#   library(patchwork)
#   png(here::here("manuscript", "figures", "steps_application1.png"),
#       width = 12, height = 8, res = 300, units = "in")
#   (p3 | p4) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
#   dev.off()
#   png(here::here("manuscript", "figures", "steps_application_all.png"),
#       width = 12, height = 10, res = 300, units = "in")
#   (p4 | p3) / p1 / p2  + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")
#   dev.off()
#
#
# }
