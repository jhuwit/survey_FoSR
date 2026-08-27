library(dplyr)
library(readr)
colnames_to_row <- function(df) {
  r1 = colnames(df) %>% t() %>% as_tibble()
  colnames(df) = paste0("V", seq_len(ncol(df)))
  out = rbind(r1, df)
  # colnames(out) <- paste0("V", seq_len(ncol(df)))
  out
}


# paletteer_d("colorBlindness::paletteMartin")

### process results from "smooth first"

sfirst = read_rds(here("results", "simulations", "all_survey_res_smoothfirst.rds"))
regular = read_rds(here("results", "simulations", "all_survey_res.rds"))


folds_sfirst = sfirst$fold |> unique()

regular =
  regular |>
  filter(fold %in% folds_sfirst)


all =
  regular |>
  mutate(type = "smooth_second") |>
  bind_rows(sfirst |> mutate(type = "smooth_first")) |>
  mutate(lMISE = log10(MISE),
         label = glue::glue("Family: {family}, L: {len}"),
         type = factor(type, levels = c("smooth_second", "smooth_first"), labels = c("Smooth after estimation", "Smooth before estimation")))

# figure comparing MISE in smoothing before vs. after
p =
  all |>
  filter(boot_type == "weighted") |>
  mutate(len = factor(glue::glue("L = {len}"), levels = c("L = 50", "L = 100")),
         family = factor(family, levels = c("gaussian", "binomial", "poisson"), labels = c("Gaussian", "Binomial", "Poisson"))) |>
  ggplot(aes(x = var, y = lMISE, fill = type)) +
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
               # aes(shape = type),
               size = 2,
               position = position_dodge(width = 0.75)) +
  theme(strip.text = element_text(size = 14)) +
  # palette.color.discrete = cols[1:2],
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 14),
                   title = element_text(size = 15)) +
  theme_sub_axis(text= element_text(size = 12),
                 title = element_text(size = 15)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  facet_grid(family ~ len, scales = "free_y") +
  scale_x_discrete(labels = c("Intercept", "X"), name = "Variable") +
  labs(y = expression(log[10]~"Mean Integrated Squared Error"))  +
  scale_fill_manual(values = c("#B6DBFFFF", "#924900FF"), name = "Smooth order")
p

png(here::here("manuscript", "figures", "smooth_order.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()

## tables
all |>
  filter(var == "x", !grepl("weight", boot_type)) |>
  mutate(cover_pw = cover_pw / len) |>
  mutate(family = factor(family, levels = c("gaussian", "binomial", "poisson"), labels = c("Gaussian", "Binomial", "Poisson"))) |>
  group_by(len, family, boot_type, type) |>
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) |>
  pivot_longer(cols = starts_with("cover")) |>
  pivot_wider(names_from = c(type, family), values_from = value) |>
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row() |>
  filter(V3 %in% c("name", "cover_pw")) |>
  select(-V3) |>
  kable(format = "latex", booktabs = TRUE)

all |>
  filter(var != "x", !grepl("weight", boot_type)) |>
  mutate(cover_pw = cover_pw / len) |>
  mutate(family = factor(family, levels = c("gaussian", "binomial", "poisson"), labels = c("Gaussian", "Binomial", "Poisson"))) |>
  group_by(len, family, boot_type, type) |>
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) |>
  pivot_longer(cols = starts_with("cover")) |>
  pivot_wider(names_from = c(type, family), values_from = value) |>
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row() |>
  filter(V3 %in% c("name", "cover_pw")) |>
  select(-V3) |>
  kable(format = "latex", booktabs = TRUE)

all |>
  filter(var == "x", !grepl("weight", boot_type)) |>
  mutate(cover_pw = cover_pw / len) |>
  mutate(family = factor(family, levels = c("gaussian", "binomial", "poisson"), labels = c("Gaussian", "Binomial", "Poisson"))) |>
  group_by(len, family, boot_type, type) |>
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) |>
  pivot_longer(cols = starts_with("cover")) |>
  pivot_wider(names_from = c(type, family), values_from = value) |>
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row() |>
  filter(V3 %in% c("name", "cover_joint")) |>
  select(-V3) |>
  kable(format = "latex", booktabs = TRUE)

all |>
  filter(var != "x", !grepl("weight", boot_type)) |>
  mutate(cover_pw = cover_pw / len) |>
  mutate(family = factor(family, levels = c("gaussian", "binomial", "poisson"), labels = c("Gaussian", "Binomial", "Poisson"))) |>
  group_by(len, family, boot_type, type) |>
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) |>
  pivot_longer(cols = starts_with("cover")) |>
  pivot_wider(names_from = c(type, family), values_from = value) |>
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  colnames_to_row() |>
  filter(V3 %in% c("name", "cover_joint")) |>
  select(-V3) |>
  kable(format = "latex", booktabs = TRUE)

### process reults from "smooth type"

rm(list = ls())



ssens = read_rds(here("results", "simulations", "all_survey_res_smoothsens.rds")) |>
  select(MISE, cover_joint, cover_pw, var, family, len, nknots_spec, splines, method) |>
  mutate(cover_pw = cover_pw / len) |>
  mutate(len = factor(glue::glue("L = {len}"), levels = c("L = 50", "L = 100")),
         family = factor(family, levels = c("gaussian", "binomial", "poisson"), labels = c("Gaussian", "Binomial", "Poisson")),
         splines = factor(splines, levels = c("tp", "cr", "ps"), labels = c("tp (default)", "cr", "ps")))
# regular = read_rds(here("results", "simulations", "all_survey_res.rds")) |>
#   filter(boot_type == "BRR", family %in% c("binomial", "gaussian", "poisson"), len %in% c(50, 100),
#          snr_b == 0.5, strata_sigma == 0.05, strata_scale == 0.125, inf_level == 10, In == 100) |>
#   filter((family !="gaussian" & is.na(snr_eps)) | (family == "gaussian" & snr_eps == 1)) |>
#   select(MISE, cover_joint, cover_pw, var, family, len)

unique(ssens$nknots_spec)
unique(ssens$splines)
unique(ssens$method)

ssens |>
  select(nknots_spec, splines, method) |>
  distinct()
# len varied from 50,100
# family: binomial, gaussian, poisson
# default: current, "tp", GCV.Cp

p = ssens |>
  filter(nknots_spec == "current", method == "GCV.Cp") |>
  # filter(family == "gaussian", var == "x", len == 50) |>
  mutate(lMISE = log10(MISE)) |>
  ggplot(aes(x = var, y = lMISE, fill = splines)) +
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
               # aes(shape = type),
               size = 2,
               position = position_dodge(width = 0.75)) +
  theme(strip.text = element_text(size = 14)) +
  # palette.color.discrete = cols[1:2],
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 14),
                   title = element_text(size = 15)) +
  theme_sub_axis(text= element_text(size = 12),
                 title = element_text(size = 15)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  facet_grid(family ~ len, scales = "free_y") +
  scale_x_discrete(labels = c("Intercept", "X"), name = "Variable") +
  labs(y = expression(log[10]~"Mean Integrated Squared Error"))  +
  scale_fill_manual(values = c("#B6DBFFFF", "#920000FF"), name = "Type of smoother", labels = c("thin plate", "cubic"))

p
png(here::here("manuscript", "figures", "smooth_type.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()

p = ssens |>
  filter(splines == "tp (default)" & method == "GCV.Cp") |>
  mutate(lMISE = log10(MISE)) |>
  ggplot(aes(x = var, y = lMISE, fill = nknots_spec)) +
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
               # aes(shape = type),
               size = 2,
               position = position_dodge(width = 0.75)) +
  theme(strip.text = element_text(size = 14)) +
  # palette.color.discrete = cols[1:2],
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 14),
                   title = element_text(size = 15)) +
  theme_sub_axis(text= element_text(size = 12),
                 title = element_text(size = 15)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  facet_grid(family ~ len, scales = "free_y") +
  scale_x_discrete(labels = c("Intercept", "X"), name = "Variable") +
  labs(y = expression(log[10]~"Mean Integrated Squared Error"))  +
  scale_fill_manual(values = c("#006DDBFF", "#B6DBFFFF", "#004949FF"), name = "Number of knots")
p

png(here::here("manuscript", "figures", "nknots.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()




p = ssens |>
  filter(splines == "tp (default)" & nknots_spec == "current") |>
  mutate(lMISE = log10(MISE)) |>
  ggplot(aes(x = var, y = lMISE, fill = method)) +
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
               # aes(shape = type),
               size = 2,
               position = position_dodge(width = 0.75)) +
  theme(strip.text = element_text(size = 14)) +
  # palette.color.discrete = cols[1:2],
  theme_sub_legend(position = "bottom",
                   text = element_text(size = 14),
                   title = element_text(size = 15)) +
  theme_sub_axis(text= element_text(size = 12),
                 title = element_text(size = 15)) +
  theme_sub_panel(grid.major.x = element_blank(),
                  grid.minor.x = element_blank()) +
  facet_grid(family ~ len, scales = "free_y") +
  scale_x_discrete(labels = c("Intercept", "X"), name = "Variable") +
  labs(y = expression(log[10]~"Mean Integrated Squared Error"))  +
  scale_fill_manual(values = c("#B6DBFFFF", "#FF6DB6FF"), name = "Smoothing paramater estimation method")
p
png(here::here("manuscript", "figures", "est_type.png"),
    width = 10, height = 8, units = "in", res = 350)
p
dev.off()
## tables

ssens |>
  filter(var == "x", method == "GCV.Cp", nknots_spec == "current") |>
  group_by(len, family, splines) |>
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) |>
  pivot_longer(cols = starts_with("cover")) |>
  pivot_wider(names_from = c(splines, family), values_from = value) |>
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  arrange(name) |>
  kable(format = "latex", booktabs = TRUE)

ssens |>
  filter(var != "x", method == "GCV.Cp", nknots_spec == "current") |>
  group_by(len, family, splines) |>
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) |>
  pivot_longer(cols = starts_with("cover")) |>
  pivot_wider(names_from = c(splines, family), values_from = value) |>
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  arrange(name) |>
  kable(format = "latex", booktabs = TRUE)

ssens |>
  filter(var == "x", method == "GCV.Cp", splines == "tp (default)") |>
  group_by(len, family, nknots_spec) |>
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) |>
  pivot_longer(cols = starts_with("cover")) |>
  pivot_wider(names_from = c(nknots_spec, family), values_from = value) |>
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  arrange(name) |>
  kable(format = "latex", booktabs = TRUE)

ssens |>
  filter(var != "x", method == "GCV.Cp", splines == "tp (default)") |>
  group_by(len, family, nknots_spec) |>
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) |>
  pivot_longer(cols = starts_with("cover")) |>
  pivot_wider(names_from = c(nknots_spec, family), values_from = value) |>
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  arrange(name) |>
  kable(format = "latex", booktabs = TRUE)

ssens |>
  filter(var == "x", nknots_spec == "current", splines == "tp (default)") |>
  group_by(len, family, method) |>
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) |>
  pivot_longer(cols = starts_with("cover")) |>
  pivot_wider(names_from = c(method, family), values_from = value) |>
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  arrange(name) |>
  kable(format = "latex", booktabs = TRUE)

ssens |>
  filter(var != "x", nknots_spec == "current", splines == "tp (default)") |>
  group_by(len, family, method) |>
  summarize(cover_pw = mean(cover_pw),
            cover_joint = mean(cover_joint)) |>
  pivot_longer(cols = starts_with("cover")) |>
  pivot_wider(names_from = c(method, family), values_from = value) |>
  mutate(across(where(is.numeric), ~format(.x, digits = 2, nsmall = 2)),
         across(everything(), as.character)) %>%
  arrange(name) |>
  kable(format = "latex", booktabs = TRUE)
