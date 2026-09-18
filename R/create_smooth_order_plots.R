library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)


######
# make plots
######
ff = here::here("results", "simulations", "smooth_order_consistency.rds")
if (!file.exists(ff)) source(here::here("R", "smooth_order_show_inconsistent.R")) else n_result = read_rds(ff)

L = 50
beta_fixed  = matrix(NA, 2, L)
beta_fixed[1, ]  = -0.15 - 0.1 * sin(2 * pi * grid) - 0.1 * cos(2 * pi * grid)
beta_fixed[2, ]  = dnorm(grid, 0.6, 0.15) / 20

bt_df =
  beta_fixed |>
  as_tibble() |>
  mutate(var = c("b0", "b1")) |>
  pivot_longer(cols = -var, names_to = "s", names_transform = ~as.numeric(sub(".*V", "", .x)))

n_result |>
  pivot_longer(cols = contains("b")) |>
  filter(grepl("bias", name)) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1") |>
  mutate(smooth_order = factor(smooth_order, labels = c("Smooth first", "Smooth second"))) |>
  ggplot(aes(x = s, y = value, color = factor(n))) +
  geom_line() +
  facet_wrap(.~smooth_order) +
  labs(x = "Functional Domain", y = "Bias", title = bquote("Bias for"~beta[1])) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "darkgrey") +
  # scale_color_viridis_d(name = "Sample size", option = "B") +
  scale_color_manual(values = c("#FF7F00FF", "#19B2FFFF", "#654CFFFF", "#E51932FF"), name = "Sample Size") +
  theme(legend.position = c(0.9, 0.2))

p1 = n_result |>
  pivot_longer(cols = contains("b")) |>
  filter(grepl("bias", name)) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1") |>
  mutate(smooth_order = factor(smooth_order, labels = c("Smooth first", "Smooth second"))) |>
  ggplot(aes(x = s, y = value, color = factor(n))) +
  geom_line() +
  facet_wrap(.~smooth_order, scales = "free_y") +
  labs(x = "Functional Domain", y = "Bias", title = bquote("Bias for"~beta[1])) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "darkgrey") +
  # scale_color_viridis_d(name = "Sample size", option = "B") +
  scale_color_manual(values = c("#FF7F00FF", "#19B2FFFF", "#654CFFFF", "#E51932FF"), name = "Sample Size") +
  theme(legend.position = c(0.9, 0.2))


n_result |>
  pivot_longer(cols = contains("b")) |>
  filter(grepl("bias", name)) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0") |>
  mutate(smooth_order = factor(smooth_order, labels = c("Smooth first", "Smooth second"))) |>
  ggplot(aes(x = s, y = value, color = factor(n))) +
  geom_line() +
  facet_wrap(.~smooth_order) +
  labs(x = "Functional Domain", y = "Bias", title = bquote("Bias for"~beta[0])) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "darkgrey") +
  # scale_color_viridis_d(name = "Sample size", option = "B") +
  scale_color_manual(values = c("#FF7F00FF", "#19B2FFFF", "#654CFFFF", "#E51932FF"), name = "Sample Size") +
  theme(legend.position = c(0.9, 0.2))

p2 = n_result |>
  pivot_longer(cols = contains("b")) |>
  filter(grepl("bias", name)) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0") |>
  mutate(smooth_order = factor(smooth_order, labels = c("Smooth first", "Smooth second"))) |>
  ggplot(aes(x = s, y = value, color = factor(n))) +
  geom_line() +
  facet_wrap(.~smooth_order, scales = "free_y") +
  labs(x = "Functional Domain", y = "Bias", title = bquote("Bias for"~beta[0])) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "darkgrey") +
  # scale_color_viridis_d(name = "Sample size", option = "B") +
  scale_color_manual(values = c("#FF7F00FF", "#19B2FFFF", "#654CFFFF", "#E51932FF"), name = "Sample Size") +
  theme(legend.position = c(0.9, 0.2))


png(here::here("manuscript", "figures", "smooth_order_consistency.png"), width = 10, height = 6, units = "in", res = 350)
(p2 / p1) + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
dev.off()


p1 = n_result |>
  filter(n == 1000) |>
  select(s, contains("mean")) |>
  pivot_longer(cols = -s)  |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0") |>
  ggplot(aes(x = s, y = value, color = smooth_order)) +
  geom_line(linewidth = .8) +
  geom_line(data = bt_df |> filter(var == "b0") |> mutate(smooth_order = "truth"), linewidth = .8, linetype = "dashed") +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF", "grey50"), name = "",
                     labels = c("Smooth first", "Smooth second", "Truth")) +
  labs(x = "Functional Domain", y = "Value", title = bquote("True and estimated"~beta[0]~"at n = 1000")) +
  theme(legend.position = c(0.8, 0.2),
        legend.title = element_blank())

p2 = n_result |>
  filter(n == 1000) |>
  select(s, contains("mean")) |>
  pivot_longer(cols = -s)  |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1") |>
  ggplot(aes(x = s, y = value, color = smooth_order)) +
  geom_line(linewidth = .8) +
  geom_line(data = bt_df |> filter(var == "b1") |> mutate(smooth_order = "truth"), linewidth = .8, linetype = "dashed") +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF", "grey50"), name = "",
                     labels = c("Smooth first", "Smooth second", "Truth")) +
  labs(x = "Functional Domain", y = "Value", title = bquote("True and estimated"~beta[1]~ "at n = 1000")) +
  theme(legend.position = c(0.6, 0.2),
        legend.title = element_blank())



png(here::here("manuscript", "figures", "smooth_order_bias.png"), width = 10, height = 6, units = "in", res = 350)
(p1 / p2) + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
dev.off()

# First, get the bias data and determine scaling
bias_data = n_result |>
  filter(n == 1000) |>
  select(s, contains("bias")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1")

# Scaling parameters: map bias range to ~[0, 0.5] so it doesn't overlap too much with coverage
bias_range = range(bias_data$value)
scale_factor = 0.4 / max(abs(bias_range))  # scale to fit in lower part of plot
offset = 0.5  # center the bias around 0.5 on the primary axis

# Coverage data
cover_data = n_result |>
  filter(n == 1000) |>
  select(s, contains("cover")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1")

# Plot
p1 = ggplot() +
  # Coverage lines (primary axis)
  geom_line(data = cover_data,
            aes(x = s, y = value, color = smooth_order),
            linewidth = 1) +
  # Bias lines (scaled to primary axis)
  geom_line(data = bias_data,
            aes(x = s, y = value * scale_factor + offset, linetype = smooth_order),
            color = "grey40", linewidth = 0.8) +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF"), name = "Coverage",
                     labels = c("Smooth sirst", "Smooth second")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Bias",
                        labels = c("Smooth first", "Smooth second")) +
  scale_y_continuous(
    name = "Pointwise coverage",
    limits = c(0, 1),
    sec.axis = sec_axis(~ (. - offset) / scale_factor, name = "Bias")
  ) +
  labs(x = "Functional Domain",
       title = bquote("Empirical coverage for" ~ beta[1] ~ "at n = 1000")) +
  geom_hline(yintercept = 0.95, color = "darkgrey", linetype = "dashed") +
  geom_hline(yintercept = offset, color = "grey40", linetype = "dotted", alpha = 0.5) +  # zero line for bias
  annotate(geom = "text", label = "Nominal = 0.95", x = 5, y = 0.95,
           color = "darkgrey", size = 3, vjust = -1) +
  theme(legend.position = "bottom")

bias_data = n_result |>
  filter(n == 1000) |>
  select(s, contains("bias")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0")

# Scaling parameters: map bias range to ~[0, 0.5] so it doesn't overlap too much with coverage
bias_range = range(bias_data$value)
scale_factor = 0.4 / max(abs(bias_range))  # scale to fit in lower part of plot
offset = 0.5  # center the bias around 0.5 on the primary axis

# Coverage data
cover_data = n_result |>
  filter(n == 1000) |>
  select(s, contains("cover")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0")

# Plot
p2 = ggplot() +
  # Coverage lines (primary axis)
  geom_line(data = cover_data,
            aes(x = s, y = value, color = smooth_order),
            linewidth = 1) +
  # Bias lines (scaled to primary axis)
  geom_line(data = bias_data,
            aes(x = s, y = value * scale_factor + offset, linetype = smooth_order),
            color = "grey40", linewidth = 0.8) +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF"), name = "Coverage",
                     labels = c("Smooth sirst", "Smooth second")) +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Bias",
                        labels = c("Smooth first", "Smooth second")) +
  scale_y_continuous(
    name = "Pointwise coverage",
    limits = c(0, 1),
    sec.axis = sec_axis(~ (. - offset) / scale_factor, name = "Bias")
  ) +
  labs(x = "Functional Domain",
       title = bquote("Empirical coverage for" ~ beta[0] ~ "at n = 1000")) +
  geom_hline(yintercept = 0.95, color = "darkgrey", linetype = "dashed") +
  geom_hline(yintercept = offset, color = "grey40", linetype = "dotted", alpha = 0.5) +  # zero line for bias
  annotate(geom = "text", label = "Nominal = 0.95", x = 5, y = 0.95,
           color = "darkgrey", size = 3, vjust = -1) +
  theme(legend.position = "bottom")

png(here::here("manuscript", "figures", "smooth_order_cover.png"), width = 10, height = 6, units = "in", res = 350)
(p2 / p1) + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
dev.off()

### next do the fix EDF part

edf_fix = read_rds(here::here("results", "simulations", "smooth_order_fix_edf.rds"))
p1 =
  edf_fix |>
  pivot_longer(cols = contains("b")) |>
  filter(grepl("bias", name)) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1") |>
  mutate(smooth_order = factor(smooth_order, levels = c("sf", "ss", "sff"), labels = c("Smooth first", "Smooth second", "Smooth first calibrated EDF"))) |>
  ggplot(aes(x = s, y = value, color = factor(n))) +
  geom_line() +
  facet_wrap(.~smooth_order, scales = "free_y") +
  labs(x = "Functional Domain", y = "Bias", title = bquote("Bias for"~beta[1])) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "darkgrey") +
  # scale_color_viridis_d(name = "Sample size", option = "B") +
  scale_color_manual(values = c("#FF7F00FF", "#19B2FFFF", "#654CFFFF"), name = "Sample Size") +
  theme(legend.position = c(0.9, 0.2))

p2 =
  edf_fix |>
  pivot_longer(cols = contains("b")) |>
  filter(grepl("bias", name)) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0") |>
  mutate(smooth_order = factor(smooth_order, levels = c("sf", "ss", "sff"), labels = c("Smooth first", "Smooth second", "Smooth first calibrated EDF"))) |>
  ggplot(aes(x = s, y = value, color = factor(n))) +
  geom_line() +
  facet_wrap(.~smooth_order, scales = "free_y") +
  labs(x = "Functional Domain", y = "Bias", title = bquote("Bias for"~beta[0])) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "darkgrey") +
  # scale_color_viridis_d(name = "Sample size", option = "B") +
  scale_color_manual(values = c("#FF7F00FF", "#19B2FFFF", "#654CFFFF"), name = "Sample Size") +
  theme(legend.position = c(0.9, 0.2))



png(here::here("manuscript", "figures", "smooth_order_consistency_calibrated.png"), width = 10, height = 6, units = "in", res = 350)
(p2 / p1) + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
dev.off()

# First, get the bias data and determine scaling
bias_data = edf_fix |>
  filter(n == 1000) |>
  select(s, contains("bias")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1")

# Scaling parameters: map bias range to ~[0, 0.5] so it doesn't overlap too much with coverage
bias_range = range(bias_data$value)
scale_factor = 0.4 / max(abs(bias_range))  # scale to fit in lower part of plot
offset = 0.5  # center the bias around 0.5 on the primary axis

# Coverage data
cover_data = edf_fix |>
  filter(n == 1000) |>
  select(s, contains("cover")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b1")

# Plot
p1 = ggplot() +
  # Coverage lines (primary axis)
  geom_line(data = cover_data,
            aes(x = s, y = value, color = smooth_order),
            linewidth = 1) +
  # Bias lines (scaled to primary axis)
  geom_line(data = bias_data,
            aes(x = s, y = value * scale_factor + offset, linetype = smooth_order),
            color = "grey40", linewidth = 0.8) +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF", "#009E73FF"), name = "Coverage",
                     labels = c("Smooth first", "Smooth second", "Smooth first calibrated EDF")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Bias",
                        labels = c("Smooth first", "Smooth second", "Smooth first calibrated EDF")) +
  scale_y_continuous(
    name = "Pointwise coverage",
    limits = c(0, 1),
    sec.axis = sec_axis(~ (. - offset) / scale_factor, name = "Bias")
  ) +
  labs(x = "Functional Domain",
       title = bquote("Empirical coverage for" ~ beta[1] ~ "at n = 1000")) +
  geom_hline(yintercept = 0.95, color = "darkgrey", linetype = "dashed") +
  geom_hline(yintercept = offset, color = "grey40", linetype = "dotted", alpha = 0.5) +  # zero line for bias
  annotate(geom = "text", label = "Nominal = 0.95", x = 5, y = 0.95,
           color = "darkgrey", size = 3, vjust = -1) +
  theme(legend.position = "bottom")

bias_data = edf_fix |>
  filter(n == 1000) |>
  select(s, contains("bias")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0")

# Scaling parameters: map bias range to ~[0, 0.5] so it doesn't overlap too much with coverage
bias_range = range(bias_data$value)
scale_factor = 0.4 / max(abs(bias_range))  # scale to fit in lower part of plot
offset = 0.5  # center the bias around 0.5 on the primary axis

# Coverage data
cover_data = edf_fix |>
  filter(n == 1000) |>
  select(s, contains("cover")) |>
  pivot_longer(cols = -s) |>
  separate_wider_delim(cols = name, delim = "_", names = c("smooth_order", "xx", "var")) |>
  filter(var == "b0")

# Plot
p2 = ggplot() +
  # Coverage lines (primary axis)
  geom_line(data = cover_data,
            aes(x = s, y = value, color = smooth_order),
            linewidth = 1) +
  # Bias lines (scaled to primary axis)
  geom_line(data = bias_data,
            aes(x = s, y = value * scale_factor + offset, linetype = smooth_order),
            color = "grey40", linewidth = 0.8) +
  scale_color_manual(values = c("#E69F00FF", "#0072B2FF", "#009E73FF"), name = "Coverage",
                     labels = c("Smooth first", "Smooth second", "Smooth first calibrated EDF")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Bias",
                        labels = c("Smooth first", "Smooth second", "Smooth first calibrated EDF")) +
  scale_y_continuous(
    name = "Pointwise coverage",
    limits = c(0, 1),
    sec.axis = sec_axis(~ (. - offset) / scale_factor, name = "Bias")
  ) +
  labs(x = "Functional Domain",
       title = bquote("Empirical coverage for" ~ beta[0] ~ "at n = 1000")) +
  geom_hline(yintercept = 0.95, color = "darkgrey", linetype = "dashed") +
  geom_hline(yintercept = offset, color = "grey40", linetype = "dotted", alpha = 0.5) +  # zero line for bias
  annotate(geom = "text", label = "Nominal = 0.95", x = 5, y = 0.95,
           color = "darkgrey", size = 3, vjust = -1) +
  theme(legend.position = "bottom")

png(here::here("manuscript", "figures", "smooth_order_cover_calibrated.png"), width = 10, height = 6, units = "in", res = 350)
(p2 / p1) + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
dev.off()

# mean EDFs

edf_fix |>
  group_by(n) |>
  summarize(across(contains("edf"), first))
