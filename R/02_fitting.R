library(tidyverse)
library(svylme)
sample_psu = read_rds(here::here("data", "sim", "sample_psu_sp.rds"))
sample = read_rds(here::here("data", "sim", "sample_sp.rds"))
sample_psu_noprob = read_rds(here::here("data", "sim", "sample_psu.rds"))
sample_noprob = read_rds(here::here("data", "sim", "sample.rds"))

sample %>%
  slice(1:100) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50)) %>%
  mutate(PSU = factor(PSU),
         y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(PSU) %>%
  mutate(mean_curve = mean(y_smooth)) %>%
  ggplot(aes(y = y_smooth, color = PSU, group = ID)) +
  geom_spaghetti(lwd = .7, alpha = .1) + # Individual curves
  geom_spaghetti(aes(y = mean_curve), lwd = 1.2)

sample_psu %>%
  slice(1:100) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50)) %>%
  mutate(PSU = factor(PSU),
         y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(PSU) %>%
  mutate(mean_curve = mean(y_smooth)) %>%
  ggplot(aes(y = y_smooth, color = PSU, group = ID)) +
  geom_spaghetti(lwd = .7, alpha = .1) + # Individual curves
  geom_spaghetti(aes(y = mean_curve), lwd = 1.2)

sample_psu_noprob %>%
  slice(1:100) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50)) %>%
  mutate(PSU = factor(PSU),
         y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(PSU) %>%
  mutate(mean_curve = mean(y_smooth)) %>%
  ggplot(aes(y = y_smooth, color = PSU, group = ID)) +
  geom_spaghetti(lwd = .7, alpha = .1) + # Individual curves
  geom_spaghetti(aes(y = mean_curve), lwd = 1.2)

sample_noprob %>%
  slice(1:100) %>%
  mutate(Y = matrix(Y, ncol = 50),
         Y_tf = tfd(Y, length = 50)) %>%
  mutate(PSU = factor(PSU),
         y_smooth = tf_smooth(Y_tf, method = "lowess")) %>%
  group_by(PSU) %>%
  mutate(mean_curve = mean(y_smooth)) %>%
  ggplot(aes(y = y_smooth, color = PSU, group = ID)) +
  geom_spaghetti(lwd = .7, alpha = .1) + # Individual curves
  geom_spaghetti(aes(y = mean_curve), lwd = 1.2)


## normalize weights??
fit_fui = function(df, weight = FALSE, svy = FALSE){
  df =
    df %>%
    mutate(wt_norm = weight / mean(weight))
  if(!weight & !svy){
    df %>%
      select(ID, X1, Y, wt_norm) %>%
      mutate(Y_tf = matrix(Y, ncol = 50),
              Y_tf = tfd(Y_tf, length = 50)) %>%
      tf_unnest(Y_tf) %>%
      select(-Y) %>%
      rename(epoch = Y_tf_arg, Y =Y_tf_value) %>%
      nest(data = -epoch) %>%
      mutate(
        model = map(.x = data,
                    .f = function(x){
                      fit = lm(Y ~ X1,
                               data = x)
                      fit}),
        result = map(model, broom::tidy)
      ) %>%
      select(epoch, result) %>%
      unnest(result)
  } else if(weight & !svy){
      df %>%
        select(ID, X1, Y, wt_norm) %>%
        mutate(Y_tf = matrix(Y, ncol = 50),
               Y_tf = tfd(Y_tf, length = 50)) %>%
        tf_unnest(Y_tf) %>%
        select(-Y) %>%
        rename(epoch = Y_tf_arg, Y =Y_tf_value) %>%
        nest(data = -epoch) %>%
        mutate(
          model = map(.x = data,
                      .f = function(x){
                        fit = lm(Y ~ X1,
                                 data = x,
                                 weights = wt_norm)
                        fit}),
          result = map(model, broom::tidy)
        ) %>%
        select(epoch, result) %>%
        unnest(result)
  } else if(weight & svy){
    df %>%
      select(ID, X1, Y, wt_norm, PSU, strata) %>%
      mutate(Y_tf = matrix(Y, ncol = 50),
             Y_tf = tfd(Y_tf, length = 50)) %>%
      tf_unnest(Y_tf) %>%
      select(-Y) %>%
      rename(epoch = Y_tf_arg, Y =Y_tf_value) %>%
      nest(data = -epoch) %>%
      mutate(
        model = map(.x = data,
                    .f = function(x){
                      des = survey::svydesign(
                        id = ~ PSU,
                        strata = ~ strata,
                        weights = ~ wt_norm,
                        data = x,
                        nest = TRUE
                      )
                      fit = svyglm(Y ~ X1,
                                   data = x,
                                   family = gaussian(),
                                   design = des)
                      fit}),
        result = map(model, broom::tidy)
      ) %>%
      select(epoch, result) %>%
      unnest(result)
  }
}


pt_est_raw = fit_fui(sample)
pt_est_wtd = fit_fui(sample, weight = TRUE)
pt_est_wtd_svy = fit_fui(sample, weight = TRUE, svy = TRUE)

fui_coef_df_raw  =
  pt_est_raw %>%
  select(epoch, term, coef = estimate) %>%
  tf_nest(coef, .id = term, .arg = epoch) %>%
  mutate(smooth_coef = tf_smooth(coef, method = "lowess")) %>%
  mutate(type = "raw")


fui_coef_df =
  pt_est_wtd %>%
  select(epoch, term, coef = estimate) %>%
  tf_nest(coef, .id = term, .arg = epoch) %>%
  mutate(smooth_coef = tf_smooth(coef, method = "lowess")) %>%
  mutate(type = "wtd")


fui_coef_df_svy  =
  pt_est_wtd_svy %>%
  select(epoch, term, coef = estimate) %>%
  tf_nest(coef, .id = term, .arg = epoch) %>%
  mutate(smooth_coef = tf_smooth(coef, method = "lowess")) %>%
  mutate(type = "svy")

# plot true functions
L = 50
grid <- seq(0, 1, length = L)
beta_true <- matrix(NA, 3, L)
beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
beta_true[2,] = dnorm(grid, .6, .15)/20

bt =
  beta_true %>%
  as.data.frame() %>%
  mutate(names = c("(Intercept)", "X1", "X2")) %>%
  pivot_longer(cols = starts_with("V")) %>%
  mutate(epoch = as.numeric(sub(".*V", "", name))) %>%
  rename(term = names) %>%
  filter(term != "X2") %>%
  select(-name)

true_beta =
  bt %>%
  select(coef = value, everything()) %>%
  tf_nest(coef, .id = term, .arg = epoch) %>%
  mutate(smooth_coef = tf_smooth(coef, method = "lowess"))


fui_coef_df %>%
  bind_rows(fui_coef_df_raw) %>%
  bind_rows(fui_coef_df_svy) %>%
  ggplot(aes(y = smooth_coef, color = type)) +
  geom_spaghetti(lwd = 1.2) +
  geom_spaghetti(aes(y = coef), alpha = .2)  +
  facet_wrap(~term, scales = "free") +
  scale_x_continuous(breaks=seq(0,50,10)) +
  geom_spaghetti(data = true_beta, aes(y = coef), lwd = 1.2, color = "black")


######## try with PSU
pt_est_raw = fit_fui(sample_psu)
pt_est_wtd = fit_fui(sample_psu, weight = TRUE)
pt_est_wtd_svy = fit_fui(sample_psu, weight = TRUE, svy = TRUE)

fui_coef_df_raw  =
  pt_est_raw %>%
  select(epoch, term, coef = estimate) %>%
  tf_nest(coef, .id = term, .arg = epoch) %>%
  mutate(smooth_coef = tf_smooth(coef, method = "lowess")) %>%
  mutate(type = "raw")


fui_coef_df =
  pt_est_wtd %>%
  select(epoch, term, coef = estimate) %>%
  tf_nest(coef, .id = term, .arg = epoch) %>%
  mutate(smooth_coef = tf_smooth(coef, method = "lowess")) %>%
  mutate(type = "wtd")


fui_coef_df_svy  =
  pt_est_wtd_svy %>%
  select(epoch, term, coef = estimate) %>%
  tf_nest(coef, .id = term, .arg = epoch) %>%
  mutate(smooth_coef = tf_smooth(coef, method = "lowess")) %>%
  mutate(type = "svy")

# plot true functions
L = 50
grid <- seq(0, 1, length = L)
beta_true <- matrix(NA, 3, L)
beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
beta_true[2,] = dnorm(grid, .6, .15)/20

bt =
  beta_true %>%
  as.data.frame() %>%
  mutate(names = c("(Intercept)", "X1", "X2")) %>%
  pivot_longer(cols = starts_with("V")) %>%
  mutate(epoch = as.numeric(sub(".*V", "", name))) %>%
  rename(term = names) %>%
  filter(term != "X2") %>%
  select(-name)

true_beta =
  bt %>%
  select(coef = value, everything()) %>%
  tf_nest(coef, .id = term, .arg = epoch) %>%
  mutate(smooth_coef = tf_smooth(coef, method = "lowess"))


fui_coef_df %>%
  bind_rows(fui_coef_df_raw) %>%
  bind_rows(fui_coef_df_svy) %>%
  ggplot(aes(y = smooth_coef, color = type)) +
  geom_spaghetti(lwd = 1.2) +
  geom_spaghetti(aes(y = coef), alpha = .2)  +
  facet_wrap(~term, scales = "free") +
  scale_x_continuous(breaks=seq(0,50,10)) +
  geom_spaghetti(data = true_beta, aes(y = coef), lwd = 1.2, color = "black")



