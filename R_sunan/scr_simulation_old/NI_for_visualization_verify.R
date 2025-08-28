######## Load Packages ########
packages <- c(
  "dplyr", "knitr", "readxl", "readr", "haven", "visdat", "sjstats",
  "openxlsx", "stargazer", "broom", "ggplot2", "sjPlot", "stringr",
  "caret", "future", "purrr", "furrr", "here", "svrep", "patchwork",
  'fastFMM'
)

# Check packages
packages_to_install <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(packages_to_install) > 0) {
  install.packages(packages_to_install) # Install and load missing packages
}
lapply(packages, library, character.only = TRUE)
print('Pakckages have been loaded')


######## Load Functions ########
source("src/00_sim_data_function.R")
source("src/01_sim_functions.R")
print('Functions have been loaded')


######## Load Data ########
nhanes_df = readRDS('data/pa_df_persub.rds') |> 
  dplyr::rename(Age = age_in_years_at_screening) |>
  dplyr::rename(weight = full_sample_2_year_mec_exam_weight) # raw data

results.lr = read_rds('data/nhanes_true_effect_lr.rds') # true coefficients
results.lr.tmp = read_rds('data/nhanes_true_effect_lr.rds') |>
  mutate(Minutes = as.numeric(sub('min_', '', Variable)) ) |>
  select(-Variable)


######## Global Parameter ########
num_strata = 30 # total strata
sampled_strata = 15 # number of strata
L = 1440 # dimension of the functional domain
B = 100 # bootstrap number
family = "gaussian" # distribution of Y
set.seed(1234)

# boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'Preston', 'weighted', 'none')
# fit_types = c('design', 'design', 'design', 'weighted', 'none')

boot_types = c('Preston')
fit_types = c('design')


######## Load Parameter ########
args <- commandArgs(trailingOnly = TRUE)[1]
sim_id <- as.numeric(args[1])
cat("Running Patterns", sim_id, "\n")

######## Run Simulations ########

######## NI_Origial_survey ########
######## Special Data Prepare #######
if (sim_id == 1) {
  X_des = cbind(intercept = 0, nhanes_df |> select(Age, weight))
  X_des$weight = scale(X_des$weight)
  Y_obs = nhanes_df[,which(grepl('min_', colnames(nhanes_df)))[1:L] ]
  
  data <- simulate_nhanes_survey_fosr_w(X_des, Y_obs, L = L, sampled_strata = sampled_strata)
  betaTilde = map(.x = fit_types, .f = get_betatilde, data = data, family = family)
  betaHat = map(.x = betaTilde, get_betahat, L = L)
  
  res = map2(
    .x = boot_types,
    .y = betaHat,
    .f = run_boots,
    data = data,
    # family = family,
    num_boots = B,
    set_seed = TRUE,
    seed = 2025,
    L = L
  )
  # cis = map2(
  #   .x = res,
  #   .y = betaHat,
  #   .f = get_cis,
  #   smooth_for_ci = TRUE,
  #   L = L
  # )
  
  tmp_dat = list()
  tmp_dat[['res']] = res
  tmp_dat[['betaHat']] = betaHat
  
  saveRDS(tmp_dat, file = "results/NI/NI_Origial_visual_only3_verify.rds")
}



######## NI_noOrigial_survey ########
######## Special Data Prepare #######
if (sim_id == 2){
  X_des = cbind(intercept = 0, nhanes_df |> select(Age))
  Y_obs = nhanes_df[,which(grepl('min_', colnames(nhanes_df)))[1:L] ]
  
  data <- simulate_nhanes_survey_fosr(X_des, Y_obs, L = L, sampled_strata = sampled_strata)
  betaTilde = map(.x = fit_types, .f = get_betatilde, data = data, family = family)
  betaHat = map(.x = betaTilde, get_betahat, L = L)
  
  res = map2(
    .x = boot_types,
    .y = betaHat,
    .f = run_boots,
    data = data,
    # family = family,
    num_boots = B,
    set_seed = TRUE,
    seed = 2025,
    L = L
  )
  # cis = map2(
  #   .x = res,
  #   .y = betaHat,
  #   .f = get_cis,
  #   smooth_for_ci = TRUE,
  #   L = L
  # )
  
  tmp_dat = list()
  tmp_dat[['res']] = res
  tmp_dat[['betaHat']] = betaHat
  
  saveRDS(tmp_dat, file = "results/NI/NI_noOrigial_visual_only3_verify.rds")
}



# 
# betaHat = read_rds('/Users/gsn/NO_noOrigial_visual.rds')
# df = as.data.frame(t(betaHat[[4]])) |>
#   rename(Beta = X, Intercept = `(Intercept)`) |>
#   mutate(Minutes = row_number() ) |>
#   merge(results.lr.tmp,
#         by = 'Minutes',
#         suffixes = c("_estimate","_lr"))
# 
# df |>
#   select(Intercept_estimate, Beta_estimate, Minutes) |>
#   mutate(type = "estimate") |>
#   rename(Intercept = Intercept_estimate, Beta = Beta_estimate) |>
#   bind_rows(
#     df |>
#       select(Intercept_lr, Beta_lr, Minutes) |>
#       mutate(type = "lr") |>
#       rename(Intercept = Intercept_lr, Beta = Beta_lr)
#   ) |>
#   tidyr::pivot_longer(cols = c(Intercept, Beta),
#                       names_to = "Metric",
#                       values_to = "Value") |>
#   ggplot(aes(x = Minutes, y = Value, color = type)) +
#   geom_point(alpha = 0.6, size = 0.5) +       # 绘制数据点
#   facet_wrap(~ Metric, scales = "free_y") +
#   labs(title = "True effects of beta and intercept throught the day",
#        x = "Minutes",
#        y = "Value") +
#   theme(text = element_text(size = 14))
