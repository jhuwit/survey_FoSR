######## Load Packages ########
packages <- c(
  "dplyr", "knitr", "readxl", "readr", "haven", "visdat", "sjstats",
  "openxlsx", "stargazer", "broom", "ggplot2", "sjPlot", "stringr",
  "caret", "future", "purrr", "furrr", "here", "svrep", "patchwork",
  'fastFMM'
)

# Check which packages are not installed
packages_to_install <- packages[!(packages %in% installed.packages()[,"Package"])]

# Install and load missing packages
if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}
lapply(packages, library, character.only = TRUE)
print('Pakckages have been loaded')

source("src/00_sim_data_function.R")
source("src/01_sim_functions.R")
print('Functions have been loaded')


######## Load Parameters ########
args <- commandArgs(trailingOnly = TRUE)[1]
sim_id <- as.numeric(args[1])

cat("Running simulation", sim_id, "\n")


######## Load Data ########
nhanes_df = readRDS('data/pa_df_persub.rds') |> 
  dplyr::rename(Age = age_in_years_at_screening) |>
  dplyr::rename(weight = full_sample_2_year_mec_exam_weight) # raw data

results.lr = read_rds('data/nhanes_true_effect_lr_smoothed_bspline.rds') # true coefficients

######## Parameter Setting ########
num_strata = 30 # total strata
sampled_strata = 15 # number of strata
L = 150 # dimension of the functional domain
B = 100 # bootstrap number
family = "gaussian" # distribution of Y
boot_types = c('Rao-Wu-Yue-Beaumont', 'Preston')
fit_types = c('design', 'design')

# boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'Preston', 'weighted', 'none')
# fit_types = c('design', 'design', 'design', 'weighted', 'none')

X_des = cbind(intercept = 0, nhanes_df |> dplyr::select(Age))
Y_obs = nhanes_df[,which(grepl('min_', colnames(nhanes_df)))[1:L] ]


######## Run Simulations ########
set.seed(sim_id)
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
cis = map2(
  .x = res,
  .y = betaHat,
  .f = get_cis,
  smooth_for_ci = TRUE,
  L = L
)

stats = map2(
  .x = cis,
  .y = boot_types,
  .f = get_coverage_stats,
  L = L,
  beta_true = t(results.lr |> select(-Variable))[, 1:L]
) %>%
  list_rbind() %>%
  mutate(n = nrow(data), n_boot = B, sim = sim_id)

saveRDS(stats, file = paste0("results/NI/case01_L150", sim_id, ".rds"))