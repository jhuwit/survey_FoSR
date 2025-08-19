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

results.lr = read_rds('data/nhanes_true_effect_lr.rds') # true coefficients


######## Parameter Setting ########
num_strata = 30 # total strata
sampled_strata = 15 # number of strata
L = 150 # dimension of the functional domain
B = 100 # bootstrap number
family = "gaussian" # distribution of Y
boot_types = c('Preston', 'weighted', 'none')
fit_types = c('design', 'weighted', 'none')

# boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'Preston', 'weighted', 'none')
# fit_types = c('design', 'design', 'design', 'weighted', 'none')

X_des = cbind(intercept = 0, nhanes_df |> select(Age))
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

saveRDS(stats, file = paste0("results/NI/sim_noOrigial_survey_L150_", sim_id, ".rds"))


# 
# sim_res <- list() ## store simulation results
# for (iter in 1:nsim) {
#   set.seed(iter)
#   x = try({
#     data <- simulate_nhanes_survey_fosr(X_des, Y_obs, L = L, sampled_strata = sampled_strata)
#     
#     betaTilde = map(.x = fit_types, .f = get_betatilde, data = data, family = family)
#     betaHat = map(.x = betaTilde, get_betahat, L = L)
#     if (parallel) {
#       plan(multisession, workers = ncores)
#       res = future_map2(
#         .x = boot_types,
#         .y = betaHat,
#         .f = run_boots,
#         data = data,
#         # family = family,
#         num_boots = B,
#         set_seed = TRUE,
#         seed = 2025,
#         L = L,
#         .options = furrr_options(seed = TRUE)
#       )
#       cis = future_map2(
#         .x = res,
#         .y = betaHat,
#         .f = get_cis,
#         L = L,
#         smooth_for_ci = FALSE,
#         .options = furrr_options(seed = TRUE)
#       )
#       plan(sequential)
#     } else{
#       res = map2(
#         .x = boot_types,
#         .y = betaHat,
#         .f = run_boots,
#         data = data,
#         # family = family,
#         num_boots = B,
#         set_seed = TRUE,
#         seed = 2025,
#         L = L)
#       cis = map2(
#         .x = res,
#         .y = betaHat,
#         .f = get_cis,
#         smooth_for_ci = TRUE,
#         L = L)
#     }
#     stats = map2(
#       .x = cis,
#       .y = boot_types,
#       .f = get_coverage_stats,
#       L = L,
#       beta_true = t(results.lr |> select(-Variable))[,1:L]
#     ) %>%
#       list_rbind() %>%
#       mutate(n = nrow(data), n_boot = B)
#     sim_res[[iter]] <- stats
#   })
#   rm(x)
#   print(iter)
# }
# 
# result =
#   sim_res %>%
#   keep(., is.data.frame) %>%
#   list_rbind(names_to = "id")
# 
# write_rds(result, 'results/sim_NI_noOrigial_survey_results.rds')



# # Load functions for sourced functions
# if (!requireNamespace("fastFMM", quietly = TRUE)) {
#   if (!requireNamespace("devtools", quietly = TRUE)) {
#     install.packages("devtools")
#   }
#   devtools::install_github("gloewing/fastFMM")
# }
# library(fastFMM)