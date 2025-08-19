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
source("src/revised_func.R")
print('Functions have been loaded')
set.seed(1234)

results.lr = read_rds('data/nhanes_true_effect_lr.rds') # true coefficients


######## Estimate Confidence Interval ########
var_name = c('Intercept', 'Beta')
L = 1440
testL = 1300:1400
dat = read_rds('results/NI/NI_noOrigial_visual_only3_verify_raw.rds')
betaTilde_boot = dat$res[[1]]
betaHat = dat$betaHat[[1]]


testL = 1200:1400
stats = get_cis_revise0415(betaTilde_boot = betaTilde_boot[,testL,],
        betaHat = betaHat[,testL],
        L=length(testL),
        smooth_for_variance = FALSE,
        smooth_for_ci = FALSE,
        nknots_min = 80,
        nknots_min_cov = 80)

visualize_notCover(beta_true = t(results.lr |> select(-Variable))[, testL],
                   mod_output = stats,
                   L_domain = testL)


L = 1440
group_size = 180
values <- 1:L
groups <- rep(1:ceiling(L / group_size), each = group_size)[1:L]
L_list <- split(values, groups)

plots = list()
for (i in 1:length(L_list)){
  testL = L_list[[i]]
  stats = get_cis(betaTilde_boot = betaTilde_boot[,testL,],
                  betaHat = betaHat[,testL],
                  L=length(testL),
                  smooth_for_variance = TRUE,
                  smooth_for_ci = FALSE,
                  nknots_min = NULL,
                  nknots_min_cov = 35)
  
  plot = visualize_notCover(beta_true = t(results.lr |> select(-Variable))[, testL],
                     mod_output = stats,
                     L_domain = testL)
  
  plots[[i]] = plot
}

combined_plot <- wrap_plots(plots, nrow = 3, ncol = 3) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")







# get_coverage_stats(stats, 
#                    name = 'Preston',
#                    beta_true = t(results.lr |> select(-Variable))[, testL],
#                    L = length(testL))
