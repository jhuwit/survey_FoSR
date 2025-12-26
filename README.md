# Function on Scalar Regression with Complex Survey Designs
Code to accompany manuscript

[Preprint](https://arxiv.org/abs/2511.05487)

[R package](https://jhuwit.github.io/svyfosr/)

Guide to repository: 
## R
+ `00_data_gen_function_ff.R`: generate functional survey data
+ `00_data_gen_function_fast.R`: generate functional survey data for 1440 simulation (faster for big data)
+ `01_sim_functions.R`: functions for running simulation
+ `02a_run_sim_nonsurvey.R`: run non-survey simulation
+ `02b_process_nonsurvey_sim_results.R`: process results from non-survey sim
+ `03a_get_population_mean.R`: get population means for results processing in simulation
+ `03b_run_sim_ff.R`: run simulation
+ `03c_process_survey_sim_results.R`: process results from survey sim
+ `04a_run_sim_empirical_1440.R`: run empirical simulation based on NHANES data
+ `04b_process_empirical_1440_sim_results.R`: process results from empirical sim
+ `05a_process_pa_data.R`: process NHANES PA data for application
+ `05b_application_1440.R`: run application to real data
+ `06a_run_mfpca_nhanes.R`: run MFPCA on NHANES data
+ `06b_run_mfpca_sim.R`: run MFPCA on simulated data
+ `06c_process_mfpca_sim_results.R`: process results from MFPCA sim
+ `07_tables_and_figures.R`: tables and figs for manuscript
+ `benchmark.R`: comparing sped-up GLMS vs original GLMS
+ `utils.R`: helpful functions for cluster 
+ `create_survey_settings.R`: create settings for survey sim
+ `create_nonsurvey_settings.R`: create settings for non-survey sim

## data 
+ `steps_covariates.rds`: data needed for empirical simulation and application 

## results/simulations
+ `all_empirical_sim_res.rds`: results from empirical simulation
+ `all_survey_res.rds`: results from survey simulation
+ `all_nonsurvey_sim_res.rds`: results from non-survey simulation

## vignettes
+ `Analytical_Sim.qmd`: quarto document describing analytical simulation
+ `Analytical_Simulation_Appendix.html`: HTML describing simulation, for manuscript
+ `Analytical_Simulation_Appendix.qmd`: quarto document to generate above html

`pipeline.sh`: bash code for running on cluster
