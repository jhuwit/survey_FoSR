# survey_FoSR
## R
+ `00_data_gen_function_ff.R`: generate functional survey data without informative sampling
+ `01_sim_functions.R`: functions for running simulation 
+ `03_run_sim_survey.R`: run simulation
+ `03_run_sim_nonsurvey.R`: run non-survey simulation
+ `03b_proc_res_nonsurvey.R`: process results from non-survey sim; make plots and tables 
+ `03b_proc_res_survey.R`: process results from survey sim; make plots and tables
+ `04_run_sim_empirical.R`: run empirical simulation based on NHANES data
+ `04b_proc_res_empirical.R`: process results from empirical sim; make plots and tables
+ `05_run_application.R`: run application to real data
+ `06_run_mfpca_nhanes.R`: run MFPCA on NHANES data
+ `06_run_mfpca_sim.R`: run MFPCA on simulated data
+ `06b_proc_res_mfpca_sim.R`: process results from MFPCA sim; make plots and tables

+ `benchmark.R`: comparing sped-up GLMS vs original GLMS
+ `utils.R`: helpful functions for cluster 
+ `create_survey_settings.R`: create settings for survey sim
+ `create_nonsurvey_settings.R`: create settings for non-survey sim
+ `process_steps_data.R`: process NHANES data to get steps data at participant level 
+ `create_nhanes_df.R`: process NHANES data to get MIMS at participant level for empirical simulation 

## data 
+ `steps_covariates.rds`: data needed for empirical simulation and application 
+ `covariates_accel_mortality_df.rds`
+ `pa_df_persub.rds`

## results
+ `application.rds`

### simulations 
+ `all_survey_sim_res.rds`
## vignettes

