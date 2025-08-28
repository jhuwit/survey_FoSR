# survey_FoSR
## R
+ `00_data_gen_function.R`: generate functional survey data without informative sampling
+ `00_data_gen_function_inf.R`: generate  functional survey data with informative sampling
+ `00_data_gen_function_ml.R`: generate multilevel functional survey data without informative sampling
+ `01_sim_functions.R`: functions for running simulation 
+ `01_sim_functions.R`: functions for running multilevel simulation 
+ `02_make_sim_grid.R`: make a dataframe of different simulation combos for running easily via folds 
+ `03_run_simulation.R`: run simulation (on cluster) for different settings
+ `03_run_simulation_inf.R`: run simulation under informative sampling (on cluster) for different settings 
+ `03_run_simulation_ml.R`: run simulation with multilevel data 

+ `utils.R`: helpful functions for cluster 


## R_sunan

This repository contains functions and results for empirical simulations derived from NHANES data.

### Folder Structure

- **`/data/`**  
  Raw NHANES data and the estimated “true” effects obtained by applying smooth functions to the raw NHANES data.

- **`/src/`**  
  Source code, including both the original functions (from Lily’s code) and revised functions for running simulations and visualizations.

- **`/src_simulation/`**  
  Functions for running non-informative sampling simulations under different settings.

- **`/results/`**  
  Saved simulation results generated using the functions in `/src_simulation/`.

- **`/src_simulation_old/`**
  Early versions of simulation functions (not used currently).

- **`/bash_old/`**
  Old bash scripts corresponding to the early simulation functions (not used currently).

### Files
- **`CR_issue_0701.qmd`**  
  Quarto document that visualizes simulation results under different settings.  
  *Main finding:* The joint-wise coverage ratio varies considerably across settings.

- **`base_smoothed.rmd`**  
  Early R Markdown file visualizing how different smoothing functions perform on NHANES data, and how they affect simulation results in the non-informative sampling setting.

- **`empirical_0510.qmd`**  
  Early document that visualizes and summarizes empirical simulation results derived from NHANES data under different sampling settings.

