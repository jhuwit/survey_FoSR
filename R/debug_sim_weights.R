library(future)
library(furrr)
library(tidyverse)
library(here)
library(devtools)
library(fastFMM)
library(haven)
library(dplyr)
library(survey)
library(progress)
library(lme4)
library(paletteer)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(tidyfun)
library(mvtnorm)
library(refund)
library(patchwork)
library(svrep)
source(here::here("R", "01_sim_functions.R"))
source(here::here("R", "00_data_gen_function_ff.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "create_settings.R"))

force = FALSE
force_iter = FALSE
nsim = 200
B = 500
ncores = parallelly::availableCores() - 1
fit_types = c('weighted', 'unweighted')
boot_types = c('Rao-Wu-Yue-Beaumont', 'BRR', 'weighted', 'unweighted')
options(survey.lonely.psu = "adjust")


ifold = get_fold()

temp = settings %>%
  filter(fold == ifold)

outfile = here::here("debug", "weight_hist", paste0("fold_", sprintf("%03d", ifold), ".pdf"))
outfile_dist = here::here("debug", "weight_summ", paste0("fold_", sprintf("%03d", ifold), ".rds"))
# create partial dir to store ea. iteration

if (!dir.exists(dirname(outfile))) dir.create(dirname(outfile), recursive = TRUE)


if(!file.exists(outfile) || force) {
  temp = settings[ifold, ] # set settings for this simulation

  lst = generate_superpopulation(
    scenario = temp$scenario,
    family = temp$family,
    I = 10e6,
    L = temp$len,
    snr_b = temp$snr_b,
    snr_eps = temp$snr_eps,
    strata_sigma = temp$strata_sigma,
    strata_scale = temp$strata_scale,
    seed = 111
  )

  weights = list()

  pdf(here::here(outfile))
  for (iter in 1:nsim) {

    set.seed(iter)
    data = sample_from_population_wor(
        X_des = lst$X_des,
        Y_obs = lst$Y_obs,
        L = temp$len,
        I_n = temp$In,
        num_strata = 30,
        stratum_assignments = lst$stratum_assignments,
        psu_assignments = lst$psu_assignments,
        dirichlet_probs = lst$dirichlet_probs,
        seed = iter,
        inf_level = temp$inf_level
      )
      weights[[iter]] = summary(data$weight) %>% unname()
      summary = summary(data$weight)
      p1 = data %>%
        ggplot(aes(x = weight)) +
        geom_histogram(col = "black")

      p2 =
        data %>%
        mutate(wt_norm = weight / mean(weight)) %>%
        ggplot(aes(x = wt_norm)) +
        geom_histogram(col = "black")

      print(p1 / p2)
      rm(data)
  }
  dev.off()
  wt_df = bind_rows(weights) %>%
    magrittr::set_colnames(c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max"))


  write_rds(wt_df, file = outfile_dist)


}


