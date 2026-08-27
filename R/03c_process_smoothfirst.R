# library(tidyverse)
library(readr)
library(dplyr)
library(purrr)
# library(kableExtra)
# library(patchwork)
source(here::here("R", "create_survey_settings.R"))

files = list.files(here::here("results", "simulations", "survey_sim_smoothfirst"), full.names = TRUE,
                   pattern = ".rds")


nums = sub(".*fold\\_(.+).rds.*", "\\1", basename(files)) %>% as.numeric

length(files)
process_file = function(x){
  f = read_rds(x)
  if(nrow(f) > 0) {
    f = read_rds(x) %>%
      mutate(fold = sub(".*fold\\_(\\d+)\\.rds.*", "\\1", x),
             fold = as.integer(fold)) %>%
      select(id, MISE, starts_with("mean"), starts_with("cover"), starts_with("var"), boot_type, n, n_boot, fold, everything()) %>%
      left_join(settings, by = join_by(fold))

    if(f$len[1] == 1440){
      if (nrow(f) >= 100 * 4 * 2){
        f = f %>%
          mutate(id2 = as.numeric(id)) %>%
          filter(id2 <= 100) %>%
          select(-id2)
        return(f)
      } else return(NULL)
    } else {
      if (nrow(f) >= 200 * 4 * 2){
        return(f)
      } else return(NULL)
    }
  } else return(NULL)
}
all_res = map(.x = files, process_file) %>% bind_rows()

all_res = all_res %>%
  mutate(rn = row_number()) %>%
  mutate(var = case_when(var %in% c("x", "X") ~ "x",
                         var %in% c("Intercept", "(Intercept)") ~ "Intercept",
                         is.na(var) & rn %%2 == 1 ~ "x",
                         is.na(var) & rn %% 2 == 0 ~ "Intercept",
                         .default = NA_character_))

all_res = all_res %>%
  filter(!(var == "Intercept" & len == 1440),
         !(len == 1440 & family != "gaussian"))


write_rds(all_res, here::here("results", "simulations", "all_survey_res_smoothfirst.rds"))
