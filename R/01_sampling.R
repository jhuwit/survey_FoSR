library(tidyverse)
raw = read_rds(here::here("data", "sim", "superpop_nopsu.rds"))
raw_psu = read_rds(here::here("data", "sim", "superpop_psu.rds"))

get_sample = function(num_strata = 30, goal_n_per_psu_strata = 50, pop_data, seed = 123, sp = 0){
  I = nrow(pop_data)
  set.seed(seed)
  dirichlet_probs = gtools::rdirichlet(1, rep(1, num_strata))
  set.seed(seed)
  stratum_assignments = sample(1:num_strata, size = I, replace = TRUE, prob = dirichlet_probs)

  pop_data =
    pop_data %>%
    mutate(strata = stratum_assignments)

  Y_mat = unclass(pop_data$Y)
  fpca_Y = refund::fpca.face(Y_mat) ### plot

  q31 = quantile(fpca_Y$scores[,1], 0.75) %>% unname()
  q32 = quantile(fpca_Y$scores[,2], 0.75) %>% unname()

  # pop_data =
  #   pop_data %>%
  #   rowwise() %>%
  #   mutate(Y_mean = mean(Y)) %>%
  #   ungroup()
  #
  # q3 = quantile(pop_data$Y_mean, 0.75) %>% unname()

  # this selection prob thing doesn't work because it's just shifting the mean, not the shape
  # base on prinicpal components


  pop_data =
    pop_data %>%
    mutate(fpc1_score = fpca_Y$scores[,1],
           fpc2_score = fpca_Y$scores[,2]) %>%
    mutate(selection_prob =
             # sp * (fpc1_score > q31) +  0.5 * (PSU == 1) + 0.5 * (PSU == 2)) %>%
            (sp/2 * (fpc1_score > q31)) + (sp/2 * fpc2_score > q32)  + 0.5 * (PSU == 1) + 0.5 * (PSU == 2)) %>%
    group_by(strata) %>%
    mutate(selection_prob = selection_prob / sum(selection_prob),
           weight = 1 / selection_prob) %>%
    ungroup()

  goal_n = goal_n_per_psu_strata * num_strata * 2
  c = goal_n / sum(pop_data$selection_prob)

  # Adjust probabilities
  pop_data$adjusted_prob = pop_data$selection_prob * c
  # Ensure adjusted probabilities sum to 1
  pop_data$adjusted_prob = pop_data$adjusted_prob / sum(pop_data$adjusted_prob)

  set.seed(seed)
  final_sample = pop_data %>%
    slice_sample(n = goal_n, weight_by = adjusted_prob, replace = FALSE)

 return(final_sample)
}

sample_psu = get_sample(pop_data = raw_psu, sp = 2.5)
sample = get_sample(30, 10000, raw, sp = 2.5)
sample_psu_noprob = get_sample(30, 10000, raw_psu)
sample_noprob = get_sample(30, 10000, raw)

write_rds(sample_psu, here::here("data", "sim", "sample_psu_sp.rds"))
write_rds(sample, here::here("data", "sim", "sample_sp.rds"))
write_rds(sample_psu_noprob, here::here("data", "sim", "sample_psu.rds"))
write_rds(sample_noprob, here::here("data", "sim", "sample.rds"))




############ old




##### what I understand from Leroux implementation
## generate strata assignment
set.seed(123) # For reproducibility
num_strata = 30
I = nrow(raw)
dirichlet_probs = gtools::rdirichlet(1, rep(1, num_strata)) # probabilities of assignment to each strata, sum to 1

stratum_assignments = sample(1:num_strata, size = I, replace = TRUE, prob = dirichlet_probs)
table(stratum_assignments) # check how many ppl assigned to each stratum

raw =
  raw %>%
  mutate(strata = stratum_assignments)

raw %>%
  group_by(strata, PSU) %>%
  count()

means = raw %>%
  rowwise() %>%
  mutate(Y_mean = mean(Y))

# summary(means$Y_mean)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -0.70180 -0.21991 -0.15196 -0.15203 -0.08404  0.34536
#

# within ea strata, generate selection probabilities
# we want the selection prob to be associated with the outcome
# use mean function of outcome for now but consider also functional PC scores for future
raw =
  raw %>%
  rowwise() %>%
  mutate(Y_mean = mean(Y)) %>%
  ungroup() %>%
  mutate(selection_prob =
           2.5 * (Y_mean > -0.08) + 0.5 * (PSU == 1) + 0.5 * (PSU == 2)) %>%
  group_by(strata) %>%
  mutate(selection_prob = selection_prob / sum(selection_prob),
         weight = 1 / selection_prob) %>%
  ungroup()


goal_n = 10000
c = goal_n / sum(raw$selection_prob)

# Adjust probabilities
raw$adjusted_prob = raw$selection_prob * c
# Ensure adjusted probabilities sum to 1
raw$adjusted_prob = raw$adjusted_prob / sum(raw$adjusted_prob)


# Sample individuals without replacement
set.seed(123) # For reproducibility
final_sample = raw %>%
  slice_sample(n = goal_n, weight_by = adjusted_prob, replace = FALSE)

final_sample %>%
  group_by(strata, PSU) %>%
  count()



write_rds(final_sample,
          here::here("data", "sim", "final_sample_leroux2.rds"),
          compress = "xz")




