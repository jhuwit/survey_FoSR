get_fold = function() {
  ifold = as.numeric(Sys.getenv("SGE_TASK_ID"))
  if (is.na(ifold)) {
    ifold = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  }
  print(paste0("fold is: ", ifold))
  ifold
}
