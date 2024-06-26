#!/usr/local/bin/Rscript
.libPaths(c(
  "/home/cwolock/R_lib",
  .libPaths()
))
suppressMessages(library(dplyr))
suppressMessages(library(survML))
suppressMessages(library(survival))
suppressMessages(library(survSuperLearner))

source("/home/cwolock/surv_vim_supplementary/sims/scenario5/do_one_bothrobust.R")
source("/home/cwolock/surv_vim_supplementary/sims/utils.R")
source("/home/cwolock/surv_vim_supplementary/sims/generate_data.R")

sim_name <- "scenario5_bothrobust_moreNs"
nreps_total <- 500
nreps_per_job <- 1

n_trains <- c(250, 500, 2500)
misspec_types <- c("none", "censoring", "event_plusf0")
robust_Vs <- c(FALSE, TRUE)
robust_fs <- c(FALSE, TRUE)

njobs_per_combo <- nreps_total/nreps_per_job

param_grid <- expand.grid(mc_id = 1:njobs_per_combo,
                          n_train = n_trains,
                          misspec_type = misspec_types,
                          robust_f = robust_fs,
			  robust_V = robust_Vs)

job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

current_dynamic_args <- param_grid[job_id, ]
set.seed(1234) # overall seed
current_seed <- as.integer((1e9*runif(job_id))[job_id])
set.seed(current_seed)
output <- replicate(nreps_per_job,
                    do_one(n_train = current_dynamic_args$n_train,
                           misspec_type = current_dynamic_args$misspec_type,
                           robust_f = current_dynamic_args$robust_f,
			   robust_V = current_dynamic_args$robust_V),
                    simplify = FALSE)
sim_output <- lapply(as.list(1:length(output)),
                     function(x) tibble::add_column(output[[x]]))
sim_output_tib <- do.call(rbind.data.frame, sim_output)
file_name <- paste0("output/", sim_name, "_", job_id, ".rds")
saveRDS(sim_output_tib, file = file_name)
