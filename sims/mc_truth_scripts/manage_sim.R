#!/usr/local/bin/Rscript
library(dplyr)

source("/home/cwolock/surv_vim_supplementary/sims/mc_truth_scripts/do_one_interactionB.R")
source("/home/cwolock/surv_vim_supplementary/sims/generate_data.R")

sim_name <- "mc_truth_variance_interactionB"
nreps_total <- 100
nreps_per_job <- 1

n_trains <- c(1000)
correlations <- c(TRUE, FALSE)

njobs_per_combo <- nreps_total/nreps_per_job

param_grid <- expand.grid(mc_id = 1:njobs_per_combo,
                          n_train = n_trains,
                          correlation = correlations)

job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

current_dynamic_args <- param_grid[job_id, ]

set.seed(1234) # overall seed
current_seed <- as.integer((1e9*runif(job_id))[job_id])
set.seed(current_seed)
output <- replicate(nreps_per_job,
                    do_one(n_train = current_dynamic_args$n_train,
                           correlation = current_dynamic_args$correlation),
                    simplify = FALSE)
sim_output <- lapply(as.list(1:length(output)),
                     function(x) tibble::add_column(output[[x]]))
sim_output_tib <- do.call(rbind.data.frame, sim_output)
file_name <- paste0("output/", sim_name, "_", job_id, ".rds")
saveRDS(sim_output_tib, file = file_name)
