#!/usr/local/bin/Rscript
library(dplyr)
library(survML)
library(SuperLearner)
library(survival)
source("/home/cwolock/surv_vim_supplementary/sims/new_survML_testing/rmst/do_one.R")
source("/home/cwolock/surv_vim_supplementary/sims/utils.R")
source("/home/cwolock/surv_vim_supplementary/sims/generate_data.R")
source("/home/cwolock/surv_vim_supplementary/sims/survSL_wrappers.R")

sim_name <- "new_survML_testing_rmst"
nreps_total <- 500
nreps_per_job <- 1

n_trains <- c(500, 750, 1000, 1250, 1500)
crossfits <- c(FALSE, TRUE)

njobs_per_combo <- nreps_total/nreps_per_job

param_grid <- expand.grid(mc_id = 1:njobs_per_combo,
                          crossfit = crossfits,
                          n_train = n_trains)

job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

current_dynamic_args <- param_grid[job_id, ]
set.seed(1234) # overall seed
current_seed <- as.integer((1e9*runif(job_id))[job_id])
set.seed(current_seed)
output <- replicate(nreps_per_job,
                    do_one(n_train = current_dynamic_args$n_train,
                           crossfit = current_dynamic_args$crossfit),
                    simplify = FALSE)
sim_output <- lapply(as.list(1:length(output)),
                     function(x) tibble::add_column(output[[x]]))
sim_output_tib <- do.call(rbind.data.frame, sim_output)
file_name <- paste0("output/", sim_name, "_", job_id, ".rds")
saveRDS(sim_output_tib, file = file_name)
