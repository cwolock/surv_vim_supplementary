#!/usr/local/bin/Rscript
library(dplyr)
library(survival)
library(survML)
library(SuperLearner)

source("/home/cwolock/surv_vim_supplementary/sims/utils.R")
source("/home/cwolock/surv_vim_supplementary/sims/boost_c_index.R")
source("/home/cwolock/surv_vim_supplementary/data_analysis/male/702_data_analysis.R")

name <- "combined_analysis"
nreps_total <- 10
nreps_per_job <- 1

approach <- c("conditional","marginal")

njobs_per_combo <- nreps_total/nreps_per_job

param_grid <- expand.grid(mc_id = 1:njobs_per_combo,
                          approach = approach)

job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

current_dynamic_args <- param_grid[job_id, ]
global_seed <- 92723
set.seed(global_seed) # overall seed
current_seed <- as.integer((1e9*runif(job_id))[job_id])
set.seed(current_seed)
output <- replicate(nreps_per_job,
                    do_one(approach = current_dynamic_args$approach,
                           seed = current_seed,
                           global_seed = global_seed),
                    simplify = FALSE)
sim_output <- lapply(as.list(1:length(output)),
                     function(x) tibble::add_column(output[[x]]))
sim_output_tib <- do.call(rbind.data.frame, sim_output)
file_name <- paste0("output/", name, "_", job_id, ".rds")
saveRDS(sim_output_tib, file = file_name)
