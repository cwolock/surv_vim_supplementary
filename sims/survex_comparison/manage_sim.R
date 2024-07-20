#!/usr/local/bin/Rscript
R.Version()
.libPaths()
#library(Rsolnp)
#library(mgcv)
library(dplyr)
library(survML)
library(SuperLearner)
library(survival)
library(randomForestSRC)
library(survSuperLearner)
library(survex)
library(tidyr)
source("/home/cwolock/surv_vim_supplementary/sims/survex_comparison/do_one.R")
source("/home/cwolock/surv_vim_supplementary/sims/utils.R")
source("/home/cwolock/surv_vim_supplementary/sims/generate_data.R")
# source("/home/cwolock/surv_vim_supplementary/sims/survSL_wrappers.R")
#source("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/survex_comparison/do_one.R")
#source("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/utils.R")
#source("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/generate_data.R")
#source("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/survSL_wrappers.R")

sim_name <- "survex_comparison_brier"
nreps_total <- 500
nreps_per_job <- 1

n_trains <- c(500, 1000,1500, 2000, 2500, 3000)
methods <- c("permutation", "exclusion")
scenarios <- c("5A", "5B", "5C", "5D")

njobs_per_combo <- nreps_total/nreps_per_job

param_grid <- expand.grid(mc_id = 1:njobs_per_combo,
                          n_train = n_trains,
                          method = methods,
                          scenario = scenarios)

job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

current_dynamic_args <- param_grid[job_id, ]
set.seed(1234) # overall seed
current_seed <- as.integer((1e9*runif(job_id))[job_id])
set.seed(current_seed)
output <- replicate(nreps_per_job,
                    do_one(n_train = current_dynamic_args$n_train,
                          method = current_dynamic_args$method,
                          scenario = current_dynamic_args$scenario),
                    simplify = FALSE)
sim_output <- lapply(as.list(1:length(output)),
                     function(x) tibble::add_column(output[[x]]))
sim_output_tib <- do.call(rbind.data.frame, sim_output)
file_name <- paste0("output/", sim_name, "_", job_id, ".rds")
saveRDS(sim_output_tib, file = file_name)
