#!/usr/local/bin/Rscript
library(dplyr)
library(survML)
library(SuperLearner)
library(survival)
library(randomForestSRC)
library(survSuperLearner)

source("/home/cwolock/surv_vim_supplementary/sims/scenario4/landmark/do_one.R")
source("/home/cwolock/surv_vim_supplementary/sims/utils.R")
source("/home/cwolock/surv_vim_supplementary/sims/generate_data.R")
source("/home/cwolock/surv_vim_supplementary/sims/survSL_wrappers.R")
sim_name <- "scenario4_landmark"
nreps_total <- 500
nreps_per_job <- 1

cens_rates <- c("30%", "40%", "50%", "60%", "70%")
nuisances <- c("survSL", "stackG", "rfsrc")
crossfits <- c(TRUE, FALSE)

njobs_per_combo <- nreps_total/nreps_per_job

param_grid <- expand.grid(mc_id = 1:njobs_per_combo,
                          crossfit = crossfits,
                          cens_rate = cens_rates,
                          nuisance = nuisances)

job_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

current_dynamic_args <- param_grid[job_id, ]
set.seed(1234) # overall seed
current_seed <- as.integer((1e9*runif(job_id))[job_id])
set.seed(current_seed)
output <- replicate(nreps_per_job,
                    do_one(cens_rate = current_dynamic_args$cens_rate,
                           nuisance = current_dynamic_args$nuisance,
                           crossfit = current_dynamic_args$crossfit),
                    simplify = FALSE)
sim_output <- lapply(as.list(1:length(output)),
                     function(x) tibble::add_column(output[[x]]))
sim_output_tib <- do.call(rbind.data.frame, sim_output)
file_name <- paste0("output/", sim_name, "_", job_id, ".rds")
saveRDS(sim_output_tib, file = file_name)
