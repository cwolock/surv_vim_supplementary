#!/usr/local/bin/Rscript

## set up directories for output, plots
output_dir <- "output/"

## set up parameter grid
name <- "male_analysis"
nreps_total <- 50
nreps_per_job <- 1
approach <- c("conditional", "marginal")

## number of monte-carlo iterations per job
nreps_per_combo <- nreps_total/nreps_per_job
## set up grid of parameters

param_grid <- expand.grid(mc_id = 1:nreps_per_combo,
                          approach = approach)

## ----------------------------------------
## load the results
## ----------------------------------------
## names of files to read in
output_nms <- paste0(name, "_", 1:dim(param_grid)[1], ".rds")
avail_nms <- list.files(output_dir, pattern = paste0(name, "_*"))
names_to_try <- output_nms[which(output_nms %in% avail_nms)]
## list of output
output_lst <- lapply(paste0(output_dir, names_to_try), readRDS)
## make it a matrix
output_df <- do.call(rbind.data.frame, output_lst)

all_seeds <- 1:nrow(param_grid)
completed_seeds <- output_df$seed

missing_seeds <- all_seeds[!(all_seeds %in% completed_seeds)]

missing_settings <- param_grid[missing_seeds,]
missing_settings$seed <- missing_seeds

saveRDS(output_df, paste0(name, ".rds"))
saveRDS(missing_settings, paste0(name, "_missing.rds"))
