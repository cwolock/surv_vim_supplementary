#!/usr/local/bin/Rscript

## load results, make nice plots

## ----------------------------------------
## load packages and user-defined functions
## ----------------------------------------
## nice plots
suppressMessages(library("ggplot2"))
suppressMessages(library("cowplot"))
## tidy stuff
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
## command line args (if specified)
suppressMessages(library("argparse"))
suppressMessages(library("grid"))
suppressMessages(library("gridExtra"))
suppressMessages(library("ggpubr"))
## ----------------------------------------
## read in command line args (if any)
## ----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--sim-name", default = "sim",
                    help = "name of simulation")
parser$add_argument("--nreps-total", type = "double", default = 1000,
                    help = "number of replicates for each set of params")
parser$add_argument("--nreps-per-job", type = "double", default = 100,
                    help = "number of replicates per job")
args <- parser$parse_args()

## set up directories for output, plots
output_dir <- "output/"
plots_dir <- "plots/"

## set up parameter grid
n_trains <- c(500, 750, 1000, 1250, 1500)
nuisances <- c("survSL")
crossfits <- c(FALSE, TRUE)
scenarios <- c("A", "C", "D")
## number of monte-carlo iterations per job
nreps_per_combo <- args$nreps_total/args$nreps_per_job
## set up grid of parameters
param_grid <- expand.grid(mc_id = 1:nreps_per_combo,
                          scenario = scenarios,
			  crossfit = crossfits,
                          n_train = n_trains,
                          nuisance = nuisances)

## ----------------------------------------
## load the results
## ----------------------------------------
## names of files to read in
output_nms <- paste0(args$sim_name, "_", 1:dim(param_grid)[1], ".rds")
avail_nms <- list.files(output_dir, pattern = paste0(args$sim_name, "_*"))
names_to_try <- output_nms[which(output_nms %in% avail_nms)]
print(length(output_nms) - length(avail_nms))
print(output_nms[which(!(output_nms %in% avail_nms))])
print(avail_nms[which(!(avail_nms %in% output_nms))])
## list of output
output_lst <- lapply(paste0(output_dir, names_to_try), readRDS)
## make it a matrix
output_df <- do.call(rbind.data.frame, output_lst)

all_seeds <- 1:nrow(param_grid)
completed_seeds <- output_df$seed

missing_seeds <- all_seeds[!(all_seeds %in% completed_seeds)]

missing_settings <- param_grid[missing_seeds,]
missing_settings$seed <- missing_seeds

saveRDS(output_df, paste0(args$sim_name, ".rds"))
saveRDS(missing_settings, paste0(args$sim_name, "_missing.rds"))
