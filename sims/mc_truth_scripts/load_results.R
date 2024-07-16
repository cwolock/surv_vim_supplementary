#!/usr/local/bin/Rscript

sim_name <- "mc_truth_variance_interactionB"
nreps_total <- 100
nreps_per_job <- 1

n_trains <- c(1000)
correlations <- c("correlation", "no_correlation")

output_dir <- "output/"

nreps_per_combo <- nreps_total/nreps_per_job
param_grid <- expand.grid(mc_id = 1:nreps_per_combo,
                          n_train = n_trains,
                          correlation = correlations)

output_nms <- paste0(sim_name, "_", 1:dim(param_grid)[1], ".rds")
avail_nms <- list.files(output_dir, pattern = paste0(sim_name, "_*"))
names_to_try <- output_nms[which(output_nms %in% avail_nms)]
output_lst <- lapply(paste0(output_dir, names_to_try), readRDS)
output_df <- do.call(rbind.data.frame, output_lst)
saveRDS(output_df, paste0(sim_name, ".rds"))
