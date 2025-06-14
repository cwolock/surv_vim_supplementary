#!/usr/local/bin/Rscript

sim_name <- "new_survML_testing_rmst"
nreps_total <- 500
nreps_per_job <- 1

output_dir <- "output/"

n_trains <- c(500, 750, 1000, 1250, 1500)
crossfits <- c(TRUE, FALSE)

njobs_per_combo <- nreps_total/nreps_per_job

param_grid <- expand.grid(mc_id = 1:njobs_per_combo,
                          crossfit = crossfits,
                          n_train = n_trains)

output_nms <- paste0(sim_name, "_", 1:dim(param_grid)[1], ".rds")
avail_nms <- list.files(output_dir, pattern = paste0(sim_name, "_*"))
names_to_try <- output_nms[which(output_nms %in% avail_nms)]
missing <- output_nms[which(!(output_nms %in% avail_nms))]
print(missing)
output_lst <- lapply(paste0(output_dir, names_to_try), readRDS)
output_df <- do.call(rbind.data.frame, output_lst)
saveRDS(output_df, paste0(sim_name, ".rds"))
