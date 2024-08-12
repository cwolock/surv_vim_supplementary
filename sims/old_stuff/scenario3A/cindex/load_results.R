#!/usr/local/bin/Rscript

sim_name <- "scenario3A_cindex"
nreps_total <- 250
nreps_per_job <- 1

output_dir <- "output/"

cens_rates <- c("30%", "40%", "50%", "60%", "70%")
nuisances <- c("survSL", "stackG", "rfsrc")
crossfits <- c(FALSE, TRUE)

njobs_per_combo <- nreps_total/nreps_per_job

param_grid <- expand.grid(mc_id = 1:njobs_per_combo,
                          crossfit = crossfits,
                          cens_rate = cens_rates,
                          nuisance = nuisances)

output_nms <- paste0(sim_name, "_", 1:dim(param_grid)[1], ".rds")
avail_nms <- list.files(output_dir, pattern = paste0(sim_name, "_*"))
names_to_try <- output_nms[which(output_nms %in% avail_nms)]
output_lst <- lapply(paste0(output_dir, names_to_try), readRDS)
output_df <- do.call(rbind.data.frame, output_lst)
saveRDS(output_df, paste0(sim_name, ".rds"))
