source("C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/figure_utils.R")

# read in truth files
truth_file <- "C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim/scratch/sims/landmark/truth.rds"
var_truth_file <- "C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/variance_truth.rds"
truth_list <- compile_truth(true_param_file = truth_file,
                            true_avar_file = var_truth_file)

# LANDMARK SIMS
landmark_dat <- readRDS("C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario4_landmark.rds")
landmark_summ <- summarize_results(landmark_dat, scenario = "4", truth_list$truth, truth_list$var_truth)

# C-INDEX SIMS
cindex_dat <- readRDS("C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario4_cindex.rds")
cindex_summ <- summarize_results(cindex_dat, scenario = "4", truth_list$truth, truth_list$var_truth)


summ <- bind_rows(landmark_summ, cindex_summ)
make_sim_plot(summ,
              scenario = "4",
              big = TRUE,
              wd = "C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/",
              fname = "scenario2_big")
