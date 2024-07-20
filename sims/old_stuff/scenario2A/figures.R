source("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/figure_utils_newA.R")

# read in truth files
# truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim/scratch/sims/landmark/truth.rds"
truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/truth_interaction.rds"
var_truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/mc_truth_variance_interaction.rds"
truth_list <- compile_truth(true_param_file = truth_file,
                            true_avar_file = var_truth_file)

# LANDMARK SIMS
# landmark_dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario1_landmark.rds")
landmark_dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario2A_landmark.rds")
landmark_summ <- summarize_results(landmark_dat, scenario = "2", truth_list$truth, truth_list$var_truth)

# C-INDEX SIMS
cindex_dat <- readRDS("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/scenario2A_cindex.rds")
# cindex_dat <- cindex_dat %>% mutate(est = one_step)
cindex_dat <- cindex_dat %>% mutate(tau = restriction_time)
cindex_summ <- summarize_results(cindex_dat, scenario = "2", truth_list$truth, truth_list$var_truth)

summ <- bind_rows(landmark_summ, cindex_summ)
make_sim_plot(summ,
              scenario = "2",
              big = TRUE,
              wd = "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/biometrika/",
              fname = "scenario2A-big_070124")

# summ <- summ %>% filter(vim == "AUC" & indx == 1)
# make_sim_plot(summ,
#               scenario = "2",
#               big = FALSE,
#               wd = "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/biometrika/",
#               fname = "scenario2A-small")
