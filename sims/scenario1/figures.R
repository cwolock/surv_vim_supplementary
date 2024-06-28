source("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/figure_utils_new.R")

# read in truth files
# truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim/scratch/sims/landmark/truth.rds"
truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/truth.rds"
var_truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/variance_truth.rds"
truth_list <- compile_truth(true_param_file = truth_file,
                            true_avar_file = var_truth_file)

# LANDMARK SIMS
landmark_dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario1_landmark.rds")
landmark_dat <- landmark_dat %>% mutate(landmark_time = tau)
# landmark_dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario1_landmark_053024.rds")
landmark_summ <- summarize_results(landmark_dat, scenario = "1", truth_list$truth, truth_list$var_truth)

# C-INDEX SIMS
cindex_dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario1_cindex_061024_subsamp.rds")
# cindex_dat <- cindex_dat %>% mutate(est = one_step)
cindex_dat <- cindex_dat %>% mutate(tau = restriction_time)
cindex_summ <- summarize_results(cindex_dat, scenario = "1", truth_list$truth, truth_list$var_truth)

summ <- bind_rows(landmark_summ, cindex_summ)
make_sim_plot(summ,
              scenario = "1",
              big = TRUE,
              wd = "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/biometrika/",
              fname = "scenario1-big_061024")

# summ <- summ %>% filter(vim == "AUC" & indx == 1)
# make_sim_plot(summ,
#               scenario = "1",
#               big = FALSE,
#               wd = "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/biometrika/",
#               fname = "scenario1-small")
