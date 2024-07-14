source("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/sims/figure_utils_newB.R")

# read in truth files
# truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim/scratch/sims/landmark/truth.rds"
truth_file <- "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/truth_interactionB.rds"
var_truth_file <- "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/mc_truth_variance_interaction.rds"
truth_list <- compile_truth(true_param_file = truth_file,
                            true_avar_file = var_truth_file)

# LANDMARK SIMS
landmark_dat <- readRDS("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/scenario1B_landmark.rds")
landmark_dat <- landmark_dat %>% mutate(tau = landmark_time)
landmark_summ <- summarize_results(landmark_dat, scenario = "1", truth_list$truth, truth_list$var_truth)

# C-INDEX SIMS
cindex_dat <- readRDS("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/scenario1A_cindex.rds")
# cindex_dat <- cindex_dat %>% mutate(est = one_step)
cindex_dat <- cindex_dat %>% mutate(tau = restriction_time)
cindex_summ <- summarize_results(cindex_dat, scenario = "1", truth_list$truth, truth_list$var_truth)

summ <- bind_rows(landmark_summ, cindex_summ)
make_sim_plot(summ,
              scenario = "1",
              big = TRUE,
              wd = "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/biometrika/",
              fname = "scenario1B-big-071324")

# summ <- summ %>% filter(vim == "AUC" & indx == 1)
# make_sim_plot(summ,
#               scenario = "1",
#               big = FALSE,
#               wd = "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/biometrika/",
#               fname = "scenario1-small")
