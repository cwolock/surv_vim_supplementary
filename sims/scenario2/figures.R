source("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/sims/figure_utils.R")

# read in truth files
truth_file <- "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/truth_interactionB.rds"
var_truth_file <- "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/mc_truth_variance_interactionB.rds"
truth_list <- compile_truth(true_param_file = truth_file,
                            true_avar_file = var_truth_file)

# LANDMARK SIMS
landmark_dat <- readRDS("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/scenario2B_landmark.rds")
landmark_summ <- summarize_results(landmark_dat, scenario = "2", truth_list$truth, truth_list$var_truth)

# C-INDEX SIMS
cindex_dat <- readRDS("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/scenario2B_cindex.rds")
cindex_dat <- cindex_dat %>% mutate(tau = restriction_time)
cindex_summ <- summarize_results(cindex_dat, scenario = "2", truth_list$truth, truth_list$var_truth)

summ <- bind_rows(landmark_summ, cindex_summ)
make_sim_plot(summ,
              scenario = "2",
              big = TRUE,
              wd = "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/biometrika/",
              fname = "scenario2-big")
