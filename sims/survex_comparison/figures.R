source("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/figure_utils_new.R")

# read in truth files
# truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim/scratch/sims/landmark/truth.rds"
truth <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/truth_permutation.rds")
truth <- truth %>% mutate(vim1 = V_full - V_01,
                          vim2 = V_full - V_02,
                          vim3 = V_full - V_03,
                          vim4 = V_full - V_04,
                          landmark_time = factor(tau))

dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/survex_comparison.rds")

summ <- dat %>% group_by(n_train, indx, method, landmark_time) %>%
  summarize(rank_right = mean(rank == true_rank),
            all_right = mean(all_right),
            avg_est = mean(est),
            runtime = mean(runtime)) %>%
  mutate(landmark_time = factor(landmark_time))

summ_small <- summ %>% filter(indx == "1")
p <- summ_small %>% ggplot(aes(x = n_train, y = all_right)) +
  geom_point(aes(shape = method, color = landmark_time)) +
  theme_bw() +
  scale_color_manual(values = c("black", "blue")) +
  ylim(c(0, 1))


# dat <- left_join(dat, truth, by = c("landmark_time"))
# dat <- dat %>% mutate(truth = case_when(indx == "1" ~ vim1,
#                                         indx == "2" ~ vim2,
#                                         indx == "3" ~ vim3,
#                                         indx == "4" ~ vim4,
#                                         indx == "5" ~ vim5))
# summ <- dat %>% group_by(n_train, indx, method, landmark_time) %>%
#   summarize(mean(est),
#             mean(runtime))
#
# # LANDMARK SIMS
# landmark_dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario1_landmark.rds")
# landmark_dat <- landmark_dat %>% mutate(est = one_step)
# landmark_dat <- landmark_dat %>% mutate(landmark_time = tau)
# # landmark_dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario1_landmark_053024.rds")
# landmark_summ <- summarize_results(landmark_dat, scenario = "1", truth_list$truth, truth_list$var_truth)
#
# # C-INDEX SIMS
# cindex_dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario1_cindex.rds")
# cindex_dat <- cindex_dat %>% mutate(est = one_step)
# # cindex_dat <- cindex_dat %>% mutate(tau = restriction_time)
# cindex_summ <- summarize_results(cindex_dat, scenario = "1", truth_list$truth, truth_list$var_truth)
#
# summ <- bind_rows(landmark_summ, cindex_summ)
# make_sim_plot(summ,
#               scenario = "1",
#               big = TRUE,
#               wd = "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/biometrika/",
#               fname = "scenario1-big_061024")

# summ <- summ %>% filter(vim == "AUC" & indx == 1)
# make_sim_plot(summ,
#               scenario = "1",
#               big = FALSE,
#               wd = "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/biometrika/",
#               fname = "scenario1-small")
