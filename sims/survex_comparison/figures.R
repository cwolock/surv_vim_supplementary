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
# dat$correlation <- c(rep(TRUE, nrow(dat)/2), rep(FALSE, nrow(dat)/2))
dat <- dat %>% filter(landmark_time == 0.5) %>%
  mutate(scaled_est = ifelse(est <= 0, 0, est))
# dat$rep <- rep(seq(1:(nrow(dat)/4)), each = 4)
#
# dat$all_right <- NA
# for (i in 1:(nrow(dat)/4)){
#   this_dat <- dat %>% filter(rep == i)
#   if (this_dat$correlation[1]){
#     all_right <- all(this_dat$rank == c(2,1,3,4))
#   } else{
#     all_right <- all(this_dat$rank == c(1,2,3,4))
#   }
#   dat[dat$rep == i,"all_right"] <- all_right
# }

summ <- dat %>% group_by(n_train, method,scenario) %>%
  summarize(nreps = n(),
            rank_right = mean(correct),
            mse = mean(est^2),
            rank_right_mc_se = sqrt(rank_right * (1 - rank_right) / nreps),
            # all_right = mean(all_right),
            avg_est = mean(scaled_est),
            runtime = mean(runtime))# %>%
  # mutate(landmark_time = factor(landmark_time))

summ_small <- summ #%>% filter(indx == "1")
p <- summ_small %>% ggplot(aes(x = n_train, y = rank_right)) +
  geom_line(aes(linetype = scenario)) +
  geom_errorbar(aes(ymin=rank_right - 1.96*rank_right_mc_se,
                    ymax=rank_right + 1.96*rank_right_mc_se),
                width=.1) +
  geom_point(size = 1) +
  facet_wrap(~ method) +
  theme_bw() +
  scale_color_manual(values = c("black", "blue")) +
  ylim(c(0, 1))

p <- summ_small %>% ggplot(aes(x = n_train, y = avg_est)) +
  geom_line(aes(linetype = scenario)) +
  # geom_errorbar(aes(ymin=rank_right - 1.96*rank_right_mc_se,
                    # ymax=rank_right + 1C.96*rank_right_mc_se),
                # width=.1) +
  geom_point(size = 1) +
  facet_wrap(~ method) +
  theme_bw() +
  scale_color_manual(values = c("black", "blue"))# +
  # ylim(c(0, 1))

p <- summ_small %>% ggplot(aes(x = n_train, y = mse)) +
  geom_line(aes(linetype = scenario)) +
  # geom_errorbar(aes(ymin=rank_right - 1.96*rank_right_mc_se,
  # ymax=rank_right + 1C.96*rank_right_mc_se),
  # width=.1) +
  geom_point(size = 1) +
  facet_wrap(~ method) +
  theme_bw() +
  scale_color_manual(values = c("black", "blue"))# +
# ylim(c(0, 1))


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
