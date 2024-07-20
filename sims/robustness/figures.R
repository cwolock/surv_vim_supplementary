source("/Users/cwolock//Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/sims/figure_utils.R")

# read in truth files
truth_file <- "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/truth_interactionB.rds"
var_truth_file <- "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/mc_truth_variance_interactionB.rds"
truth_list <- compile_truth(true_param_file = truth_file,
                            true_avar_file = var_truth_file)

dat <- readRDS("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/scenario5B.rds")
dat <- dat %>% filter(n_train %in% c(250, 500, 1000, 2500, 5000)) %>%
  mutate(tau = landmark_time)

dat <- dat %>% mutate(correlation = TRUE)
dat <- left_join(dat, truth_list$truth, by = c("tau", "vim", "correlation"))

dat <- dat %>% mutate(param = case_when(
  indx == 1 ~ V_full - V_01,
  indx == 2 ~ V_full - V_02,
  indx == 6 ~ V_full - V_06,
  indx == "1,6" ~ V_full - V_016
))

dat <- dat %>% mutate(err = (est - param))

this_vim <- "AUC"
this_indx <- "1,6"
this_t <- 0.5#

this_dat <- dat %>% filter(vim == this_vim & tau == this_t & indx == this_indx)
summ <- this_dat %>%
  group_by(misspec_type, n_train, indx, robust_V, robust_f) %>%
  summarize(nreps = n(),
            mean_err = mean(err),
            mse = mean(err^2),
            sd_err = sd(err),
            se_err = sqrt( mean((err - mean_err) ^ 2) / (nreps-1)),
            runtime = mean(runtime)) %>%
  mutate(misspec_type = factor(misspec_type,
                               levels = c("none", "censoring", "event_plusf0"),
                               labels = c("None", "Censoring", "Event")),
         robust_V = factor(robust_V, levels = c(FALSE, TRUE),
                           labels = c("Direct", "Indirect")),
         robust_f = factor(robust_f, levels = c(FALSE, TRUE),
                           labels = c("Conditional surv. function", "Doubly-robust pseudo-outcome")))

p <- summ %>% ggplot(aes(x = factor(n_train), y = mse)) +
  geom_point(aes(shape = robust_V), size = 2.5, color = "grey40") +
  geom_line(aes(linetype = robust_f,
                group = interaction(robust_f, robust_V)),
            color = "grey40") +
  scale_color_manual(values = c("red", "blue")) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  scale_shape_manual(values = c(0,8)) +
  facet_wrap(~ misspec_type) +
  theme_bw() +
  geom_hline(yintercept = 0) +
  ggtitle("Misspecification type") +
  xlab("Sample size") +
  ylab("Empirical mean squared error") +
  theme( strip.background = element_blank(),
         strip.placement = "outside",
         plot.title = element_text(hjust = 0.5, size = 16,family = "Times New Roman"),
         legend.title = element_text(size = 16,family = "Times New Roman"),
         legend.text = element_text(size = 14,family = "Times New Roman"),
         legend.position="bottom",
         panel.spacing.x = unit(0.3, "cm"),
         axis.text = element_text(size = 14,family = "Times New Roman"),
         strip.text = element_text(size = 14, family = "Times New Roman"),
         axis.title = element_text(size = 14, family = "Times New Roman"),
         panel.grid.major.x = element_blank(),
         panel.grid.minor = element_blank()) +
  guides(linetype=guide_legend(title="Oracle estimator:", nrow = 1, ncol = 2),
         shape = guide_legend(title = "Debiasing:", nrow = 1, ncol = 2))

ggsave(filename = "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/biometrika/robust_mseB_071724.pdf",
       plot = p, device = "pdf",
       width = 11, height = 5, dpi = 300, units = "in")