source("/Users/cwolock//Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/figure_utils_newA.R")

# read in truth files
truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/truth_interaction.rds"
var_truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/mc_truth_variance_interaction.rds"
truth_list <- compile_truth(true_param_file = truth_file,
                            true_avar_file = var_truth_file)

dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario5A.rds")
dat <- dat %>% filter(n_train %in% c(250, 500, 1000, 2500, 5000)) %>%
  mutate(tau = landmark_time)

dat <- dat %>% mutate(correlation = TRUE)
dat <- left_join(dat, truth_list$truth, by = c("tau", "vim", "correlation"))

dat <- dat %>% mutate(param = case_when(
  indx == 1 ~ V_full - V_01,
  indx == 2 ~ V_full - V_02,
  indx == 5 ~ V_full - V_05,
  indx == "1,5" ~ V_full - V_015
))

dat <- dat %>% mutate(err = (est - param))

this_vim <- "AUC" #indx_vim_t_combos$vim[i]
this_indx <- "1,5"#indx_vim_t_combos$indx[i]
this_t <- 0.5#indx_vim_t_combos$t[i]

this_dat <- dat %>% filter(vim == this_vim & tau == this_t & indx == this_indx)
blah <- this_dat %>%
  group_by(misspec_type, n_train, indx, robust_V, robust_f) %>%
  summarize(nreps = n(),
            mean_err = mean(err),
            mse = mean(err^2),
            sd_err = sd(err),
            se_err = sqrt( mean((err - mean_err) ^ 2) / (nreps-1))) %>%
  mutate(misspec_type = factor(misspec_type,
                               levels = c("none", "censoring", "event_plusf0"),
                               labels = c("None", "Censoring", "Event")),
         robust_V = factor(robust_V, levels = c(FALSE, TRUE),
                           labels = c("Direct", "Indirect")),
         robust_f = factor(robust_f, levels = c(FALSE, TRUE),
                           labels = c("Conditional surv. function", "doubly-robust pseudo-outcome")))

# blah <- this_dat %>% filter(t == 0.5) %>%
#   mutate(misspec_type = factor(misspec_type,
#                                levels = c("none", "censoring", "event_plusf0"),
#                                labels = c("None", "Censoring", "Event")),
#          robust_V = factor(robust_V, levels = c(FALSE, TRUE),
#                            labels = c("Direct", "Indirect")),
#          robust_f = factor(robust_f, levels = c(FALSE, TRUE),
#                            labels = c("Conditional surv. function", "DR pseudo-outcome")))


# p1 <- blah %>% filter(misspec_type == "None") %>%
#   ggplot(aes(y = err, x = factor(n_train))) +
#   geom_histogram() +
#   facet_grid(robust_V ~robust_f) +
#   theme_bw() +
#   geom_hline(yintercept = 0)
# p2 <- blah %>% filter(misspec_type == "Censoring") %>%
#   ggplot(aes(y = err, x = factor(n_train))) +
#   geom_boxplot() +
#   facet_grid(robust_V ~robust_f) +
#   theme_bw() +
#   geom_hline(yintercept = 0)
# p3 <- blah %>% filter(misspec_type == "Event") %>%
#   ggplot(aes(y = err, x = factor(n_train))) +
#   geom_boxplot() +
#   facet_grid(robust_V ~robust_f) +
#   theme_bw() +
#   geom_hline(yintercept = 0)



p <- blah %>% ggplot(aes(x = factor(n_train), y = mean_err )) +
  geom_point(aes(shape = robust_V), size = 2.5, color = "grey40") +
  geom_line(aes(linetype = robust_f,
                group = interaction(robust_f, robust_V)),
            color = "grey40") +
  # geom_errorbar(aes(ymin = mean_err - sd_err,
  # ymax=mean_err + sd_err),
  # width=100) +
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
         axis.title = element_text(size = 14, family = "Times New Roman")) +
  guides(linetype=guide_legend(title="Oracle estimator:", nrow = 1, ncol = 2),
         shape = guide_legend(title = "Debiasing:", nrow = 1, ncol = 2))

ggsave(filename = "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/biometrika/robust_mse.pdf",
       plot = p, device = "pdf",
       width = 12, height = 4, dpi = 300, units = "in")

# summarize results
# dat <- dat %>% mutate(err = (estimate - param))

# dat <- dat %>% filter(tau == 0.5)
#
# dat_robust <- dat %>% filter(robust) %>%
#   select(n_train, err, misspec_type) %>%
#   mutate(id = 1:n(),
#          robust_err = err) %>%
#   select(-err)
# dat_nonrobust <- dat %>% filter(!robust) %>%
#   select(n_train, err, misspec_type) %>%
#   mutate(id = 1:n(),
#          nonrobust_err = err) %>%
#   select(-err)
# dat <- inner_join(dat_robust, dat_nonrobust, by = c("id", "misspec_type", "n_train")) %>%
#   mutate(misspec_type = factor(misspec_type, levels = c("none",
#                                                         "censoring",
#                                                         "event"),
#                                labels = c("None",
#                                           "Censoring",
#                                           "Event")))
#
# fills <- c("Indirect one-step" = "grey33", "Direct one-step" = "grey75")
# p <- dat %>%
#   ggplot() +
#   geom_histogram(aes(x = robust_err, fill = "Indirect one-step", y = after_stat(count / sum(count))),
#                  bins = 50) +
#   geom_histogram(aes(x = nonrobust_err, fill = "Direct one-step", y = -after_stat(count / sum(count))),
#                  bins = 50) +
#   scale_fill_manual(values = fills, breaks=c("Indirect one-step", "Direct one-step")) +
#   facet_grid(factor(n_train)~misspec_type) +
#   geom_vline(xintercept = 0) +
#   theme_bw() +
#   ggtitle("Misspecification type") +
#   xlab("Error") +
#   ylab("Frequency") +
#   scale_y_continuous(breaks = c(-0.015, -0.01, -0.005, 0, 0.005, 0.01, 0.015),
#                      labels = c("0.015", "0.01", "0.005", "0", "0.005", "0.01", "0.015"),
#                      sec.axis = sec_axis(~ . , name = "Sample size",
#                                          labels = NULL, breaks = NULL)) +
#   theme( strip.background = element_blank(),
#          strip.placement = "outside",
#          plot.title = element_text(hjust = 0.5, size = 12),
#          legend.position="bottom",
#          panel.spacing.x = unit(0.5, "cm"),
#          axis.text = element_text(size = 9)) +
# #   guides(fill=guide_legend(title="Estimator type"))
#
# ggsave(filename = "C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/biometrika/robust.pdf",
#        plot = p,
#        device = "pdf",
#        width = 12,
#        height = 5,
#        dpi = 300,
#        units = "in")
