source("C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/figure_utils.R")

# read in truth files
truth_file <- "C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim/scratch/sims/landmark/truth.rds"
var_truth_file <- "C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/variance_truth.rds"
truth_list <- compile_truth(true_param_file = truth_file,
                            true_avar_file = var_truth_file)

dat <- readRDS("C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario5.rds")

dat <- dat %>% pivot_longer(cols = c("one_step"),
                            names_to = "estimator",
                            values_to = "estimate")

dat <- dat %>% mutate(correlation = FALSE)
dat <- left_join(dat, truth_list$truth, by = c("tau", "vim", "correlation"))

dat <- dat %>% mutate(param = case_when(
  indx == 1 ~ V_full - V_01,
  indx == 2 ~ V_full - V_02,
  indx == 4 ~ 0,
  indx == "1,4" ~ V_full - V_014
))

# summarize results
dat <- dat %>% mutate(err = (estimate - param))

dat <- dat %>% filter(tau == 0.5)

dat_robust <- dat %>% filter(robust) %>%
  select(n_train, err, misspec_type) %>%
  mutate(id = 1:n(),
         robust_err = err) %>%
  select(-err)
dat_nonrobust <- dat %>% filter(!robust) %>%
  select(n_train, err, misspec_type) %>%
  mutate(id = 1:n(),
         nonrobust_err = err) %>%
  select(-err)
dat <- inner_join(dat_robust, dat_nonrobust, by = c("id", "misspec_type", "n_train")) %>%
  mutate(misspec_type = factor(misspec_type, levels = c("none",
                                                        "censoring",
                                                        "event"),
                               labels = c("None",
                                          "Censoring",
                                          "Event")))

fills <- c("Indirect one-step" = "grey33", "Direct one-step" = "grey75")
p <- dat %>%
  ggplot() +
  geom_histogram(aes(x = robust_err, fill = "Indirect one-step", y = after_stat(count / sum(count))),
                 bins = 50) +
  geom_histogram(aes(x = nonrobust_err, fill = "Direct one-step", y = -after_stat(count / sum(count))),
                 bins = 50) +
  scale_fill_manual(values = fills, breaks=c("Indirect one-step", "Direct one-step")) +
  facet_grid(factor(n_train)~misspec_type) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  ggtitle("Misspecification type") +
  xlab("Error") +
  ylab("Frequency") +
  scale_y_continuous(breaks = c(-0.015, -0.01, -0.005, 0, 0.005, 0.01, 0.015),
                     labels = c("0.015", "0.01", "0.005", "0", "0.005", "0.01", "0.015"),
                     sec.axis = sec_axis(~ . , name = "Sample size",
                                         labels = NULL, breaks = NULL)) +
  theme( strip.background = element_blank(),
         strip.placement = "outside",
         plot.title = element_text(hjust = 0.5, size = 12),
         legend.position="bottom",
         panel.spacing.x = unit(0.5, "cm"),
         axis.text = element_text(size = 9)) +
  guides(fill=guide_legend(title="Estimator type"))

ggsave(filename = "C:/Users/cwolo/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/scenario5.pdf",
       plot = p,
       device = "pdf",
       width = 8,
       height = 11,
       dpi = 300,
       units = "in")
