source("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/sims/figure_utils.R")

dat <- readRDS("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/survex_comparison_brier.rds")
dat <- dat %>% filter(landmark_time == 0.5) %>%
  mutate(method = factor(method, levels = c("permutation", "exclusion"),
                         labels = c("Permutation", "Exclusion")),
         scenario = factor(scenario, levels = c("permA", "permB", "permC", "permD"),
                           labels = c("0", "0.3", "0.6", "0.9")),
         n_train = factor(n_train))

summ <- dat %>% group_by(n_train, method,scenario) %>%
  summarize(nreps = n(),
            rank_right = mean(correct),
            rank_right_mc_se = sqrt(rank_right * (1 - rank_right) / nreps),
            runtime = mean(runtime),
            bias = mean(est))

p <- summ %>% ggplot(aes(x = n_train, y = rank_right, group = interaction(scenario, method))) +
  geom_line(aes(linetype = scenario)) +
  geom_errorbar(aes(ymin=rank_right - 1.96*rank_right_mc_se,
                    ymax=rank_right + 1.96*rank_right_mc_se),
                width=.1) +
  geom_point(size = 1) +
  facet_wrap(~ method) +
  theme_bw() +
  scale_color_manual(values = c("black", "blue")) +
  ylim(c(0, 1)) +
  ylab("Proportion of replicates ranked correctly") +
  xlab("Sample size") +
  ggtitle("Variable importance type") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        plot.title = element_text(hjust = 0.5, size = 16,family = "Times New Roman"),
        legend.title = element_text(size = 16,family = "Times New Roman"),
        legend.text = element_text(size = 14,family = "Times New Roman"),
        legend.position="bottom",
        axis.text = element_text(size = 14,family = "Times New Roman"),
        strip.text = element_text(size = 14, family = "Times New Roman"),
        axis.title = element_text(size = 14, family = "Times New Roman"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(linetype=guide_legend(title="Correlation:", nrow = 1, ncol = 4))

ggsave(filename = "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/biometrika/survex-comparison-071724.pdf",
       plot = p, device = "pdf",
       width = 10, height = 5, dpi = 300, units = "in")
