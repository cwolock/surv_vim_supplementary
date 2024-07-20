source("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/sims/figure_utils.R")

# read in truth files
truth_file <- "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/truth_interactionB.rds"
var_truth_file <- "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/mc_truth_variance_interactionB.rds"
truth_list <- compile_truth(true_param_file = truth_file,
                            true_avar_file = var_truth_file)

dat <- readRDS("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/cindex_subsampling.rds")

dat <- dat %>% 
  mutate(tau = restriction_time)

dat <- dat %>% mutate(correlation = TRUE)
dat <- left_join(dat, truth_list$truth, by = c("tau", "vim", "correlation"))

dat <- dat %>% mutate(param = case_when(
  indx == 1 ~ V_full - V_01,
  indx == 2 ~ V_full - V_02,
  indx == 6 ~ V_full - V_06,
  indx == "1,6" ~ V_full - V_016
))

# summarize results
dat <- dat %>% mutate(err = (est - param),
                      cov = ifelse(cil <= param & ciu >= param, 1, 0),
                      reject = ifelse(p < 0.05, 1, 0),
                      width = (est - cil)*2,
                      n_eff = n_train/2,
                      n_eff = factor(n_eff),
                      n_train = factor(n_train),
                      subsample = factor(subsample,
                                         levels = c(0.25, 0.33, 0.5, 1),
                                         labels = c("1/4", "1/3", "1/2", "1")))

messed_up_landmark <- dat %>% filter(abs(err) > 1 | is.na(err))

summ <- dat %>% group_by(tau, n_train, n_eff, nuisance, vim, indx, crossfit, subsample) %>%
  summarize(runtime = mean(runtime),
            bias = mean(err),
            nreps = n(),
            variance = var(est),
            bias_mc_se = sqrt( mean((err - bias) ^ 2) / (nreps-1)),
            coverage = mean(cov),
            power = mean(reject),
            cov_mc_se = sqrt(coverage * (1 - coverage)/nreps),
            power_mc_se = sqrt(power * (1 - power)/nreps),
            ci_width = mean(width),
            width_mc_se = sqrt( mean((width - ci_width) ^ 2) / (nreps-1)), # note that this is not really a "MC standard error" just the standard deviation
            # since the CI width isn't really estimating a parameter.
            var_mc_se = sqrt(2/(nreps-1)),#variance/sqrt(2*(nreps-1)),#
            .groups = "drop")

subsample_linetypes <- c("solid", "longdash", "dashed", "dotted")
point_size <- 1
title_text_size <- 12
legend_text_size <- 12
axis_text_size <- 12
fig_width <- 12
fig_height <- 8
strip_text_size <- 10


bias_plot <- summ %>% 
  ggplot(aes(x = n_train, y = bias)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_point(size = point_size) +
  geom_line(aes(group = subsample, linetype = subsample)) +
  # geom_errorbar(aes(ymin=bias-1.96*bias_mc_se, ymax=bias + 1.96*bias_mc_se), width=.1) +
  scale_linetype_manual(values = subsample_linetypes) +
  ylim(c(-0.015, 0.015)) + 
  ylab("Empirical bias") +
  xlab("Sample size") +
  labs(linetype = "Method:", color = "Nuisance:") +
  ggtitle("(a)") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "cm"),
    panel.spacing.y = unit(0.2, "cm"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank())

var_plot <- summ %>% 
  ggplot(aes(x = n_train, y = variance)) +
  geom_point(size = point_size) +
  geom_line(aes(group = subsample, linetype = subsample)) +
  # geom_errorbar(aes(ymin=variance-1.96*var_mc_se, ymax=variance + 1.96*var_mc_se), width=.1) +
  scale_linetype_manual(values = subsample_linetypes) +
  ylab("Empirical variance") +
  xlab("Sample size") +
  labs(linetype = "Method:", color = "Nuisance:") +
  ggtitle("(b)") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "cm"),
    panel.spacing.y = unit(0.2, "cm"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank())

cover_plot <- summ %>% 
  ggplot(aes(x = n_train, y = coverage)) +
  geom_hline(yintercept = 0.95, linetype = "solid", color = "black") +
  geom_point(size = point_size) +
  geom_line(aes(group = subsample, linetype = subsample)) +
  scale_linetype_manual(values = subsample_linetypes) +
  ylim(c(0.5, 1)) + 
  ylab("Empirical coverage") +
  xlab("Sample size") +
  labs(linetype = "Method:", color = "Nuisance:") +
  ggtitle("(c)") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "cm"),
    panel.spacing.y = unit(0.2, "cm"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank())

width_plot <- summ %>% 
  ggplot(aes(x = n_train, y = ci_width)) +
  geom_point(size = point_size) +
  geom_line(aes(group = subsample, linetype = subsample)) +
  scale_linetype_manual(values = subsample_linetypes) +
  ylab("Confidence interval width") +
  xlab("Sample size") +
  labs(linetype = "Method:", color = "Nuisance:") +
  ggtitle("(d)") +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0, "cm"),
    panel.spacing.y = unit(0.2, "cm"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor.y = element_blank())

four_panel_plot <- plot_grid(
  bias_plot + theme(legend.position = "none",
                      title = element_text(size = title_text_size, family = "Times New Roman"),
                      axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                      axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                      plot.margin = unit(c(0, 0.3, 0.1, 0), "cm"),
                      strip.text = element_blank()),
  var_plot + theme(legend.position = "none",
                     title = element_text(size = title_text_size, family = "Times New Roman"),
                     axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                     axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                     plot.margin = unit(c(0, 0.3, 0.1, 0), "cm"),
                     strip.text = element_blank()),
  cover_plot + theme(legend.position = "none",
                       title = element_text(size = title_text_size, family = "Times New Roman"),
                       axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                       axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                       plot.margin = unit(c(0, 0.3, 0.1, 0), "cm"),
                       strip.text = element_blank()),
  width_plot + theme(legend.position = "none",
                       title = element_text(size = title_text_size, family = "Times New Roman"),
                       axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                       axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                       plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                       strip.text = element_text(size = strip_text_size, family = "Times New Roman")),
  labels = NULL, nrow = 1, ncol = 4
)

legend <- get_legend(
  bias_plot +
    guides(linetype = guide_legend(title = "Subsampling proportion:", nrow = 1, ncol = 4)) +
    theme(legend.direction = "horizontal",
          legend.position = "bottom",
          legend.title = element_text(size = legend_text_size, family = "Times New Roman"),
          legend.text = element_text(size = legend_text_size, family = "Times New Roman"))
)
full_plot <- plot_grid(four_panel_plot, legend, ncol = 1, nrow = 2,
                         rel_heights = c(1, .1))
ggsave(filename = "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/biometrika/cindex-subsampling-072024.pdf",
       plot = full_plot, device = "pdf",
       width = 12, height = 4, dpi = 300, units = "in")

