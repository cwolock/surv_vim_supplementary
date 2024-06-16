#!/usr/local/bin/Rscript

# packages
library("cowplot")
library("grid")
library("gridExtra")
library("ggpubr")
library("tidyverse")
library("extrafont")

# font_import(prompt = FALSE)
loadfonts(device = "all")
myFont <- "Times New Roman"

# plotting parameters
point_size <- 1
title_text_size <- 12
legend_text_size <- 12
axis_text_size <- 12
big_fig_width <- 12
big_fig_height <- 8
small_fig_width <- 12
small_fig_height <- 6
# big_fig_width <- 6
# big_fig_height <- 4
# small_fig_width <- 6
# small_fig_height <- 4
nuisance_cols <- c("red", "blue", "green")
xfit_linetypes <- c("solid", "longdash")
strip_text_size <- 10

compile_truth <- function(true_param_file, true_avar_file){
  # get true VIM values
  truth <- readRDS(true_param_file)
  truth$tau <- truth$t

  # this is a MC estimate of the true asymptotic variance of the estimator using the IF
  # with known nuisances plugged in
  var_truth <- readRDS(true_avar_file)
  var_truth$tau <- var_truth$t
  var_truth <- var_truth %>% group_by(vim, tau, scenario) %>%
    summarize(var_1 = mean(vim_1),
              var_2 = mean(vim_2),
              var_1_split = mean(vim_1_split),
              var_4_split = mean(vim_4_split),
              var_14_split = mean(vim_14_split),
              .groups = "drop") %>%
    pivot_longer(cols = c(var_1, var_2, var_1_split, var_4_split, var_14_split),
                 names_to = "indx",
                 values_to = "true_avar") %>%
    filter(!(!grepl("split", indx) & scenario == "4")) %>%
    mutate(
      scenario = case_when(
        grepl("split", indx) & scenario == "1" ~ "2",
        !grepl("split", indx) & scenario == "1" ~ "1",
        scenario == "4" ~ "4"
      ), indx = case_when(
        indx == "var_4" | indx == "var_4_split" ~ "4",
        indx == "var_14" | indx == "var_14_split" ~ "1,4",
        indx == "var_1" | indx == "var_1_split" ~ "1",
        indx == "var_2" | indx == "var_2_split" ~ "2"
      ))

  return(list(truth = truth, var_truth = var_truth))
}

summarize_results <- function(dat, scenario, truth, var_truth){

  dat <- dat %>% mutate(scenario = scenario)

  if (scenario != "3"){
    dat <- dat %>% mutate(cens_rate = "50%")
  } else {
    dat <- dat %>% mutate(scenario = "3", n_train = 1000)
  }

  dat <- dat %>% #pivot_longer(cols = c("one_step"),
                #              names_to = "estimator",
                #              values_to = "estimate") %>%
    mutate(tau = ifelse(vim == "cindex", tau, landmark_time),
           nuisance = factor(dat$nuisance,
                             levels = c("rfsrc",
                                        "stackG",
                                        "survSL"),
                             labels = c("random surv. forest",
                                        "global surv. stacking",
                                        "surv. Super Learner")),
           n_eff = ifelse(scenario == "1",
                          dat$n_train, dat$n_train/2)) # for sample splitting)

  dat <- dat %>% mutate(correlation = ifelse(scenario == "4", TRUE, FALSE))
  dat <- left_join(dat, truth, by = c("tau", "vim", "correlation"))

  dat <- dat %>% mutate(param = case_when(
    indx == 1 ~ V_full - V_01,
    indx == 2 ~ V_full - V_02,
    indx == 5 ~ V_full - V_05,
    indx == "1,5" ~ V_full - V_015
  ))

  # summarize results
  dat <- dat %>% mutate(err = (est - param),
                        cov = ifelse(cil <= param & ciu >= param, 1, 0),
                        reject = ifelse(p < 0.05, 1, 0),
                        width = (est - cil)*2)

  messed_up_landmark <- dat %>% filter(abs(err) > 1 | is.na(err))

  summ <- dat %>% group_by(tau, n_train, n_eff, nuisance, vim, indx, crossfit, scenario, cens_rate) %>%
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
              var_mc_se = sqrt(2/(nreps-1)),#variance/sqrt(2*(nreps-1)),# # this may be wrong
              .groups = "drop")
  if (scenario == "3"){
    summ <- summ %>% mutate(scaled_var = 100*variance,
                            var_mc_se = 1000*variance/sqrt(2*(nreps -1)))
  } else{
    summ <- summ %>% mutate(bias = sqrt(n_eff)*bias,
                            variance = n_eff * variance,
                            bias_mc_se = sqrt(n_eff)*bias_mc_se)
    summ <- left_join(summ, var_truth, by = c("tau", "vim", "indx", "scenario"))
    summ <- summ %>% mutate(scaled_var = variance / true_avar)
  }
  return(summ)
}

make_sim_plot <- function(summ, scenario, big = TRUE, wd, fname){
  plot_tib <- summ %>%
    mutate(n = factor(n_train),
           Method = factor(crossfit,
                           levels = c(FALSE, TRUE),
                           labels = c("Not cross-fit", "Cross-fit"))) %>%
    unite(vim, c("vim", "tau"), sep="_") %>%
    mutate(vim = factor(vim,
                        levels = c("AUC_0.5", "AUC_0.9",
                                   "brier_0.5", "brier_0.9",
                                   "rsquared_0.5", "rsquared_0.9",
                                   "cindex_0.9"),
                        labels = c("AUC at 0.5", "AUC at 0.9",
                                   "Brier score at 0.5", "Brier score at 0.9",
                                   "R-squared at 0.5", "R-squared at 0.9",
                                   "C-index")))

  if (!big){

    plot_tib <- plot_tib %>% mutate(indx = factor(indx,
                                                  levels = c("1", "4"),
                                                  labels = c("(a)~X[1]~(non-null)", "(b)~X[4]~(null)")))
    Switch <- FALSE
    xvar <- "n"
    xlab <- "Sample size"
    scales <- "fixed"

    this_bias_lim <- plot_tib %>%
      mutate(scaled_bias_low = bias - 2*bias_mc_se,
             scaled_bias_hi = bias + 2*bias_mc_se) %>%
      summarize(bias_l = min(scaled_bias_low), bias_h = max(scaled_bias_hi))
    this_bias_lim_lower <- this_bias_lim$bias_l
    this_bias_lim_upper <- this_bias_lim$bias_h

    this_width_lim <- plot_tib %>%
      mutate(width_low = ci_width - 2*width_mc_se,
             width_hi = ci_width + 2*width_mc_se) %>%
      summarize(width_l = min(width_low), width_h = max(width_hi))
    this_width_lim_lower <- this_width_lim$width_l
    this_width_lim_upper <- this_width_lim$width_h

    this_cov_lim <- plot_tib %>%
      mutate(cov_low = coverage - 2 * cov_mc_se,
             cov_hi = coverage + 2*cov_mc_se) %>%
      summarize(cov_l = min(cov_low), cov_h = max(cov_hi))
    this_cov_lim_lower <- this_cov_lim$cov_l
    this_cov_lim_upper <- this_cov_lim$cov_h

    this_var_lim <- plot_tib %>%
      mutate(var_low = scaled_var - 2*var_mc_se,
             var_hi = scaled_var + 2*var_mc_se) %>%
      summarize(var_l = min(var_low), var_h = max(var_hi))
    this_var_lim_lower <- this_var_lim$var_l
    this_var_lim_upper <- this_var_lim$var_h

    this_power_lim <- plot_tib %>%
      mutate(power_low = power - 2*power_mc_se,
             power_hi = power + 2*power_mc_se) %>%
      summarize(power_l = min(power_low), power_h = max(power_hi))
    this_power_lim_lower <- this_power_lim$power_l
    this_power_lim_upper <- this_power_lim$power_h

    ylab_bias <- ifelse(Switch, "Empirical bias", expression(paste(sqrt(n), " x empirical bias")))
    bias_plot_1 <- plot_tib %>% filter(indx == "(a)~X[1]~(non-null)") %>%
      ggplot(aes(x = eval(str2lang(xvar)), y = bias, color = nuisance)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      geom_point(size = point_size) +
      geom_errorbar(aes(ymin=bias-1.96*bias_mc_se, ymax=bias + 1.96*bias_mc_se), width=.1) +
      geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
      scale_color_manual(values = nuisance_cols) +
      scale_linetype_manual(values = xfit_linetypes) +
      ylim(c(this_bias_lim_lower, this_bias_lim_upper)) +
      ylab(ylab_bias) +
      xlab(xlab) +
      labs(linetype = "Method:", color = "Nuisance:") +
      # ggtitle("Bias") +
      theme_bw() +
      theme(#axis.ticks.length.x = unit(0, "cm"),
            # axis.text.x = element_blank(),
            # axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey85"))
    bias_plot_4 <- plot_tib %>% filter(indx == "(b)~X[4]~(null)") %>%
      ggplot(aes(x = eval(str2lang(xvar)), y = bias, color = nuisance)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      geom_point(size = point_size) +
      geom_errorbar(aes(ymin=bias-1.96*bias_mc_se, ymax=bias + 1.96*bias_mc_se), width=.1) +
      geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
      scale_color_manual(values = nuisance_cols) +
      scale_linetype_manual(values = xfit_linetypes) +
      ylim(c(this_bias_lim_lower, this_bias_lim_upper)) +
      ylab(ylab_bias) +
      xlab(xlab) +
      # ggtitle("Bias") +
      labs(linetype = "Method:", color = "Nuisance:") +
      theme_bw() +
      theme(#axis.ticks.length.x = unit(0, "cm"),
            plot.title = element_text(hjust = 0.5),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey85"))

    cover_plot_1 <- plot_tib %>% filter(indx == "(a)~X[1]~(non-null)") %>%
      ggplot(aes(x = eval(str2lang(xvar)), y = coverage, color = nuisance)) +
      geom_hline(yintercept = 0.95, linetype = "solid", color = "black") +
      geom_point(size = point_size) +
      geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
      geom_errorbar(aes(ymin=ifelse(coverage-1.96*cov_mc_se <0, 0, coverage-1.96*cov_mc_se),
                        ymax=ifelse(coverage+1.96*cov_mc_se >1, 1, coverage+1.96*cov_mc_se)),
                    width=.1) +
      scale_color_manual(values = nuisance_cols) +
      scale_linetype_manual(values = xfit_linetypes) +
      ylim(c(0, 1)) +
      ylab("Empirical coverage") +
      labs(linetype = "Method:", color = "Nuisance:") +
      xlab(xlab) +
      # ggtitle("Coverage") +
      theme_bw() +
      theme(#axis.ticks.length.x = unit(0, "cm"),
            # axis.text.x = element_blank(),
            # axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey85"))
    cover_plot_4 <- plot_tib %>% filter(indx == "(b)~X[4]~(null)") %>%
      ggplot(aes(x = eval(str2lang(xvar)), y = coverage, color = nuisance)) +
      geom_hline(yintercept = 0.95, linetype = "solid", color = "black") +
      geom_point(size = point_size) +
      geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
      geom_errorbar(aes(ymin=ifelse(coverage-1.96*cov_mc_se <0, 0, coverage-1.96*cov_mc_se),
                        ymax=ifelse(coverage+1.96*cov_mc_se >1, 1, coverage+1.96*cov_mc_se)),
                    width=.1) +
      scale_color_manual(values = nuisance_cols) +
      scale_linetype_manual(values = xfit_linetypes) +
      ylim(c(0, 1)) +
      ylab("Empirical coverage") +
      xlab(xlab) +
      # ggtitle("Coverage") +
      labs(linetype = "Method:", color = "Nuisance:") +
      theme_bw() +
      theme(#axis.ticks.length.x = unit(0, "cm"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_line(color = "grey85"))

    power_plot_4 <- plot_tib %>% filter(indx == "(b)~X[4]~(null)") %>%
      ggplot(aes(x = eval(str2lang(xvar)), y = power, color = nuisance)) +
      geom_point(size = point_size) +
      geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
      geom_errorbar(aes(ymin=ifelse(power-1.96*power_mc_se <0, 0, power-1.96*power_mc_se),
                        ymax=ifelse(power+1.96*power_mc_se >1, 1, power+1.96*power_mc_se)),
                    width=.1) +
      geom_hline(yintercept = 0.05, linetype = "solid", color = "black") +
      scale_color_manual(values = nuisance_cols) +
      scale_linetype_manual(values = xfit_linetypes) +
      ylim(c(this_power_lim_lower, this_power_lim_upper)) +
      facet_wrap(~ indx, labeller = label_parsed, strip.position = "right") +
      ylab("Empirical type I error") +
      xlab(xlab) +
      # ggtitle("Type I error") +
      labs(linetype = "Method:", color = "Nuisance:") +
      theme_bw() +
      theme(#axis.ticks.length.x = unit(0, "cm"),
            strip.background = element_blank(),
            strip.text = element_text(size = 16, family = "Times New Roman", face = "bold"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_line(color = "grey85"))

    ylab_var <- ifelse(Switch, "1000 x empirical var.", expression(paste(n, " x empirical var./(true asym. var.)")))
    var_plot_1 <- plot_tib  %>% filter(indx == "(a)~X[1]~(non-null)") %>%
      ggplot(aes(x = eval(str2lang(xvar)), y = scaled_var, color = nuisance)) +
      geom_point(size = point_size) +
      geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
      geom_errorbar(aes(ymin=scaled_var - 1.96*var_mc_se,
                        ymax=scaled_var + 1.96*var_mc_se),
                    width=.1) +
      {if(!Switch) geom_hline(yintercept = 1, linetype = "solid", color = "black")} +
      scale_color_manual(values = nuisance_cols)+
      scale_linetype_manual(values = xfit_linetypes) +
      ylim(c(this_var_lim_lower, ifelse(Switch, NA, this_var_lim_upper))) +
      ylab(ylab_var) +
      labs(linetype = "Method:", color = "Nuisance:") +
      # ggtitle("Variance") +
      theme_bw() +
      xlab(xlab) +
      theme(#axis.ticks.length.x = unit(0, "cm"),
            # axis.text.x = element_blank(),
            # axis.title.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey85"))
    var_plot_4 <- plot_tib  %>% filter(indx == "(b)~X[4]~(null)") %>%
      ggplot(aes(x = eval(str2lang(xvar)), y = scaled_var, color = nuisance)) +
      geom_point(size = point_size) +
      geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
      geom_errorbar(aes(ymin=scaled_var - 1.96*var_mc_se,
                        ymax=scaled_var + 1.96*var_mc_se),
                    width=.1) +
      {if(!Switch) geom_hline(yintercept = 1, linetype = "solid", color = "black")} +
      scale_color_manual(values = nuisance_cols)+
      scale_linetype_manual(values = xfit_linetypes) +
      ylim(c(this_var_lim_lower, ifelse(Switch, NA, this_var_lim_upper))) +
      ylab(ylab_var) +
      xlab(xlab) +
      # ggtitle("Variance") +
      labs(linetype = "Method:", color = "Nuisance:") +
      theme_bw() +
      theme(#axis.ticks.length.x = unit(0, "cm"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_line(color = "grey85"))

    width_plot_1 <-plot_tib %>% filter(indx == "(a)~X[1]~(non-null)") %>%
      ggplot(aes(x = eval(str2lang(xvar)), y = ci_width, color = nuisance)) +
      geom_point(size = point_size) +
      geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
      geom_errorbar(aes(ymin = ci_width - 1.96*width_mc_se,
                        ymax = ci_width + 1.96*width_mc_se),
                    width=.1) +
      scale_color_manual(values = nuisance_cols) +
      scale_linetype_manual(values = xfit_linetypes) +
      facet_wrap(~ indx, labeller = label_parsed,strip.position = "right") +
      ylim(c(this_width_lim_lower, this_width_lim_upper)) +
      ylab("Confidence interval width") +
      labs(linetype = "Method:", color = "Nuisance:") +
      # ggtitle("Width") +
      theme_bw() +
      xlab(xlab) +
      theme(#axis.ticks.length.x = unit(0, "cm"),
            # axis.text.x = element_blank(),
            # axis.title.x = element_blank(),
            strip.background = element_blank(),
            strip.text = element_text(size = 16, family = "Times New Roman", face = "bold"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_line(color = "grey85"))

    eight_panel_plot <- plot_grid(
      bias_plot_1 + theme(legend.position = "none",
                          title = element_text(size = title_text_size, family = "Times New Roman"),
                          axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                          axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                          plot.margin = unit(c(0.1, 0.3, 0.2, 0), "cm"),
                          strip.text = element_blank()),
      var_plot_1 + theme(legend.position = "none",
                         title = element_text(size = title_text_size, family = "Times New Roman"),
                         axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                         axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                         plot.margin = unit(c(0.1, 0.3, 0.2, 0), "cm"),
                         strip.text = element_blank()),
      cover_plot_1 + theme(legend.position = "none",
                           title = element_text(size = title_text_size, family = "Times New Roman"),
                           # axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                           axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                           axis.title.y = element_text(size = axis_text_size, margin = margin(t = 0, r = 5, b = 0, l = 0)),
                           plot.margin = unit(c(0.1, 0.3, 0.2, 0), "cm"),
                           strip.text = element_blank()),
      width_plot_1 + theme(legend.position = "none",
                           title = element_text(size = title_text_size, family = "Times New Roman"),
                           axis.title.x= element_text(size = axis_text_size, family = "Times New Roman"),
                           axis.title.y = element_text(size = axis_text_size, family = "Times New Roman",
                                                       margin = margin(t = 0, r = 5, b = 0, l = 0)),
                           axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                           plot.margin = unit(c(0.1, 0.1, 0.2, 0), "cm")),
      bias_plot_4 + theme(legend.position = "none",
                          title = element_text(size = title_text_size, family = "Times New Roman"),
                          axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                          axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                          plot.margin = unit(c(0.3, 0.3, 0, 0), "cm"),
                          strip.text = element_blank()),
      var_plot_4 + theme(legend.position = "none",
                         title = element_text(size = title_text_size, family = "Times New Roman"),
                         axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                         axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                         plot.margin = unit(c(0.3, 0.3, 0, 0), "cm"),
                         strip.text = element_blank()),
      cover_plot_4 + theme(legend.position = "none",
                           title = element_text(size = title_text_size, family = "Times New Roman"),
                           axis.title.x = element_text(size = axis_text_size, family = "Times New Roman"),
                           axis.title.y = element_text(size = axis_text_size, family = "Times New Roman",
                                                       margin = margin(t = 0, r = 5, b = 0, l = 0)),
                           axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                           plot.margin = unit(c(0.3, 0.3, 0, 0), "cm"),
                           strip.text = element_blank()),
      power_plot_4 + theme(legend.position = "none",
                           title = element_text(size = title_text_size, family = "Times New Roman"),
                           axis.title.x = element_text(size = axis_text_size, family = "Times New Roman"),
                           axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                           axis.title.y = element_text(size = axis_text_size, family = "Times New Roman",
                                                       margin = margin(t = 0, r = 5, b = 0, l = 0)),
                           plot.margin = unit(c(0.3, 0.1, 0, 0), "cm")),
      labels = NULL, nrow = 2, ncol = 4, rel_heights = c(0.5, 0.5), rel_widths = c(0.26, 0.26, 0.28, 0.31)
    )

    legend_j <- get_legend(
      bias_plot_1 +
        guides(linetype = guide_legend(nrow = 1, ncol = 2),
               color = guide_legend(nrow = 1, ncol = 3)) +
        theme(legend.direction = "horizontal",
              legend.position = "bottom",
              legend.title = element_text(size = legend_text_size, family = "Times New Roman"),
              legend.text = element_text(size = legend_text_size, family = "Times New Roman"))
    )
    full_plot_j <- plot_grid(eight_panel_plot, legend_j, ncol = 1, nrow = 2,
                             rel_heights = c(1, .1))
    ggsave(filename = paste0(wd, fname, ".pdf"),
           plot = full_plot_j, device = "pdf",
           width = ifelse(big, big_fig_width, small_fig_width),
           height = ifelse(big, big_fig_height, small_fig_height),
           dpi = 300, units = "in")
  } else {
    indxs <- unique(plot_tib$indx)

    Switch <- scenario == "3"
    xvar <- ifelse(Switch, "cens_rate", "n")
    xlab <- ifelse(Switch, "Censoring rate", "Sample size")
    scales <- ifelse(Switch, "free_y", "fixed")

    for (i in 1:length(indxs)){
      this_indx <- indxs[i]

      # if (sum(grepl("Brier", plot_tib$vim)) > 0){
      #   this_bias_lim_brier <- plot_tib %>% filter(indx == this_indx & grepl("Brier", vim)) %>%
      #     mutate(scaled_bias_low = bias - 2*bias_mc_se,
      #            scaled_bias_hi = bias + 2*bias_mc_se) %>%
      #     summarize(bias_l = min(scaled_bias_low), bias_h = max(scaled_bias_hi))
      #   this_bias_lim_lower_brier <- this_bias_lim$bias_l
      #   this_bias_lim_upper_brier <- this_bias_lim$bias_h
      # }
      this_bias_lim <- plot_tib %>% filter(indx == this_indx) %>%
        mutate(scaled_bias_low = bias - 2*bias_mc_se,
               scaled_bias_hi = bias + 2*bias_mc_se) %>%
        summarize(bias_l = min(scaled_bias_low), bias_h = max(scaled_bias_hi))
      this_bias_lim_lower <- this_bias_lim$bias_l
      this_bias_lim_upper <- this_bias_lim$bias_h

      this_width_lim <- plot_tib %>% filter(indx == this_indx) %>%
        mutate(width_low = ci_width - 2*width_mc_se,
               width_hi = ci_width + 2*width_mc_se) %>%
        summarize(width_l = min(width_low), width_h = max(width_hi))
      this_width_lim_lower <- this_width_lim$width_l
      this_width_lim_upper <- this_width_lim$width_h

      this_cov_lim <- plot_tib %>% filter(indx == this_indx) %>%
        mutate(cov_low = coverage - 2 * cov_mc_se,
               cov_hi = coverage + 2*cov_mc_se) %>%
        summarize(cov_l = min(cov_low), cov_h = max(cov_hi))
      this_cov_lim_lower <- this_cov_lim$cov_l
      this_cov_lim_upper <- this_cov_lim$cov_h

      this_var_lim <- plot_tib %>% filter(indx == this_indx) %>%
        mutate(var_low = scaled_var - 2*var_mc_se,
               var_hi = scaled_var + 2*var_mc_se) %>%
        summarize(var_l = min(var_low), var_h = max(var_hi))
      this_var_lim_lower <- this_var_lim$var_l
      this_var_lim_upper <- this_var_lim$var_h

      this_power_lim <- plot_tib %>% filter(indx == this_indx) %>%
        mutate(power_low = power - 2*power_mc_se,
               power_hi = power + 2*power_mc_se) %>%
        summarize(power_l = min(power_low), power_h = max(power_hi))
      this_power_lim_lower <- this_power_lim$power_l
      this_power_lim_upper <- this_power_lim$power_h

      num_rows <- plot_tib %>% filter(indx == this_indx) %>%
        pull(vim)
      num_rows <- length(unique(num_rows))

      ylab_bias <- ifelse(Switch, "Empirical bias", expression(paste(sqrt(n), " x empirical bias")))
      bias_plot_j <- plot_tib %>% filter(indx == this_indx) %>%
        ggplot(aes(x = eval(str2lang(xvar)), y = bias, color = nuisance)) +
        geom_hline(yintercept = 0, linetype = "solid", color = "black") +
        geom_point(size = point_size) +
        geom_errorbar(aes(ymin=bias-1.96*bias_mc_se, ymax=bias + 1.96*bias_mc_se), width=.1) +
        geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
        scale_color_manual(values = nuisance_cols) +
        scale_linetype_manual(values = xfit_linetypes) +
        ylim(c(this_bias_lim_lower, this_bias_lim_upper)) +
        ylab(ylab_bias) +
        xlab(xlab) +
        labs(linetype = "Method:", color = "Nuisance:") +
        # ggtitle("A. BIAS") +
        ggtitle("(a)") +
        facet_wrap(~vim, nrow = num_rows, ncol = 1, strip.position = "right") +
        theme_bw() +
        theme(#axis.ticks.length.x = unit(0, "cm"),
              panel.spacing.x = unit(0, "cm"),
              panel.spacing.y = unit(0.2, "cm"),
              strip.background = element_blank(),
              strip.placement = "outside",
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "grey85"))

      cover_plot_j <- plot_tib %>% filter(indx == this_indx) %>%
        ggplot(aes(x = eval(str2lang(xvar)), y = coverage, color = nuisance)) +
        geom_hline(yintercept = 0.95, linetype = "solid", color = "black") +
        geom_point(size = point_size) +
        geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
        geom_errorbar(aes(ymin=ifelse(coverage-1.96*cov_mc_se <0, 0, coverage-1.96*cov_mc_se),
                          ymax=ifelse(coverage+1.96*cov_mc_se >1, 1, coverage+1.96*cov_mc_se)),
                      width=.1) +
        scale_color_manual(values = nuisance_cols) +
        scale_linetype_manual(values = xfit_linetypes) +
        ylim(c(0, 1)) +
        ylab("Empirical coverage") +
        xlab(xlab) +
        labs(linetype = "Method:", color = "Nuisance:") +
        # ggtitle("C. COVERAGE") +
        ggtitle("(c)") +
        facet_wrap(~vim, nrow = num_rows, ncol = 1, strip.position = "right") +
        theme_bw() +
        theme(#axis.ticks.length.x = unit(0, "cm"),
              panel.spacing.x = unit(0, "cm"),
              panel.spacing.y = unit(0.2, "cm"),
              strip.background = element_blank(),
              strip.placement = "outside",
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "grey85"))

      power_plot_j <- plot_tib %>% filter(indx == this_indx) %>%
        ggplot(aes(x = eval(str2lang(xvar)), y = power, color = nuisance)) +
        geom_point(size = point_size) +
        geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
        geom_errorbar(aes(ymin=ifelse(power-1.96*power_mc_se <0, 0, power-1.96*power_mc_se),
                          ymax=ifelse(power+1.96*power_mc_se >1, 1, power+1.96*power_mc_se)),
                      width=.1) +
        geom_hline(yintercept = 0.05, linetype = "solid", color = "black") +
        scale_color_manual(values = nuisance_cols) +
        scale_linetype_manual(values = xfit_linetypes) +
        ylim(c(this_power_lim_lower, this_power_lim_upper)) +
        ylab("Empirical type I error") +
        xlab(xlab) +
        labs(linetype = "Method:", color = "Nuisance:") +
        # ggtitle("D. TYPE I ERROR") +
        ggtitle("(d)") +
        facet_wrap(~vim, nrow = num_rows, ncol = 1, strip.position = "right") +
        theme_bw() +
        theme(#axis.ticks.length.x = unit(0, "cm"),
              panel.spacing.x = unit(0, "cm"),
              panel.spacing.y = unit(0.2, "cm"),
              strip.background = element_blank(),
              strip.placement = "outside",
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "grey85"))

      ylab_var <- ifelse(Switch, "1000 x empirical var.", expression(paste(n, " x empirical var./(true asymp. var.)")))
      var_plot_j <- plot_tib  %>% filter(indx == this_indx) %>%
        ggplot(aes(x = eval(str2lang(xvar)), y = scaled_var, color = nuisance)) +
        geom_point(size = point_size) +
        geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
        geom_errorbar(aes(ymin=scaled_var - 1.96*var_mc_se,
                          ymax=scaled_var + 1.96*var_mc_se),
                      width=.1) +
        {if(!Switch) geom_hline(yintercept = 1, linetype = "solid", color = "black")} +
        # geom_hline(yintercept = 1, linetype = "solid", color = "black") +
        scale_color_manual(values = nuisance_cols)+
        scale_linetype_manual(values = xfit_linetypes) +
        # ylim(c(this_var_lim_lower, this_var_lim_upper)) +
        ylim(c(this_var_lim_lower, ifelse(Switch, NA, this_var_lim_upper))) +
        ylab(ylab_var) +
        xlab(xlab) +
        labs(linetype = "Method:", color = "Nuisance:") +
        # ggtitle("B. VARIANCE") +
        ggtitle("(b)") +
        facet_wrap(~vim, nrow = num_rows, ncol = 1, strip.position = "right", scales = scales) +
        theme_bw() +
        theme(#axis.ticks.length.x = unit(0, "cm"),
              panel.spacing.x = unit(0, "cm"),
              panel.spacing.y = unit(0.2, "cm"),
              strip.background = element_blank(),
              strip.placement = "outside",
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "grey85"))

      width_plot_j <-plot_tib %>% filter(indx == this_indx) %>%
        ggplot(aes(x = eval(str2lang(xvar)), y = ci_width, color = nuisance)) +
        geom_point(size = point_size) +
        geom_line(aes(group = interaction(nuisance, Method), linetype = Method)) +
        geom_errorbar(aes(ymin = ci_width - 1.96*width_mc_se,
                          ymax = ci_width + 1.96*width_mc_se),
                      width=.1) +
        scale_color_manual(values = nuisance_cols) +
        scale_linetype_manual(values = xfit_linetypes) +
        ylim(c(this_width_lim_lower, this_width_lim_upper)) +
        ylab("Confidence interval width") +
        xlab(xlab) +
        labs(linetype = "Method:", color = "Nuisance:") +
        # ggtitle("D. WIDTH") +
        ggtitle("(d)") +
        facet_wrap(~vim, nrow = num_rows, ncol = 1, strip.position = "right") +
        theme_bw() +
        theme(#axis.ticks.length.x = unit(0, "cm"),
              panel.spacing.x = unit(0, "cm"),
              panel.spacing.y = unit(0.2, "cm"),
              strip.background = element_blank(),
              strip.placement = "outside",
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "grey85"))

      if (this_indx == "4"){
        four_panel_plot_j <- plot_grid(
          bias_plot_j + theme(legend.position = "none",
                              title = element_text(size = title_text_size, family = "Times New Roman"),
                              axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                              axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                              plot.margin = unit(c(0, 0.3, 0.1, 0), "cm"),
                              strip.text = element_blank()),
          var_plot_j + theme(legend.position = "none",
                             title = element_text(size = title_text_size, family = "Times New Roman"),
                             axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                             axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                             plot.margin = unit(c(0, 0.3, 0.1, 0), "cm"),
                             strip.text = element_blank()),
          cover_plot_j + theme(legend.position = "none",
                               title = element_text(size = title_text_size, family = "Times New Roman"),
                               axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                               axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                               plot.margin = unit(c(0, 0.3, 0.1, 0), "cm"),
                               strip.text = element_blank()),
          power_plot_j + theme(legend.position = "none",
                               title = element_text(size = title_text_size, family = "Times New Roman"),
                               axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                               axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                               plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                               strip.text = element_text(size = strip_text_size, family = "Times New Roman")),
          labels = NULL, nrow = 1, ncol = 4
        )
      } else{
        four_panel_plot_j <- plot_grid(
          bias_plot_j + theme(legend.position = "none",
                              title = element_text(size = title_text_size, family = "Times New Roman"),
                              axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                              axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                              plot.margin = unit(c(0, 0.3, 0.1, 0), "cm"),
                              strip.text = element_blank()),
          var_plot_j + theme(legend.position = "none",
                             title = element_text(size = title_text_size, family = "Times New Roman"),
                             axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                             axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                             plot.margin = unit(c(0, 0.3, 0.1, 0), "cm"),
                             strip.text = element_blank()),
          cover_plot_j + theme(legend.position = "none",
                               title = element_text(size = title_text_size, family = "Times New Roman"),
                               axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                               axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                               plot.margin = unit(c(0, 0.3, 0.1, 0), "cm"),
                               strip.text = element_blank()),
          width_plot_j + theme(legend.position = "none",
                               title = element_text(size = title_text_size, family = "Times New Roman"),
                               axis.title = element_text(size = axis_text_size, family = "Times New Roman"),
                               axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
                               plot.margin = unit(c(0, 0, 0.1, 0), "cm"),
                               strip.text = element_text(size = strip_text_size, family = "Times New Roman")),
          labels = NULL, nrow = 1, ncol = 4
        )
      }

      legend_j <- get_legend(
        bias_plot_j +
          guides(linetype = guide_legend(nrow = 1, ncol = 2),
                 color = guide_legend(nrow = 1, ncol = 3)) +
          theme(legend.direction = "horizontal",
                legend.position = "bottom",
                legend.title = element_text(size = legend_text_size, family = "Times New Roman"),
                legend.text = element_text(size = legend_text_size, family = "Times New Roman"))
      )
      full_plot_j <- plot_grid(four_panel_plot_j, legend_j, ncol = 1, nrow = 2,
                               rel_heights = c(1, .1))
      ggsave(filename = paste0(wd, fname, "-",
                               this_indx, ".pdf"),
             plot = full_plot_j, device = "pdf",
             width = ifelse(big, big_fig_width, small_fig_width),
             height = ifelse(big, big_fig_height, small_fig_height),
             dpi = 300, units = "in")
    }
  }


}
