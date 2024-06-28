source("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/figure_utils_newA.R")

# read in truth files
truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/truth_interaction.rds"
var_truth_file <- "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/variance_truth.rds"
truth_list <- compile_truth(true_param_file = truth_file,
                            true_avar_file = var_truth_file)

dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/cindex_subsampling.rds")

dat <- dat %>% #pivot_longer(cols = c("one_step"),
  #              names_to = "estimator",
  #              values_to = "estimate") %>%
  mutate(tau = restriction_time)

dat <- dat %>% mutate(correlation = TRUE)
dat <- left_join(dat, truth_list$truth, by = c("tau", "vim", "correlation"))

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

summ <- dat %>% group_by(tau, n_train, nuisance, vim, indx, crossfit, subsample) %>%
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
summ <- summ %>% mutate(bias = sqrt(n_eff)*bias,
                        variance = n_eff * variance,
                        bias_mc_se = sqrt(n_eff)*bias_mc_se)
