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
point_size <- 2.5
title_text_size <- 15
axis_text_size <- 13
axis_title_size <- 15
big_fig_width <- 12
big_fig_height <- 7
strip_text_size <- 15

########################################
### combine results from multiple splits
########################################

aggregate_p <- function(ps, method = "bonferroni"){
  K <- length(ps)
  if (method == "bonferroni"){
    p_agg <- min(ps) * K
  } else if (method == "hommel"){
    constant <- sum(apply(matrix(1:K), MARGIN = 1, FUN = function(x) 1/x))
    sorted_ps <- sort(ps)
    corrected_ps <- apply(matrix(1:K),
                          MARGIN = 1,
                          FUN = function(x) K/x * sorted_ps[x])
    p_agg <- constant * min(corrected_ps)
  } else if (method == "arithmetic"){
    p_agg <- mean(ps) * 2
  } else if (method == "geometric"){
    p_agg <- exp(mean(log(ps))) * exp(1)
  } else if (method == "harmonic"){
    p_agg <- 1/mean(1/ps) * exp(1) * log(K)
  } else if (method == "compound_bg"){
    bonf <- min(ps) * K
    geo <- exp(mean(log(ps))) * exp(1)
    p_agg <- 2 * min(c(bonf, geo))
  } else if (method == "compound_ba"){
    bonf <- min(ps) * K
    arith <- mean(ps) * 2
    p_agg <- 2 * min(c(bonf, arith))
  }
  return(min(c(p_agg, 1)))
}

multisplit <- function(f, f1, method, this_vim = NULL){
  # f <- f %>% filter(vim == this_vim)
  # f1 <- f1 %>% filter(vim == this_vim)
  # look at multisplitting
  ci_grid <- seq(0, 0.6, by = 0.001)
  seeds <- unique(f$seed)
  n1 <- nrow(f1)
  n2 <- length(ci_grid)
  n3 <- length(seeds)
  p_array_twosided <- array(NA, dim = c(n1, n2, n3))
  p_array_onesided <- array(NA, dim = c(n1, n2, n3))
  for (i in 1:n3){
    dat <- f %>% filter(seed == seeds[i])
    for (j in 1:n2){
      for (k in 1:n1){
        z <- (dat$est[k] - ci_grid[j]) / sqrt(dat$var_est[k]/dat$n_eff[k])
        p <- pnorm(abs(z), lower.tail = FALSE)*2
        p_array_twosided[k, j, i] <- p
        p <- pnorm(z, lower.tail = FALSE)
        p_array_onesided[k, j, i] <- p
      }
    }
  }
  
  p_array_twosided_corrected <- array(NA, dim = c(n1, n2))
  p_array_onesided_corrected <- array(NA, dim = c(n1, n2))
  
  for (j in 1:n2){
    for (k in 1:n1){
      # this is just a bonferroni
      # p_two <- min(p.adjust(p_array_twosided[k,j,],method="bonferroni"))
      # p_one <- min(p.adjust(p_array_onesided[k,j,],method="bonferroni"))
      print(j)
      print(k)
      p_two <- aggregate_p(p_array_twosided[k,j,],method=method)
      p_one <- aggregate_p(p_array_onesided[k,j,],method=method)
      p_array_twosided_corrected[k,j] <- p_two
      p_array_onesided_corrected[k,j] <- p_one
    }
  }
  
  cil <- rep(NA, n1)
  ciu <- rep(NA, n1)
  ci_1sided <- rep(NA, n1)
  pvals <- rep(NA, n1)
  
  for (k in 1:n1){
    print(k)
    curr_ps_twosided <- p_array_twosided_corrected[k,]
    curr_ps_onesided <- p_array_onesided_corrected[k,]
    pvals[k] <- curr_ps_onesided[1]
    if (curr_ps_twosided[1] < 0.05){
      lower_index <- min(which(curr_ps_twosided > 0.05))
      cil[k] <- ci_grid[lower_index]
      reject_indices <- which(curr_ps_twosided < 0.05)
      upper_index <-  min(reject_indices[reject_indices > lower_index])
      ciu[k]<- ci_grid[upper_index]
    } else{
      cil[k] <- 0
      ciu[k]<- ci_grid[min(which(curr_ps_twosided < 0.05))]
    }
    
    if (curr_ps_onesided[1] < 0.05){
      ci_1sided[k] <- ci_grid[min(which(curr_ps_onesided > 0.05))]
    } else{
      ci_1sided[k] <- 0
    }
    
  }
  
  f_avg <- f %>% group_by(landmark_time, indx_name, vim, SAB) %>% summarize(est = mean(est))
  
  combined_ss <- f1 %>% select(landmark_time, vim, indx, indx_name, SAB) %>%
    mutate(cil = cil, ciu = ciu, ci_1sided = ci_1sided, pval = pvals)
  
  combined_ss <- left_join(combined_ss, f_avg, by = c("landmark_time", "indx_name", "vim", "SAB"))
  
  combined_ss <- combined_ss %>%
    filter(indx_name != "sex_bmi" & indx_name != "all_but_geo" & indx_name != "BRS") %>%
    mutate(cil = ifelse(cil >= 0, cil, 0),
           indx_number = case_when(
             indx_name == "sex" ~ "sex",
             indx_name == "age" ~ "age",
             indx_name == "bmi" ~ "BMI",
             indx_name == "BRS_sexhealth" ~ "sex. health",
             indx_name == "BRS_sexbehave" ~ "sex. behavior",
             indx_name == "BRS_social" ~ "housing",
             indx_name == "geo" ~ "geography"
           ))
  
  return(combined_ss)
}

method = "compound_bg"
unscale <- FALSE

### marginal combined
setwd("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/data_analysis/combined")
# f <- readRDS("analysis_112223.rds")
# f <- readRDS("combined_analysis_061324.rds")
f <- readRDS("combined_analysis_interactions_oldlib_smallersigma_subsample1750_oldseed.rds")
names(f)[names(f) == "tau"] <- "landmark_time"
f <- f %>% filter(approach == "marginal") %>% mutate(SAB = "Combined cohort")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(n_eff = 2641)
if (unscale){
  f <- f %>% mutate(est = ciu - 1.96* sqrt(var_est/n_eff))
}
combined_ss_both_marginal <- multisplit(f, f1, method)
combined_ss_both_marginal <- combined_ss_both_marginal

### conditional combined
# f <- readRDS("analysis_112223.rds")
f <- readRDS("combined_analysis_interactions_oldlib_smallersigma_subsample1750_oldseed.rds")
names(f)[names(f) == "tau"] <- "landmark_time"
f <- f %>% filter(approach == "conditional") %>% mutate(SAB = "Combined cohort")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(n_eff = 2641)
if (unscale){
  f <- f %>% mutate(est = ciu - 1.96* sqrt(var_est/n_eff))
}
combined_ss_both_conditional <- multisplit(f, f1, method)
combined_ss_both_conditional <- combined_ss_both_conditional

### marginal female
setwd("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/data_analysis/female")
# f <- readRDS("analysis_112223.rds")
f <- readRDS("female_analysis_interactions_oldlib_smallersigma_subsample1750.rds")
names(f)[names(f) == "tau"] <- "landmark_time"
f <- f %>% filter(approach == "marginal") %>% mutate(SAB = "Female cohort")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(n_eff = 1846)
if (unscale){
  f <- f %>% mutate(est = ciu - 1.96* sqrt(var_est/n_eff))
}
female_ss_both_marginal <- multisplit(f, f1, method)
female_ss_both_marginal <- female_ss_both_marginal
#
# ### female conditional
# f <- readRDS("analysis_112223.rds")
f <- readRDS("female_analysis_interactions_oldlib_smallersigma_subsample1750.rds")
names(f)[names(f) == "tau"] <- "landmark_time"
f <- f %>% filter(approach == "conditional") %>% mutate(SAB = "Female cohort")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(n_eff = 1846)
if (unscale){
  f <- f %>% mutate(est = ciu - 1.96* sqrt(var_est/n_eff))
}
female_ss_both_conditional <- multisplit(f, f1, method)
female_ss_both_conditional <- female_ss_both_conditional

### pooled male and female marginal
# setwd("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/data_analysis/male")
# # f <- readRDS("analysis_112223.rds")
# f <- readRDS("pooled_analysis_interactions_stackGoldlib.rds")
# names(f)[names(f) == "tau"] <- "landmark_time"
# f <- f %>% filter(approach == "marginal")
# f1 <- f[f$seed == unique(f$seed)[1],]
# f <- f %>% mutate(n_eff = ifelse(SAB == "male", 795, 1846))
# if (unscale){
#   f <- f %>% mutate(est = ciu - 1.96* sqrt(var_est/n_eff))
# }
# pooled_ss_both_marginal <- multisplit(f, f1, method)
# pooled_ss_both_marginal <- pooled_ss_both_marginal %>% mutate(SAB = ifelse(SAB == "male", "Male cohort", "Female cohort"))
#
# setwd("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/data_analysis/male")
# # f <- readRDS("analysis_112223.rds")
# f <- readRDS("pooled_analysis_interactions_stackGoldlib.rds")
# names(f)[names(f) == "tau"] <- "landmark_time"
# f <- f %>% filter(approach == "conditional")
# f1 <- f[f$seed == unique(f$seed)[1],]
# f <- f %>% mutate(n_eff = ifelse(SAB == "male", 795, 1846))
# if (unscale){
#   f <- f %>% mutate(est = ciu - 1.96* sqrt(var_est/n_eff))
# }
# pooled_ss_both_conditional <- multisplit(f, f1, method)
# pooled_ss_both_conditional <- pooled_ss_both_conditional %>% mutate(SAB = ifelse(SAB == "male", "Male cohort", "Female cohort"))

#

### marginal male
setwd("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/data_analysis/male")
# f <- readRDS("analysis_112223.rds")
f <- readRDS("male_analysis_interactions_oldlib_smallersigma_subsample1750_oldseed.rds")
names(f)[names(f) == "tau"] <- "landmark_time"
f <- f %>% filter(approach == "marginal") %>% mutate(SAB = "Male cohort")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(n_eff = 795)
if (unscale){
  f <- f %>% mutate(est = ciu - 1.96* sqrt(var_est/n_eff))
}
male_ss_both_marginal <- multisplit(f, f1, method)
male_ss_both_marginal <- male_ss_both_marginal

### male conditional
# f <- readRDS("analysis_112223.rds")
f <- readRDS("male_analysis_interactions_oldlib_smallersigma_subsample1750_oldseed.rds")
names(f)[names(f) == "tau"] <- "landmark_time"
f <- f %>% filter(approach == "conditional")%>% mutate(SAB = "Male cohort")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(n_eff = 795)
if (unscale){
  f <- f %>% mutate(est = ciu - 1.96* sqrt(var_est/n_eff))
}
male_ss_both_conditional  <- multisplit(f, f1, method)
male_ss_both_conditional  <- male_ss_both_conditional

######################
#### MAKE FIGURES ####
######################

make_plot_combined <- function(combined, type){
  xlab <- "Estimated variable importance"
  ylab <- "Variable group"
  
  xmax_auc <- 0.5
  xmin <- 0
  
  for (i in 1:nrow(combined)){
    if (combined$pval[i] <= 0.05){
      combined$indx_number[i] <- paste0("*", combined$indx_number[i])
    }
  }
  
  
  vimp_plot_combined <- combined %>%
    filter(SAB == "Combined cohort") %>%
    mutate(vim = case_when(vim == "AUC" & landmark_time == 545 ~ "AUC at 18 mo.",
                           vim == "AUC" & landmark_time == 730 ~ "AUC at 24 mo.",
                           vim == "AUC" & landmark_time == 912 ~ "AUC at 30 mo.",
                           vim == "cindex" ~ "C-index")) %>%
    arrange(vim, est) %>%
    mutate(Order = row_number()) %>%
    {ggplot(., aes(x = est, y = Order)) +
        geom_errorbarh(aes(xmin = cil, xmax = ciu), height = 0.0) +
        geom_point(size = point_size) +
        theme_bw() +
        xlab(xlab) +
        ylab(ylab) +
        xlim(c(xmin,xmax_auc)) +
        scale_y_continuous(
          breaks = .$Order,
          labels = .$indx_number,
          expand = c(ifelse(type == "marginal", 0.15, 0.1) ,0)
        ) +
        facet_wrap(~vim, dir = "v", strip.position = "right", scales = "free_y", nrow = 4) +
        ggtitle("Combined cohort")+
        theme(title = element_text(size = title_text_size, family = "Times New Roman"),
              axis.title.x = element_text(size = axis_title_size, family = "Times New Roman"),
              axis.title.y = element_text(size = axis_title_size, margin = margin(t = 0, r = 7, b = 0, l = 0)),
              axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank())
    }
  
  vimp_plot_female <- combined %>%
    filter(SAB == "Female cohort") %>%
    mutate(vim = case_when(vim == "AUC" & landmark_time == 545 ~ "AUC at 18 mo.",
                           vim == "AUC" & landmark_time == 730 ~ "AUC at 24 mo.",
                           vim == "AUC" & landmark_time == 912 ~ "AUC at 30 mo.",
                           vim == "cindex" ~ "C-index")) %>%
    arrange(vim, est) %>%
    mutate(Order = row_number()) %>%
    {ggplot(., aes(x = est, y = Order)) +
        geom_errorbarh(aes(xmin = cil, xmax = ciu), height = 0.0) +
        geom_point(size = point_size) +
        theme_bw() +
        xlab(xlab) +
        ylab(ylab) +
        xlim(c(xmin,xmax_auc)) +
        scale_y_continuous(
          breaks = .$Order,
          labels = .$indx_number,
          expand = c(ifelse(type == "marginal", 0.2,0.12), 0)
        ) +
        facet_wrap(~vim, dir = "v", strip.position = "right", scales = "free_y", nrow = 4) +
        ggtitle("Female cohort")+
        theme(title = element_text(size = title_text_size, family = "Times New Roman"),
              axis.title.x = element_text(size = axis_title_size, family = "Times New Roman"),
              axis.title.y = element_blank(),#element_text(size = axis_text_size, margin = margin(t = 0, r = 7, b = 0, l = 0)),
              axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank())
    }
  
  vimp_plot_male <- combined %>%
    filter(SAB == "Male cohort") %>%
    mutate(vim = case_when(vim == "AUC" & landmark_time == 545 ~ "AUC at 18 mo.",
                           vim == "AUC" & landmark_time == 730 ~ "AUC at 24 mo.",
                           vim == "AUC" & landmark_time == 912 ~ "AUC at 30 mo.",
                           vim == "cindex" ~ "C-index")) %>%
    arrange(vim, est) %>%
    mutate(Order = row_number()) %>%
    {ggplot(., aes(x = est, y = Order)) +
        geom_errorbarh(aes(xmin = cil, xmax = ciu), height = 0.0) +
        geom_point(size = point_size) +
        theme_bw() +
        xlab(xlab) +
        ylab(ylab) +
        xlim(c(xmin,xmax_auc)) +
        scale_y_continuous(
          breaks = .$Order,
          labels = .$indx_number,
          expand = c(ifelse(type == "marginal", 0.2,0.12), 0)
        ) +
        facet_wrap(~vim, dir = "v", strip.position = "right", scales = "free_y", nrow = 4) +
        ggtitle("Male cohort")+
        theme(title = element_text(size = title_text_size, family = "Times New Roman"),
              axis.title.x = element_text(size = axis_title_size, family = "Times New Roman"),
              axis.title.y = element_blank(),#element_text(size = axis_text_size, margin = margin(t = 0, r = 7, b = 0, l = 0)),
              axis.text = element_text(size = axis_text_size, family = "Times New Roman"),
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_text(size = strip_text_size, family = "Times New Roman"),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank())
    }
  
  multipanel <- ggarrange(vimp_plot_combined,
                          vimp_plot_female,
                          vimp_plot_male,
                          nrow = 1,
                          ncol = 3,
                          widths = c(0.34, 0.32, 0.35))
  return(multipanel)
}

# combine results and make plots
combined_ss_marginal <- bind_rows(combined_ss_both_marginal,
                                  # pooled_ss_both_marginal
                                  male_ss_both_marginal,
                                  female_ss_both_marginal) %>%
  mutate(SAB = factor(SAB, levels = c("Combined cohort", "Female cohort", "Male cohort")))

# FOR ENAR
# combined_ss_marginal <- combined_ss_marginal %>% filter(SAB != "Female cohort")

p_marg <- make_plot_combined(combined_ss_marginal, type = "marginal")


wd <- "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/biometrika/"
fname <- "702-marginal-vax-072124-10split-tnr-sub1750"
ggsave(filename = paste0(wd, fname, ".pdf"),
       plot = p_marg, device = "pdf",
       width = big_fig_width, height = big_fig_height, dpi = 300, units = "in")

combined_ss_conditional <- bind_rows(combined_ss_both_conditional,
                                     # pooled_ss_both_conditional) %>%
                                     female_ss_both_conditional,
                                     male_ss_both_conditional) %>%
  mutate(SAB = factor(SAB, levels = c("Combined cohort", "Female cohort", "Male cohort")))

p_cond <- make_plot_combined(combined_ss_conditional, type = "conditional")

fname <- "702-conditional-vax-072124-10split-tnr-sub1750"
ggsave(filename = paste0(wd, fname, ".pdf"),
       plot = p_cond, device = "pdf",
       width = big_fig_width, height = big_fig_height, dpi = 300, units = "in")


# look at significance
sig_marginal <- combined_ss_marginal %>%
  filter(pval < 0.05)
sig_conditional <- combined_ss_conditional %>%
  filter(pval < 0.05)

