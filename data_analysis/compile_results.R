#!/usr/local/bin/Rscript

# packages
library("cowplot")
library("grid")
library("gridExtra")
library("ggpubr")
library("tidyverse")

# plotting parameters
point_size <- 3
title_text_size <- 14
axis_text_size <- 12
big_fig_width <- 12
big_fig_height <- 8
strip_text_size <- 10

########################################
### combine results from multiple splits
########################################

multisplit <- function(f, f1){
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
        z <- (dat$estimate[k] - ci_grid[j]) / sqrt(dat$var_est[k]/dat$n_eff[k])
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
      p_two <- min(p.adjust(p_array_twosided[k,j,],method="bonferroni"))
      p_one <- min(p.adjust(p_array_onesided[k,j,],method="bonferroni"))
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
  
  f_avg <- f %>% group_by(tau, indx_name, vim) %>% summarize(estimate = mean(estimate)) %>%
    mutate(estimate = ifelse(estimate >= 0, estimate, 0))
  
  combined_ss <- f1 %>% select(tau, vim, indx, indx_name) %>%
    mutate(cil = cil, ciu = ciu, ci_1sided = ci_1sided, pval = pvals)
  
  combined_ss <- left_join(combined_ss, f_avg, by = c("tau", "indx_name", "vim"))
  
  combined_ss <- combined_ss %>% 
    filter(indx_name != "sex_bmi" & indx_name != "all_but_geo" & indx_name != "BRS") %>%
    mutate(cil = ifelse(cil >= 0, cil, 0),
           indx_number = case_when(
             indx_name == "sex" ~ "1",
             indx_name == "age" ~ "2",
             indx_name == "bmi" ~ "3",
             indx_name == "BRS_sexhealth" ~ "4",
             indx_name == "BRS_sexbehave" ~ "5",
             indx_name == "BRS_social" ~ "6",
             indx_name == "geo" ~ "7"
           ))
  
  return(combined_ss)
}

### marginal combined
setwd("/home/cwolock/surv_vim_supplementary/data_analysis/combined")
f <- readRDS("analysis_112223.rds")
f <- f %>% filter(approach == "marginal")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(estimate = one_step, n_eff = 2641)
combined_ss_both_marginal <- multisplit(f, f1)
combined_ss_both_marginal <- combined_ss_both_marginal %>% mutate(SAB = "Combined cohort")

### conditional combined
f <- readRDS("analysis_112223.rds")
f <- f %>% filter(approach == "conditional")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(estimate = one_step, n_eff = 2641)
combined_ss_both_conditional <- multisplit(f, f1)
combined_ss_both_conditional <- combined_ss_both_conditional %>% mutate(SAB = "Combined cohort")

### marginal female
setwd("/home/cwolock/surv_vim_supplementary/data_analysis/female")
f <- readRDS("analysis_112223.rds")
f <- f %>% filter(approach == "marginal")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(estimate = one_step, n_eff = 1846)
female_ss_both_marginal <- multisplit(f, f1)
female_ss_both_marginal <- female_ss_both_marginal %>% mutate(SAB = "Female cohort")

### female conditional
f <- readRDS("analysis_112223.rds")
f <- f %>% filter(approach == "conditional")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(estimate = one_step, n_eff = 1846)
female_ss_both_conditional <- multisplit(f, f1)
female_ss_both_conditional <- female_ss_both_conditional %>% mutate(SAB = "Female cohort")

### marginal male
setwd("/home/cwolock/surv_vim_supplementary/data_analysis/male")
f <- readRDS("analysis_112223.rds")
f <- f %>% filter(approach == "marginal")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(estimate = one_step, n_eff = 795)
male_ss_both_marginal <- multisplit(f, f1)
male_ss_both_marginal <- male_ss_both_marginal %>% mutate(SAB = "Male cohort")

### male conditional
f <- readRDS("analysis_112223.rds")
f <- f %>% filter(approach == "conditional")
f1 <- f[f$seed == unique(f$seed)[1],]
f <- f %>% mutate(estimate = one_step, n_eff = 795)
male_ss_both_conditional  <- multisplit(f, f1)
male_ss_both_conditional  <- male_ss_both_conditional  %>% mutate(SAB = "Male cohort")

######################
#### MAKE FIGURES ####
######################

make_plot_combined <- function(combined){
  xlab <- "Estimated variable importance"
  ylab <- "Variable group"
  
  xmax_auc <- 0.5
  xmin <- 0
  
  vimp_plot_auc_545 <- combined %>%
    filter(tau == 545 & vim == "AUC") %>%
    arrange(SAB, estimate) %>%
    mutate(Order = row_number()) %>%
    {ggplot(., aes(x = estimate, y = Order)) +
        geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
        geom_point(size = point_size) +
        theme_bw() +
        xlab(xlab) +
        ylab(ylab) +
        xlim(c(xmin,xmax_auc)) +
        scale_y_continuous(
          breaks = .$Order,
          labels = .$indx_number,
        ) + 
        facet_wrap(~SAB, dir = "v", strip.position = "right", scales = "free_y") + 
        ggtitle("AUC (18 mo. follow-up)")+
        theme(title = element_text(size = title_text_size),
              axis.title = element_text(size = axis_text_size),
              axis.text = element_text(size = axis_text_size),
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_blank())
    }
  
  vimp_plot_auc_730 <- combined %>%
    filter(tau == 730 & vim == "AUC") %>%
    arrange(SAB, estimate) %>%
    mutate(Order = row_number()) %>%
    {ggplot(., aes(x = estimate, y = Order)) +
        geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
        geom_point(size = point_size) +
        theme_bw() +
        xlab(xlab) +
        ylab(ylab) +
        xlim(c(xmin,xmax_auc)) +
        scale_y_continuous(
          breaks = .$Order,
          labels = .$indx_number,
        ) + 
        facet_wrap(~SAB, dir = "v", strip.position = "right", scales = "free_y") + 
        ggtitle("AUC (24 mo. follow-up)")+
        theme(title = element_text(size = title_text_size),
              axis.title.y = element_blank(), 
              axis.title = element_text(size = axis_text_size),
              axis.text = element_text(size = axis_text_size),
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_blank())
    }
  
  vimp_plot_auc_912 <- combined %>%
    filter(tau == 912 & vim == "AUC") %>%
    arrange(SAB, estimate) %>%
    mutate(Order = row_number()) %>%
    {ggplot(., aes(x = estimate, y = Order)) +
        geom_errorbarh(aes(xmin = cil, xmax = ciu)) +
        geom_point(size = point_size) +
        theme_bw() +
        xlab(xlab) +
        ylab(ylab) +
        xlim(c(xmin,xmax_auc)) +
        scale_y_continuous(
          breaks = .$Order,
          labels = .$indx_number,
        ) + 
        facet_wrap(~SAB, dir = "v", strip.position = "right", scales = "free_y") + 
        ggtitle("AUC (30 mo. follow-up)")+
        theme(title = element_text(size = title_text_size),
              axis.title.y = element_blank(), 
              axis.title = element_text(size = axis_text_size),
              axis.text = element_text(size = axis_text_size),
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_text(size = strip_text_size))
    }
  
  multipanel <- ggarrange(vimp_plot_auc_545,
                          vimp_plot_auc_730,
                          vimp_plot_auc_912,
                          nrow = 1,
                          ncol = 3)
  return(multipanel)
}

# combine results and make plots
combined_ss_marginal <- bind_rows(combined_ss_both_marginal, 
                                  combined_ss_male_marginal, 
                                  combined_ss_female_marginal) %>%
  mutate(SAB = factor(SAB, levels = c("Combined cohort", "Female cohort", "Male cohort")))

p <- make_plot_combined(combined_ss_marginal)


wd <- "/home/cwolock/surv_vim_supplementary/scratch/biometrika/"
fname <- "702_marginal_vax_112223_10split"
ggsave(filename = paste0(wd, fname, ".eps"),
       plot = p, device = "eps",
       width = big_fig_width, height = big_fig_height, dpi = 300, units = "in")

combined_ss_conditional <- bind_rows(combined_ss_both_conditional, 
                                     combined_ss_male_conditional, 
                                     combined_ss_female_conditional) %>%
  mutate(SAB = factor(SAB, levels = c("Combined cohort", "Female cohort", "Male cohort")))

p <- make_plot_combined(combined_ss_conditional)

fname <- "702_conditional_vax_112223_10split"
ggsave(filename = paste0(wd, fname, ".pdf"),
       plot = p, device = "pdf",
       width = big_fig_width, height = big_fig_height, dpi = 300, units = "in")


# look at significance
sig_marginal <- combined_ss_marginal %>% 
  filter((pval < 0.05/5 & SAB != "Combined cohort") | (pval < 0.05/6 & SAB == "Combined cohort"))
sig_conditional <- combined_ss_conditional %>% 
  filter((pval < 0.05/6 & SAB != "Combined cohort") | (pval < 0.05/7 & SAB == "Combined cohort"))
