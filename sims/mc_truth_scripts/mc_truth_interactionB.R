library(tidyverse)
# source("/home/cwolock/surv_vim_supplementary/sims/generate_data.R")
source("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/sims/generate_data.R")
n_train <- 2e7
set.seed(1234)
sdy <- 1
##################
### no correlation
##################
p <- 5
dat <- generate_data(n = n_train, scenario = "2", sdy = sdy)
x <- as.matrix(dat[,1:p])
beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
beta_int <- c(0.1, -0.1, 0.1)
interceptt <- 0
landmark_times <- c(0.5, 0.9)
full_brier <- rep(NA, length(landmark_times))
brier_01 <- rep(NA, length(landmark_times))
brier_02 <- rep(NA, length(landmark_times))
full_rsquared <- rep(NA, length(landmark_times))
rsquared_01 <- rep(NA, length(landmark_times))
rsquared_02 <- rep(NA, length(landmark_times))
full_auc <- rep(NA, length(landmark_times))
auc_01 <- rep(NA, length(landmark_times))
auc_02 <- rep(NA, length(landmark_times))

for (t in landmark_times){
  y <- ifelse(dat$t > t, 1, 0)
  f_0 <- 1-pnorm(-(x %*% beta_t + x[,1]*x[,2]*beta_int[1] + x[,3]*x[,4]*beta_int[2]+x[,1]*x[,5]*beta_int[3]) -interceptt + log(t),
                 mean = 0,
                 sd = sdy)
  f_01 <- 1-pnorm(-(x[, 2] * beta_t[2] + x[,3]*x[,4]*beta_int[2] ) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt((beta_t[1] + beta_int[1]*x[,2] + beta_int[3]*x[,5])^2 + sdy^2))
  f_02 <- 1-pnorm(-(x[, 1] * beta_t[1] + x[,3]*x[,4]*beta_int[2] + x[,1]*x[,5]*beta_int[3]) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt((beta_t[2] + beta_int[1]*x[,1])^2 + sdy^2))

  mse_t <- vimp::measure_mse(f_0, y)$point_est
  auc_t <- cvAUC::AUC(f_0, y)
  mse_t1 <- vimp::measure_mse(f_01, y)$point_est
  mse_t2 <- vimp::measure_mse(f_02, y)$point_est
  auc_t1 <- cvAUC::AUC(f_01, y)
  auc_t2 <- cvAUC::AUC(f_02, y)

  full_brier[which(landmark_times == t)] <- -mse_t
  brier_01[which(landmark_times == t)] <- -mse_t1
  brier_02[which(landmark_times == t)] <- -mse_t2
  full_rsquared[which(landmark_times == t)] <- 1-mse_t/var(y)
  rsquared_01[which(landmark_times == t)] <- 1-mse_t1/var(y)
  rsquared_02[which(landmark_times == t)] <- 1-mse_t2/var(y)
  full_auc[which(landmark_times == t)] <- auc_t
  auc_01[which(landmark_times == t)] <- auc_t1
  auc_02[which(landmark_times == t)] <- auc_t2
}

output1_landmark <- data.frame(vim = rep(c("brier", "AUC", "rsquared"), each = length(landmark_times)),
                               tau = c(landmark_times, landmark_times, landmark_times),
                               V_full = c(full_brier, full_auc, full_rsquared),
                               V_01 = c(brier_01, auc_01, rsquared_01),
                               V_02 = c(brier_02, auc_02, rsquared_02),
                               V_06 = c(full_brier, full_auc, full_rsquared),
                               n_mc = rep(n_train, 3*length(landmark_times)))
output1_landmark <- output1_landmark %>% mutate(V_016 = V_01, correlation = FALSE)

taus <- c(0.9)
# c index stuff
dat_test <- generate_data(n = n_train, scenario = "2")
outcome <- dat_test$t
preds <- -beta_t[1]*dat_test[,1] - beta_t[2]*dat_test[,2] - beta_int[1]*dat_test[,1]*dat_test[,2] - beta_int[2]*dat_test[,3]*dat_test[,4] - beta_int[3]*dat_test[,1]*dat_test[,5]
preds_01 <- -beta_t[2]*dat_test[,2] - beta_int[2]*dat_test[,3]*dat_test[,4]
preds_02 <- -beta_t[1]*dat_test[,1] - beta_int[2]*dat_test[,3]*dat_test[,4] - beta_int[3]*dat_test[,1]*dat_test[,5]
outcome2 <- dat$t
preds2 <- -beta_t[1]*dat[,1] - beta_t[2]*dat[,2]- beta_int[1]*dat[,1]*dat[,2] - beta_int[2]*dat[,3]*dat[,4] - beta_int[3]*dat[,1]*dat[,5]
preds2_01 <- -beta_t[2]*dat[,2] - beta_int[2]*dat[,3]*dat[,4]
preds2_02 <- -beta_t[1]*dat[,1] - beta_int[2]*dat[,3]*dat[,4] - beta_int[3]*dat[,1]*dat[,5]

full_c <- rep(NA, length(taus))
c_01 <- rep(NA, length(taus))
c_02 <- rep(NA, length(taus))
for (tau in taus){
  denom <- mean((outcome < outcome2 & outcome <= tau) | (outcome2 < outcome & outcome2 <= tau))

  full_c[which(taus == tau)] <- (mean(ifelse(preds > preds2, 1, 0)*(outcome < outcome2 & outcome <= tau)) +
                                   mean(ifelse(preds2 > preds, 1, 0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_01[which(taus == tau)] <- (mean(ifelse(preds_01 > preds2_01, 1, 0)*(outcome < outcome2 & outcome <= tau)) +
                                 mean(ifelse(preds2_01 > preds_01, 1, 0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_02[which(taus == tau)] <- (mean(ifelse(preds_02 > preds2_02,1,0)*(outcome < outcome2 & outcome <= tau)) +
                                 mean(ifelse(preds2_02 > preds_02,1,0)*(outcome2 < outcome & outcome2 <= tau))) / denom

}

output1_global <- data.frame(vim = "cindex",
                             tau = taus,
                             V_full = full_c,
                             V_01 = c_01,
                             V_02 = c_02,
                             V_06 = full_c,
                             n_mc = n_train)
output1_global <- output1_global %>% mutate(V_016 = V_01, correlation = FALSE)
output1 <- bind_rows(output1_landmark, output1_global)

###############
### correlation
###############
dat <- generate_data(n = n_train, scenario = "1")
p <- 6
x <- as.matrix(dat[,1:p])
beta_t <- matrix(c(0.5, -0.3, rep(0, p- 2)))
beta_int <- c(0.1, -0.1, 0.1)
rho16 <- 0.7
rho23 <- -0.3
interceptt <- 0

full_brier <- rep(NA, length(landmark_times))
brier_01 <- rep(NA, length(landmark_times))
brier_02 <- rep(NA, length(landmark_times))
full_rsquared <- rep(NA, length(landmark_times))
rsquared_01 <- rep(NA, length(landmark_times))
rsquared_02 <- rep(NA, length(landmark_times))
full_auc <- rep(NA, length(landmark_times))
auc_01 <- rep(NA, length(landmark_times))
auc_02 <- rep(NA, length(landmark_times))
brier_016 <- rep(NA, length(landmark_times))
rsquared_016 <- rep(NA, length(landmark_times))
auc_016 <- rep(NA, length(landmark_times))

for (t in landmark_times){
  y <- ifelse(dat$t > t, 1, 0)
  f_0 <- 1-pnorm(-(x %*% beta_t + x[,1]*x[,2]*beta_int[1] + x[,3]*x[,4]*beta_int[2] + x[,1]*x[,5]*beta_int[3]) -interceptt + log(t),
                 mean = 0,
                 sd = sdy)
  f_01 <- 1-pnorm(-(x[,2]*beta_t[2] + beta_int[2]*x[,3]*x[,4]) - interceptt + log(t),
                  mean = (beta_t[1] + beta_int[1]*x[,2] + beta_int[3]*x[,5])*rho16*x[,6],
                  sd = sqrt((beta_t[1] + beta_int[1]*x[,2] + beta_int[3]*x[,5])^2*(1-rho16^2) + sdy^2))
  f_02 <- 1-pnorm(-(x[,1]*beta_t[1] + beta_int[2]*x[,3]*x[,4] + beta_int[3]*x[,1]*x[,5]) - interceptt + log(t),
                  mean = (beta_t[2] + beta_int[1]*x[,1])*rho23*x[,3],
                  sd = sqrt((beta_t[2] + beta_int[1]*x[,1])^2*(1-rho23^2) + sdy^2))
  f_016 <- 1-pnorm(-(x[,2]*beta_t[2] + beta_int[2]*x[,3]*x[,4]) - interceptt + log(t),
                   mean = 0,
                   sd = sqrt((beta_t[1] + beta_int[1]*x[,2] + beta_int[3]*x[,5])^2 + sdy^2))
  mse_t <- vimp::measure_mse(f_0, y)$point_est
  auc_t <- cvAUC::AUC(f_0, y)
  mse_t1 <- vimp::measure_mse(f_01, y)$point_est
  mse_t2 <- vimp::measure_mse(f_02, y)$point_est
  auc_t1 <- cvAUC::AUC(f_01, y)
  auc_t2 <- cvAUC::AUC(f_02, y)

  mse_t16 <- vimp::measure_mse(f_016, y)$point_est
  auc_t16 <- cvAUC::AUC(f_016, y)

  full_brier[which(landmark_times == t)] <- -mse_t
  brier_01[which(landmark_times == t)] <- -mse_t1
  brier_02[which(landmark_times == t)] <- -mse_t2
  full_rsquared[which(landmark_times == t)] <- 1-mse_t/var(y)
  rsquared_01[which(landmark_times == t)] <- 1-mse_t1/var(y)
  rsquared_02[which(landmark_times == t)] <- 1-mse_t2/var(y)
  full_auc[which(landmark_times == t)] <- auc_t
  auc_01[which(landmark_times == t)] <- auc_t1
  auc_02[which(landmark_times == t)] <- auc_t2

  brier_016[which(landmark_times == t)] <- -mse_t16
  rsquared_016[which(landmark_times == t)] <- 1-mse_t16/var(y)
  auc_016[which(landmark_times == t)] <- auc_t16
}

output2_landmark <- data.frame(vim = rep(c("brier", "AUC", "rsquared"), each = length(landmark_times)),
                               tau = c(landmark_times, landmark_times, landmark_times),
                               V_full = c(full_brier, full_auc, full_rsquared),
                               V_01 = c(brier_01, auc_01, rsquared_01),
                               V_02 = c(brier_02, auc_02, rsquared_02),
                               V_016 = c(brier_016, auc_016, rsquared_016),
                               V_06 = c(full_brier, full_auc, full_rsquared),
                               n_mc = rep(n_train, 3*length(landmark_times)))
output2_landmark <- output2_landmark %>% mutate(correlation = TRUE)

dat_test <- generate_data(n = n_train, scenario = "1")
outcome <- dat_test$t
preds <- -beta_t[1]*dat_test[,1] -beta_t[2]*dat_test[,2]- beta_int[1]*dat_test[,1]*dat_test[,2] - beta_int[2]*dat_test[,3]*dat_test[,4] - beta_int[3]*dat_test[,1]*dat_test[,5]
preds_01 <- -beta_t[2]*dat_test[,2] - (beta_t[1] + beta_int[1]*dat_test[,2] + beta_int[3]*dat_test[,5])*rho16*dat_test[,6] - beta_int[2]*dat_test[,3]*dat_test[,4]
preds_02 <- -beta_t[1]*dat_test[,1] - (beta_t[2] + beta_int[1]*dat_test[,1])*rho23*dat_test[,3]- beta_int[2]*dat_test[,3]*dat_test[,4] - beta_int[3]*dat_test[,1]*dat_test[,5]
preds_016 <- -beta_t[2]*dat_test[,2] - beta_int[2]*dat_test[,3]*dat_test[,4]
outcome2 <- dat$t
preds2 <- -beta_t[1]*dat[,1] -beta_t[2]*dat[,2]- beta_int[1]*dat[,1]*dat[,2] - beta_int[2]*dat[,3]*dat[,4] - beta_int[3]*dat[,1]*dat[,5]
preds2_01 <- -beta_t[2]*dat[,2] - (beta_t[1] + beta_int[1]*dat[,2] + beta_int[3]*dat[,5])*rho16*dat[,6]- beta_int[2]*dat[,3]*dat[,4]
preds2_02 <- -beta_t[1]*dat[,1] -(beta_t[2] + beta_int[1]*dat[,1])*rho23*dat[,3] - beta_int[2]*dat[,3]*dat[,4] - beta_int[3]*dat[,1]*dat[,5]
preds2_016 <- -beta_t[2]*dat[,2]- beta_int[2]*dat[,3]*dat[,4]

full_c <- rep(NA, length(taus))
c_01 <- rep(NA, length(taus))
c_02 <- rep(NA, length(taus))
c_016 <- rep(NA, length(taus))

for (tau in taus){
  denom <- mean((outcome < outcome2 & outcome <= tau) | (outcome2 < outcome & outcome2 <= tau))

  full_c[which(taus == tau)] <- (mean(ifelse(preds > preds2, 1, 0)*(outcome < outcome2 & outcome <= tau)) +
                                   mean(ifelse(preds2 > preds, 1, 0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_01[which(taus == tau)]  <- (mean(ifelse(preds_01 > preds2_01, 1, 0)*(outcome < outcome2 & outcome <= tau)) +
                                  mean(ifelse(preds2_01 > preds_01, 1, 0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_02[which(taus == tau)]  <- (mean(ifelse(preds_02 > preds2_02,1,0)*(outcome < outcome2 & outcome <= tau)) +
                                  mean(ifelse(preds2_02 > preds_02,1,0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_016[which(taus == tau)]  <- (mean(ifelse(preds_016 > preds2_016, 1, 0)*(outcome < outcome2 & outcome <= tau)) +
                                   mean(ifelse(preds2_016 > preds_016, 1, 0)*(outcome2 < outcome & outcome2 <= tau))) / denom
}

output2_global <- data.frame(vim = "cindex",
                             tau = taus,
                             V_full = full_c,
                             V_01 = c_01,
                             V_02 = c_02,
                             V_016 = c_016,
                             V_06 = full_c,
                             n_mc = n_train)
output2_global <- output2_global %>% mutate(correlation = TRUE)
output2 <- bind_rows(output2_landmark, output2_global)

output <- bind_rows(output1, output2)
saveRDS(output, "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/truth_interactionB.rds")
