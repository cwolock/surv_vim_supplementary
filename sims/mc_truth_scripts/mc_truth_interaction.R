library(tidyverse)
# source("/home/cwolock/surv_vim_supplementary/sims/generate_data.R")
source("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/generate_data.R")
n_train <- 1e6
set.seed(1234)
sdy <- 1
##################
### no correlation
##################
p <- 4
dat <- generate_data(n = n_train, scenario = "1A", sdy = sdy)
x <- as.matrix(dat[,1:p])
beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
beta_int1 <- 0.2
beta_int2 <- -0.1
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
  f_0 <- 1-pnorm(-(x %*% beta_t + x[,1]*x[,2]*beta_int1 + x[,3]*x[,4]*beta_int2) -interceptt + log(t),
                 mean = 0,
                 sd = sdy)
  f_01 <- 1-pnorm(-(x[, 2] * beta_t[2] + x[,3]*x[,4]*beta_int2) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt((beta_t[1] + beta_int1*x[,2])^2 + sdy^2))
  f_02 <- 1-pnorm(-(x[, 1] * beta_t[1] + x[,3]*x[,4]*beta_int2) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt((beta_t[2] + beta_int1*x[,1])^2 + sdy^2))

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
                               V_05 = c(full_brier, full_auc, full_rsquared),
                               n_mc = rep(n_train, 3*length(landmark_times)))
output1_landmark <- output1_landmark %>% mutate(V_015 = V_01, correlation = FALSE)

taus <- c(0.9)
# approx_times <- seq(0, max(taus), by = .01)
# S_0 <- matrix(NA, nrow = n_train, ncol = length(approx_times))
# S_01 <- matrix(NA, nrow = n_train, ncol = length(approx_times))
# S_02 <- matrix(NA, nrow = n_train, ncol = length(approx_times))
# for (t in approx_times){
#   S_0[,which(approx_times == t)] <- 1-pnorm(-(x %*% beta_t) -interceptt + log(t), sd = sdy)
#   S_01[,which(approx_times == t)] <- 1-pnorm(-(x[, 2] * beta_t[2]) - interceptt + log(t), sd = sqrt(beta_t[1]^2 + sdy^2))
#   S_02[,which(approx_times == t)] <- 1-pnorm(-(x[, 1] * beta_t[1]) - interceptt + log(t), sd = sqrt(beta_t[2]^2 + sdy^2))
# }
# c index stuff
dat_test <- generate_data(n = n_train, scenario = "1A")
outcome <- dat_test$t
preds <- -beta_t[1]*dat_test[,1] - beta_t[2]*dat_test[,2] - beta_int1*dat_test[,1]*dat_test[,2] - beta_int2*dat_test[,3]*dat_test[,4]
preds_01 <- -beta_t[2]*dat_test[,2] - beta_int2*dat_test[,3]*dat_test[,4]
preds_02 <- -beta_t[1]*dat_test[,1] - beta_int2*dat_test[,3]*dat_test[,4]
outcome2 <- dat$t
preds2 <- -beta_t[1]*dat[,1] - beta_t[2]*dat[,2]- beta_int1*dat[,1]*dat[,2] - beta_int2*dat[,3]*dat[,4]
preds2_01 <- -beta_t[2]*dat[,2] - beta_int2*dat[,3]*dat[,4]
preds2_02 <- -beta_t[1]*dat[,1] - beta_int2*dat[,3]*dat[,4]

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

  # f_0 <- calc_rmst(S_0, tau, approx_times)
  # f_01 <- calc_rmst(S_01, tau, approx_times)
  # f_02 <- calc_rmst(S_02, tau, approx_times)
  # f_03 <- f_04 <- f_0
}

output1_global <- data.frame(vim = "cindex",
                             tau = taus,
                             V_full = full_c,
                             V_01 = c_01,
                             V_02 = c_02,
                             V_05 = full_c,
                             n_mc = n_train)
output1_global <- output1_global %>% mutate(V_015 = V_01, correlation = FALSE)
output1 <- bind_rows(output1_landmark, output1_global)

###############
### correlation
###############
dat <- generate_data(n = n_train, scenario = "4A")
p <- 5
x <- as.matrix(dat[,1:p])
beta_t <- matrix(c(0.5, -0.3, rep(0, p- 2)))
beta_int1 <- 0.2
beta_int2 <- -0.1
rho15 <- 0.7
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
brier_015 <- rep(NA, length(landmark_times))
rsquared_015 <- rep(NA, length(landmark_times))
auc_015 <- rep(NA, length(landmark_times))

for (t in landmark_times){
  y <- ifelse(dat$t > t, 1, 0)
  f_0 <- 1-pnorm(-(x %*% beta_t + x[,1]*x[,2]*beta_int1 + x[,3]*x[,4]*beta_int2) -interceptt + log(t),
                 mean = 0,
                 sd = sdy)
  f_01 <- 1-pnorm(-(x[,2]*beta_t[2] + beta_int2*x[,3]*x[,4]) - interceptt + log(t),
                  mean = (beta_t[1] + beta_int1*x[,2])*rho15*x[,5],
                  sd = sqrt((beta_t[1] + beta_int1*x[,2])^2*(1-rho15^2) + sdy^2))
  f_02 <- 1-pnorm(-(x[,1]*beta_t[1] + beta_int2*x[,3]*x[,4]) - interceptt + log(t),
                  mean = (beta_t[2] + beta_int1*x[,1])*rho23*x[,3],
                  sd = sqrt((beta_t[2] + beta_int1*x[,1])^2*(1-rho23^2) + sdy^2))
  f_015 <- 1-pnorm(-(x[,2]*beta_t[2] + beta_int2*x[,3]*x[,4]) - interceptt + log(t),
                   mean = 0,
                   sd = sqrt((beta_t[1] + beta_int1*x[,2])^2 + sdy^2))
  mse_t <- vimp::measure_mse(f_0, y)$point_est
  auc_t <- cvAUC::AUC(f_0, y)
  mse_t1 <- vimp::measure_mse(f_01, y)$point_est
  mse_t2 <- vimp::measure_mse(f_02, y)$point_est
  auc_t1 <- cvAUC::AUC(f_01, y)
  auc_t2 <- cvAUC::AUC(f_02, y)

  mse_t15 <- vimp::measure_mse(f_015, y)$point_est
  auc_t15 <- cvAUC::AUC(f_015, y)

  full_brier[which(landmark_times == t)] <- -mse_t
  brier_01[which(landmark_times == t)] <- -mse_t1
  brier_02[which(landmark_times == t)] <- -mse_t2
  full_rsquared[which(landmark_times == t)] <- 1-mse_t/var(y)
  rsquared_01[which(landmark_times == t)] <- 1-mse_t1/var(y)
  rsquared_02[which(landmark_times == t)] <- 1-mse_t2/var(y)
  full_auc[which(landmark_times == t)] <- auc_t
  auc_01[which(landmark_times == t)] <- auc_t1
  auc_02[which(landmark_times == t)] <- auc_t2

  brier_015[which(landmark_times == t)] <- -mse_t15
  rsquared_015[which(landmark_times == t)] <- 1-mse_t15/var(y)
  auc_015[which(landmark_times == t)] <- auc_t15
}

output2_landmark <- data.frame(vim = rep(c("brier", "AUC", "rsquared"), each = length(landmark_times)),
                               tau = c(landmark_times, landmark_times, landmark_times),
                               V_full = c(full_brier, full_auc, full_rsquared),
                               V_01 = c(brier_01, auc_01, rsquared_01),
                               V_02 = c(brier_02, auc_02, rsquared_02),
                               V_015 = c(brier_015, auc_015, rsquared_015),
                               V_05 = c(full_brier, full_auc, full_rsquared),
                               n_mc = rep(n_train, 3*length(landmark_times)))
output2_landmark <- output2_landmark %>% mutate(correlation = TRUE)

dat_test <- generate_data(n = n_train, scenario = "4A")
outcome <- dat_test$t
preds <- -beta_t[1]*dat_test[,1] -beta_t[2]*dat_test[,2]- beta_int1*dat_test[,1]*dat_test[,2] - beta_int2*dat_test[,3]*dat_test[,4]
preds_01 <- -beta_t[2]*dat_test[,2] - (beta_t[1] + beta_int1*dat_test[,2])*rho15*dat_test[,5] - beta_int2*dat_test[,3]*dat_test[,4]
preds_02 <- -beta_t[1]*dat_test[,1] - (beta_t[2] + beta_int1*dat_test[,1])*rho23*dat_test[,3]- beta_int2*dat_test[,3]*dat_test[,4]
preds_015 <- -beta_t[2]*dat_test[,2] - beta_int2*dat_test[,3]*dat_test[,4]
outcome2 <- dat$t
preds2 <- -beta_t[1]*dat[,1] -beta_t[2]*dat[,2]- beta_int1*dat[,1]*dat[,2] - beta_int2*dat[,3]*dat[,4]
preds2_01 <- -beta_t[2]*dat[,2] - (beta_t[1] + beta_int1*dat[,2])*rho15*dat[,5]- beta_int2*dat[,3]*dat[,4]
preds2_02 <- -beta_t[1]*dat[,1] -(beta_t[2] + beta_int1*dat[,1])*rho23*dat[,3] - beta_int2*dat[,3]*dat[,4]
preds2_015 <- -beta_t[2]*dat[,2]- beta_int2*dat[,3]*dat[,4]

full_c <- rep(NA, length(taus))
c_01 <- rep(NA, length(taus))
c_02 <- rep(NA, length(taus))
c_015 <- rep(NA, length(taus))

for (tau in taus){
  denom <- mean((outcome < outcome2 & outcome <= tau) | (outcome2 < outcome & outcome2 <= tau))

  full_c[which(taus == tau)] <- (mean(ifelse(preds > preds2, 1, 0)*(outcome < outcome2 & outcome <= tau)) +
                                   mean(ifelse(preds2 > preds, 1, 0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_01[which(taus == tau)]  <- (mean(ifelse(preds_01 > preds2_01, 1, 0)*(outcome < outcome2 & outcome <= tau)) +
                                  mean(ifelse(preds2_01 > preds_01, 1, 0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_02[which(taus == tau)]  <- (mean(ifelse(preds_02 > preds2_02,1,0)*(outcome < outcome2 & outcome <= tau)) +
                                  mean(ifelse(preds2_02 > preds_02,1,0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_015[which(taus == tau)]  <- (mean(ifelse(preds_015 > preds2_015, 1, 0)*(outcome < outcome2 & outcome <= tau)) +
                                   mean(ifelse(preds2_015 > preds_015, 1, 0)*(outcome2 < outcome & outcome2 <= tau))) / denom
}

output2_global <- data.frame(vim = "cindex",
                             tau = taus,
                             V_full = full_c,
                             V_01 = c_01,
                             V_02 = c_02,
                             V_015 = c_015,
                             V_05 = full_c,
                             n_mc = n_train)
output2_global <- output2_global %>% mutate(correlation = TRUE)
output2 <- bind_rows(output2_landmark, output2_global)

output2 <- bind_rows(output1, output2)
# saveRDS(output, "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/truth_interaction.rds")
