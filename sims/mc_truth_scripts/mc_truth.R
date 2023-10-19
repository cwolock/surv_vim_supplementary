library(tidyverse)
source("/home/cwolock/surv_vim_supplementary/sims/generate_data_AFT.R")
n_train <- 1e7
set.seed(1234)

##################
### no correlation
##################
dat <- generate_data(n = n_train, scenario = "1", sdy = 1)
x <- as.matrix(dat[,1:2])
beta_t <- matrix(c(0.5, -0.3))
interceptt <- 0
landmark_times <- c(0.5, 0.9)
full_brier <- rep(NA, length(landmark_times))
brier_01 <- rep(NA, length(landmark_times))
brier_02 <- rep(NA, length(landmark_times))
full_auc <- rep(NA, length(landmark_times))
auc_01 <- rep(NA, length(landmark_times))
auc_02 <- rep(NA, length(landmark_times))

for (t in landmark_times){
  y <- ifelse(dat$t > t, 1, 0)
  f_0 <- 1-pnorm(-(x %*% beta_t) -interceptt + log(t), sd = sdy)
  f_01 <- 1-pnorm(-(x[, 2] * beta_t[2]) - interceptt + log(t), sd = sqrt(beta_t[1]^2 + sdy^2))
  f_02 <- 1-pnorm(-(x[, 1] * beta_t[1]) - interceptt + log(t), sd = sqrt(beta_t[2]^2 + sdy^2))
  f_03 <- f_04 <- f_0
  mse_t <- vimp::measure_mse(f_0, y)$point_est
  auc_t <- cvAUC::AUC(f_0, y)
  mse_t1 <- vimp::measure_mse(f_01, y)$point_est
  mse_t2 <- vimp::measure_mse(f_02, y)$point_est
  auc_t1 <- cvAUC::AUC(f_01, y)
  auc_t2 <- cvAUC::AUC(f_02, y)

  full_brier[which(landmark_times == t)] <- -mse_t
  brier_01[which(landmark_times == t)] <- -mse_t1
  brier_02[which(landmark_times == t)] <- -mse_t2
  full_auc[which(landmark_times == t)] <- auc_t
  auc_01[which(landmark_times == t)] <- auc_t1
  auc_02[which(landmark_times == t)] <- auc_t2
}

output1_landmark <- data.frame(vim = rep(c("brier", "AUC"), each = length(landmark_times)),
                               t = c(landmark_times, landmark_times),
                               V_full = c(full_brier, full_auc),
                               V_01 = c(brier_01, auc_01),
                               V_02 = c(brier_02, auc_02),
                               V_03 = rep(0, 2*length(landmark_times)),
                               V_04 = rep(0, 2*length(landmark_times)),
                               V_05 = rep(0, 2*length(landmark_times)),
                               n_mc = rep(n_train, 2*length(landmark_times)))
output1_landmark <- output1_landmark %>% mutate(V_014 = V_01, V_023 = V_02, correlation = FALSE)

taus <- c(0.9)
approx_times <- seq(0, max(taus), by = .01)
S_0 <- matrix(NA, nrow = n_train, ncol = length(approx_times))
S_01 <- matrix(NA, nrow = n_train, ncol = length(approx_times))
S_02 <- matrix(NA, nrow = n_train, ncol = length(approx_times))
for (t in approx_times){
  S_0[,which(approx_times == t)] <- 1-pnorm(-(x %*% beta_t) -interceptt + log(t), sd = sdy)
  S_01[,which(approx_times == t)] <- 1-pnorm(-(x[, 2] * beta_t[2]) - interceptt + log(t), sd = sqrt(beta_t[1]^2 + sdy^2))
  S_02[,which(approx_times == t)] <- 1-pnorm(-(x[, 1] * beta_t[1]) - interceptt + log(t), sd = sqrt(beta_t[2]^2 + sdy^2))
}
# c index stuff
dat_test <- generate_data(n = n_train, scenario = "A")
outcome <- dat_test$t
preds <- -beta_t[1]*dat_test[,1] - beta_t[2]*dat_test[,2]
preds_01 <- -beta_t[2]*dat_test[,2]
preds_02 <- -beta_t[1]*dat_test[,1]
outcome2 <- dat$t
preds2 <- -beta_t[1]*dat[,1] - beta_t[2]*dat[,2]
preds2_01 <- -beta_t[2]*dat[,2]
preds2_02 <- -beta_t[1]*dat[,1]

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

  f_0 <- calc_rmst(S_0, tau, approx_times)
  f_01 <- calc_rmst(S_01, tau, approx_times)
  f_02 <- calc_rmst(S_02, tau, approx_times)
  f_03 <- f_04 <- f_0
}

output1_global <- data.frame(vim = "cindex",
                             t = taus,
                             V_full = full_c,
                             V_01 = c_01,
                             V_02 = c_02,
                             V_03 = 0,
                             V_04 = 0,
                             V_05 = 0,
                             n_mc = n_train)
output1_global <- output1_global %>% mutate(V_014 = V_01, V_023 = V_02, correlation = FALSE)
output1 <- bind_rows(output1_landmark, output1_global)

###############
### correlation
###############
dat <- generate_data(n = n_train, scenario = "D")
p <- 5
x <- as.matrix(dat[,1:p])
beta_t <- matrix(c(0.5, -0.3, rep(0, p- 2)))
rho14 <- 0.7
rho23 <- -0.3
interceptt <- 0

full_brier <- rep(NA, length(landmark_times))
brier_01 <- rep(NA, length(landmark_times))
brier_02 <- rep(NA, length(landmark_times))
full_auc <- rep(NA, length(landmark_times))
auc_01 <- rep(NA, length(landmark_times))
auc_02 <- rep(NA, length(landmark_times))
brier_014 <- rep(NA, length(landmark_times))
brier_023 <- rep(NA, length(landmark_times))
full_auc <- rep(NA, length(landmark_times))
auc_014 <- rep(NA, length(landmark_times))
auc_023 <- rep(NA, length(landmark_times))

for (t in landmark_times){
  y <- ifelse(dat$t > t, 1, 0)
  f_0 <- 1-pnorm(-(x %*% beta_t) -interceptt + log(t), sd = sdy)
  f_01 <- 1-pnorm(-(x[, 2] * beta_t[2]) - interceptt + log(t),
                  mean = beta_t[1]*rho14*x[,4],
                  sd = sqrt(beta_t[1]^2*(1-rho14^2) + sdy^2))
  f_02 <- 1-pnorm(-(x[, 1] * beta_t[1]) - interceptt + log(t),
                  mean = beta_t[2]*rho23*x[,3],
                  sd = sqrt(beta_t[2]^2*(1-rho23^2) + sdy^2))
  f_014 <- 1-pnorm(-(x[, 2] * beta_t[2]) - interceptt + log(t), sd = sqrt(beta_t[1]^2 + sdy^2))
  f_023 <- 1-pnorm(-(x[, 1] * beta_t[1]) - interceptt + log(t), sd = sqrt(beta_t[2]^2 + sdy^2))
  f_03 <- f_04 <- f_0
  mse_t <- vimp::measure_mse(f_0, y)$point_est
  auc_t <- cvAUC::AUC(f_0, y)
  mse_t1 <- vimp::measure_mse(f_01, y)$point_est
  mse_t2 <- vimp::measure_mse(f_02, y)$point_est
  auc_t1 <- cvAUC::AUC(f_01, y)
  auc_t2 <- cvAUC::AUC(f_02, y)

  mse_t14 <- vimp::measure_mse(f_014, y)$point_est
  mse_t23 <- vimp::measure_mse(f_023, y)$point_est
  auc_t14 <- cvAUC::AUC(f_014, y)
  auc_t23 <- cvAUC::AUC(f_023, y)

  full_brier[which(landmark_times == t)] <- -mse_t
  brier_01[which(landmark_times == t)] <- -mse_t1
  brier_02[which(landmark_times == t)] <- -mse_t2
  full_auc[which(landmark_times == t)] <- auc_t
  auc_01[which(landmark_times == t)] <- auc_t1
  auc_02[which(landmark_times == t)] <- auc_t2

  brier_014[which(landmark_times == t)] <- -mse_t14
  brier_023[which(landmark_times == t)] <- -mse_t23
  auc_014[which(landmark_times == t)] <- auc_t14
  auc_023[which(landmark_times == t)] <- auc_t23
}

output2_landmark <- data.frame(vim = rep(c("brier", "AUC"), each = length(landmark_times)),
                               t = c(landmark_times, landmark_times),
                               V_full = c(full_brier, full_auc),
                               V_01 = c(brier_01, auc_01),
                               V_02 = c(brier_02, auc_02),
                               V_014 = c(brier_014, auc_014),
                               V_023 = c(brier_023, auc_023),
                               V_03 = rep(0, 2*length(landmark_times)),
                               V_04 = rep(0, 2*length(landmark_times)),
                               V_05 = rep(0, 2*length(landmark_times)),
                               n_mc = rep(n_train, 2*length(landmark_times)))
output2_landmark <- output2_landmark %>% mutate(correlation = TRUE)

dat_test <- generate_data(n = n_train, scenario = "D")
outcome <- dat_test$t
preds <- -beta_t[1]*dat_test[,1] -beta_t[2]*dat_test[,2]
preds_01 <- -beta_t[2]*dat_test[,2] - beta_t[1]*rho14*dat_test[,4]
preds_02 <- -beta_t[1]*dat_test[,1] -beta_t[2]*rho23*dat_test[,3]
preds_014 <- -beta_t[2]*dat_test[,2]
preds_023 <- -beta_t[1]*dat_test[,1]
outcome2 <- dat$t
preds2 <- -beta_t[1]*dat[,1] -beta_t[2]*dat[,2]
preds2_01 <- -beta_t[2]*dat[,2] - beta_t[1]*rho14*dat[,4]
preds2_02 <- -beta_t[1]*dat[,1] -beta_t[2]*rho23*dat[,3]
preds2_014 <- -beta_t[2]*dat[,2]
preds2_023 <- -beta_t[1]*dat[,1]

full_c <- rep(NA, length(taus))
c_01 <- rep(NA, length(taus))
c_02 <- rep(NA, length(taus))
c_014 <- rep(NA, length(taus))
c_023 <- rep(NA, length(taus))

for (tau in taus){
  denom <- mean((outcome < outcome2 & outcome <= tau) | (outcome2 < outcome & outcome2 <= tau))

  full_c[which(taus == tau)] <- (mean(ifelse(preds > preds2, 1, 0)*(outcome < outcome2 & outcome <= tau)) +
                                   mean(ifelse(preds2 > preds, 1, 0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_01[which(taus == tau)]  <- (mean(ifelse(preds_01 > preds2_01, 1, 0)*(outcome < outcome2 & outcome <= tau)) +
                                  mean(ifelse(preds2_01 > preds_01, 1, 0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_02[which(taus == tau)]  <- (mean(ifelse(preds_02 > preds2_02,1,0)*(outcome < outcome2 & outcome <= tau)) +
                                  mean(ifelse(preds2_02 > preds_02,1,0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_014[which(taus == tau)]  <- (mean(ifelse(preds_014 > preds2_014, 1, 0)*(outcome < outcome2 & outcome <= tau)) +
                                   mean(ifelse(preds2_014 > preds_014, 1, 0)*(outcome2 < outcome & outcome2 <= tau))) / denom
  c_023[which(taus == tau)]  <- (mean(ifelse(preds_023 > preds2_023,1,0)*(outcome < outcome2 & outcome <= tau)) +
                                   mean(ifelse(preds2_023 > preds_023,1,0)*(outcome2 < outcome & outcome2 <= tau))) / denom
}

output2_global <- data.frame(vim = "cindex",
                             t = taus,
                             V_full = full_c,
                             V_01 = c_01,
                             V_02 = c_02,
                             V_014 = c_014,
                             V_023 = c_023,
                             V_03 = 0,
                             V_04 = 0,
                             V_05 = 0,
                             n_mc = n_train)
output2_global <- output2_global %>% mutate(correlation = TRUE)
output2 <- bind_rows(output2_landmark, output2_global)

output <- bind_rows(output1, output2)
saveRDS(output, "truth.rds")
