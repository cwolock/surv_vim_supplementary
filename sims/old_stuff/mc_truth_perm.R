library(tidyverse)
# source("/home/cwolock/surv_vim_supplementary/sims/generate_data.R")
source("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/generate_data.R")
n_train <- 2e7
set.seed(1234)
sdy <- 1
landmark_times <- c(0.5, 0.9)

###############
### no correlation
###############
dat <- generate_data(n = n_train, scenario = "1C", sdy = sdy)
p <- 4
x <- as.matrix(dat[,1:p])
beta_t <- matrix(c(0.5, -0.3, -0.15, 0))
rho14 <- 0
rho23 <- 0
interceptt <- 0

full_auc <- rep(NA, length(landmark_times))
auc_01 <- rep(NA, length(landmark_times))
auc_02 <- rep(NA, length(landmark_times))
auc_03 <- rep(NA, length(landmark_times))
auc_04 <- rep(NA, length(landmark_times))

for (t in landmark_times){
  y <- ifelse(dat$t > t, 1, 0)
  f_0 <- 1-pnorm(-(x %*% beta_t) -interceptt + log(t),
                 mean = 0,
                 sd = sdy)
  f_01 <- 1-pnorm(-(x[,2]*beta_t[2] + beta_t[3]*x[,3] ) - interceptt + log(t),
                  mean = beta_t[1]*rho14*x[,4],
                  sd = sqrt(beta_t[1]^2*(1-rho14^2) + sdy^2))
  f_02 <- 1-pnorm(-(x[,1]*beta_t[1] + beta_t[3]*x[,3] ) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt(beta_t[2]^2 + sdy^2))
  f_03 <- 1-pnorm(-(x[,1]*beta_t[1] + x[,2]*beta_t[2] ) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt((beta_t[3])^2+ sdy^2))
  f_04 <- 1-pnorm(-(x %*% beta_t) -interceptt + log(t),
                  mean = 0,
                  sd = sdy)
  auc_t <- cvAUC::AUC(f_0, y)
  auc_t1 <- cvAUC::AUC(f_01, y)
  auc_t2 <- cvAUC::AUC(f_02, y)
  auc_t3 <- cvAUC::AUC(f_03, y)
  auc_t4 <- cvAUC::AUC(f_04, y)

  full_auc[which(landmark_times == t)] <- auc_t
  auc_01[which(landmark_times == t)] <- auc_t1
  auc_02[which(landmark_times == t)] <- auc_t2
  auc_03[which(landmark_times == t)] <- auc_t3
  auc_04[which(landmark_times == t)] <- auc_t4
}

output1_landmark <- data.frame(vim = rep("AUC", length(landmark_times)),
                               tau = landmark_times,
                               V_full = full_auc,
                               V_01 = auc_01,
                               V_02 = auc_02,
                               V_03 = auc_03,
                               V_04 = auc_04,
                               n_mc = rep(n_train, length(landmark_times)))
output1_landmark <- output1_landmark %>% mutate(correlation = FALSE)


###############
### correlation
###############
dat <- generate_data(n = n_train, scenario = "1B", sdy = sdy)
p <- 4
x <- as.matrix(dat[,1:p])
beta_t <- matrix(c(0.5, -0.3, -0.15, 0))
rho14 <- 0.9
rho23 <- 0
interceptt <- 0

full_auc <- rep(NA, length(landmark_times))
auc_01 <- rep(NA, length(landmark_times))
auc_02 <- rep(NA, length(landmark_times))
auc_03 <- rep(NA, length(landmark_times))
auc_04 <- rep(NA, length(landmark_times))

for (t in landmark_times){
  y <- ifelse(dat$t > t, 1, 0)
  f_0 <- 1-pnorm(-(x %*% beta_t) -interceptt + log(t),
                 mean = 0,
                 sd = sdy)
  f_01 <- 1-pnorm(-(x[,2]*beta_t[2] + beta_t[3]*x[,3] ) - interceptt + log(t),
                  mean = beta_t[1]*rho14*x[,4],
                  sd = sqrt(beta_t[1]^2*(1-rho14^2) + sdy^2))
  f_02 <- 1-pnorm(-(x[,1]*beta_t[1] + beta_t[3]*x[,3] ) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt(beta_t[2]^2 + sdy^2))
  f_03 <- 1-pnorm(-(x[,1]*beta_t[1] + x[,2]*beta_t[2] ) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt((beta_t[3])^2+ sdy^2))
  f_04 <- 1-pnorm(-(x %*% beta_t) -interceptt + log(t),
                  mean = 0,
                  sd = sdy)
  auc_t <- cvAUC::AUC(f_0, y)
  auc_t1 <- cvAUC::AUC(f_01, y)
  auc_t2 <- cvAUC::AUC(f_02, y)
  auc_t3 <- cvAUC::AUC(f_03, y)
  auc_t4 <- cvAUC::AUC(f_04, y)

  full_auc[which(landmark_times == t)] <- auc_t
  auc_01[which(landmark_times == t)] <- auc_t1
  auc_02[which(landmark_times == t)] <- auc_t2
  auc_03[which(landmark_times == t)] <- auc_t3
  auc_04[which(landmark_times == t)] <- auc_t4
}

output2_landmark <- data.frame(vim = rep("AUC", length(landmark_times)),
                               tau = landmark_times,
                               V_full = full_auc,
                               V_01 = auc_01,
                               V_02 = auc_02,
                               V_03 = auc_03,
                               V_04 = auc_04,
                               n_mc = rep(n_train, length(landmark_times)))
output2_landmark <- output2_landmark %>% mutate(correlation = TRUE)

output <- rbind(output1_landmark, output2_landmark)


saveRDS(output, "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/truth_permutation.rds")
