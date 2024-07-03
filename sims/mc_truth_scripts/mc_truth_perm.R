library(tidyverse)
# source("/home/cwolock/surv_vim_supplementary/sims/generate_data.R")
source("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/sims/generate_data.R")
n_train <- 1e7
set.seed(1234)
sdy <- 1

###############
### correlation
###############
dat <- generate_data(n = n_train, scenario = "1B")
p <- 5
x <- as.matrix(dat[,1:p])
beta_t <- matrix(c(0.5, -0.3, rep(0, p- 2)))
beta_int1 <- 0.2
beta_int2 <- -0.1
rho15 <- 0.9
rho23 <- 0
interceptt <- 0

full_auc <- rep(NA, length(landmark_times))
auc_01 <- rep(NA, length(landmark_times))
auc_02 <- rep(NA, length(landmark_times))
auc_03 <- rep(NA, length(landmark_times))
auc_04 <- rep(NA, length(landmark_times))
auc_05 <- rep(NA, length(landmark_times))

for (t in landmark_times){
  y <- ifelse(dat$t > t, 1, 0)
  f_0 <- 1-pnorm(-(x %*% beta_t + x[,1]*x[,2]*beta_int1 + x[,3]*x[,4]*beta_int2) -interceptt + log(t),
                 mean = 0,
                 sd = sdy)
  f_01 <- 1-pnorm(-(x[,2]*beta_t[2] + beta_int2*x[,3]*x[,4]) - interceptt + log(t),
                  mean = (beta_t[1] + beta_int1*x[,2])*rho15*x[,5],
                  sd = sqrt((beta_t[1] + beta_int1*x[,2])^2*(1-rho15^2) + sdy^2))
  f_02 <- 1-pnorm(-(x[,1]*beta_t[1] + beta_int2*x[,3]*x[,4]) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt((beta_t[2] + beta_int1*x[,1])^2 + sdy^2))
  f_03 <- 1-pnorm(-(x[,1]*beta_t[1] + x[,2]*beta_t[2] + beta_int1*x[,1]*x[,2]) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt((beta_int2*x[,4])^2+ sdy^2))
  f_04 <- 1-pnorm(-(x[,1]*beta_t[1] + x[,2]*beta_t[2] + beta_int1*x[,1]*x[,2]) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt((beta_int2*x[,3])^2 + sdy^2))
  f_05 <- 1-pnorm(-(x %*% beta_t + x[,1]*x[,2]*beta_int1 + x[,3]*x[,4]*beta_int2) -interceptt + log(t),
                  mean = 0,
                  sd = sdy)
  auc_t <- cvAUC::AUC(f_0, y)
  auc_t1 <- cvAUC::AUC(f_01, y)
  auc_t2 <- cvAUC::AUC(f_02, y)
  auc_t3 <- cvAUC::AUC(f_03, y)
  auc_t4 <- cvAUC::AUC(f_04, y)
  auc_t5 <- cvAUC::AUC(f_05, y)

  full_auc[which(landmark_times == t)] <- auc_t
  auc_01[which(landmark_times == t)] <- auc_t1
  auc_02[which(landmark_times == t)] <- auc_t2
  auc_03[which(landmark_times == t)] <- auc_t3
  auc_04[which(landmark_times == t)] <- auc_t4
  auc_05[which(landmark_times == t)] <- auc_t5
}

output2_landmark <- data.frame(vim = rep("AUC", length(landmark_times)),
                               tau = landmark_times,
                               V_full = full_auc,
                               V_01 = auc_01,
                               V_02 = auc_02,
                               V_03 = auc_03,
                               V_04 = auc_04,
                               V_05 = auc_05,
                               n_mc = rep(n_train, length(landmark_times)))
output2_landmark <- output2_landmark %>% mutate(correlation = TRUE)


saveRDS(output2_landmark, "/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/scratch/truth_permutation.rds")
