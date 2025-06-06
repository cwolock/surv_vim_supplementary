library(tidyverse)
# source("/home/cwolock/surv_vim_supplementary/sims/generate_data.R")
source("/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/sims/generate_data.R")
n_train <- 2e6
set.seed(1234)
sdy <- 1
##################
### no correlation
##################
tau <- 0.9
p <- 5
dat <- generate_data(n = n_train, scenario = "2", sdy = sdy)
x <- as.matrix(dat[,1:p])
beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
beta_int <- c(0.1, -0.1, 0.1)
interceptt <- 0
grid <- 0.002
approx_times <- seq(0, tau, by = grid)

# full_rmst <- rep(NA, length(landmark_times))
# # rmst_01 <- rep(NA, length(landmark_times))
# rmst_02 <- rep(NA, length(landmark_times))

f_0 <- matrix(NA, nrow = nrow(x), ncol = length(approx_times))
f_01 <- matrix(NA, nrow = nrow(x), ncol = length(approx_times))
f_02 <- matrix(NA, nrow = nrow(x), ncol = length(approx_times))

for (i in 1:length(approx_times)){
  t <- approx_times[i]
  f_0_t <- 1-pnorm(-(x %*% beta_t + x[,1]*x[,2]*beta_int[1] + x[,3]*x[,4]*beta_int[2]+x[,1]*x[,5]*beta_int[3]) -interceptt + log(t),
                 mean = 0,
                 sd = sdy)
  f_01_t <- 1-pnorm(-(x[, 2] * beta_t[2] + x[,3]*x[,4]*beta_int[2] ) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt((beta_t[1] + beta_int[1]*x[,2] + beta_int[3]*x[,5])^2 + sdy^2))
  f_02_t <- 1-pnorm(-(x[, 1] * beta_t[1] + x[,3]*x[,4]*beta_int[2] + x[,1]*x[,5]*beta_int[3]) - interceptt + log(t),
                  mean = 0,
                  sd = sqrt((beta_t[2] + beta_int[1]*x[,1])^2 + sdy^2))
  f_0[,i] <- matrix(f_0_t, nrow = nrow(x), ncol = 1)
  f_01[,i] <- matrix(f_01_t, nrow = nrow(x), ncol = 1)
  f_02[,i] <- matrix(f_02_t, nrow = nrow(x), ncol = 1)
}

f_0 <- rowSums(f_0) * grid
f_01 <- rowSums(f_01) * grid
f_02 <- rowSums(f_02) * grid

restricted_ts <- pmin(dat$t, tau)


mse_t <- vimp::measure_mse(f_0, restricted_ts)$point_est
mse_t1 <- vimp::measure_mse(f_01, restricted_ts)$point_est
mse_t2 <- vimp::measure_mse(f_02, restricted_ts)$point_est

output <- data.frame(vim = "survival_time_MSE",
                     tau = tau,
                     V_full = -mse_t,
                     V_01 = -mse_t1,
                     V_02 = -mse_t2,
                     V_06 = -mse_t,
                     n_mc = n_train)

saveRDS(output, "/Users/cwolock/Dropbox/UW/RESEARCH/paper_supplements/surv_vim_supplementary/scratch/truth_interactionB_rmst.rds")
