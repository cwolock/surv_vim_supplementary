do_one <- function(n_train, correlation){
  scenario <- ifelse(correlation, "2B", "4B")
  sdy <- 1
  # training data
  dat <- generate_data(n = n_train, scenario = scenario, sdy = sdy)
  x <- as.matrix(dat[,1:6])
  dat2 <- generate_data(n = n_train, scenario = scenario, sdy = sdy)
  x2 <- as.matrix(dat2[,1:6])
  p <- 6
  beta_t <- matrix(c(0.5, -0.3, rep(0, p- 2)))
  beta_c <- matrix(c(-0.2, 0.2, rep(0, p- 2)))
  interceptt <- 0
  interceptc <- 0
  beta_int <- c(0.1, -0.1, 0.1)

  rho16 <- 0.7
  rho23 <- -0.3

  # benchmarks
  approx_times <- seq(0, 0.9, by = 0.01)
  landmark_times <- c(0.5,0.9)
  c_landmark_times <- c(0.9)

  brier_01 <- rep(NA, length(landmark_times))
  brier_02 <- rep(NA, length(landmark_times))
  auc_01 <- rep(NA, length(landmark_times))
  auc_02 <- rep(NA, length(landmark_times))
  brier_01_split <- rep(NA, length(landmark_times))
  brier_06_split <- rep(NA, length(landmark_times))
  brier_016_split <- rep(NA, length(landmark_times))
  auc_01_split <- rep(NA, length(landmark_times))
  auc_06_split <- rep(NA, length(landmark_times))
  auc_016_split <- rep(NA, length(landmark_times))

  S0 <- matrix(NA, nrow = nrow(dat), ncol = length(approx_times))
  G0 <- matrix(NA, nrow = nrow(dat), ncol = length(approx_times))
  S02 <- matrix(NA, nrow = nrow(dat), ncol = length(approx_times))
  G02 <- matrix(NA, nrow = nrow(dat), ncol = length(approx_times))
  for (t in approx_times){
    S <- 1-pnorm(-(x %*% beta_t + x[,1]*x[,2]*beta_int[1] + x[,3]*x[,4]*beta_int[2] + x[,1]*x[,5]*beta_int[3]) -interceptt + log(t), sd = sdy)
    G <- 1-pnorm(-(as.matrix(x) %*% beta_c) -interceptc + log(t), sd = sdy)
    S0[,which(approx_times == t)] <- S
    G0[,which(approx_times == t)] <- G
    S2 <- 1-pnorm(-(as.matrix(x2) %*% beta_t + x2[,1]*x2[,2]*beta_int[1] + x2[,3]*x2[,4]*beta_int[2] + x2[,1]*x2[,5]*beta_int[3]) -interceptt + log(t), sd = sdy)
    G2 <- 1-pnorm(-(as.matrix(x2) %*% beta_c) -interceptc + log(t), sd = sdy)
    S02[,which(approx_times == t)] <- S2
    G02[,which(approx_times == t)] <- G2
  }

  for (t in landmark_times){
    if (scenario == "4"){
      f <- 1-pnorm(-(x %*% beta_t + x[,1]*x[,2]*beta_int[1] + x[,3]*x[,4]*beta_int[2] + x[,1]*x[,5]*beta_int[3]) -interceptt + log(t),
                     mean = 0,
                     sd = sdy)
      f2 <- 1-pnorm(-(x2 %*% beta_t + x2[,1]*x2[,2]*beta_int[1] + x2[,3]*x2[,4]*beta_int[2] + x2[,1]*x2[,5]*beta_int[3]) -interceptt + log(t),
                   mean = 0,
                   sd = sdy)
      f_01 <- 1-pnorm(-(x[,2]*beta_t[2] + beta_int[2]*x[,3]*x[,4]) - interceptt + log(t),
                      mean = (beta_t[1] + beta_int[1]*x[,2] + beta_int[3]*x[,5])*rho16*x[,6],
                      sd = sqrt((beta_t[1] + beta_int[1]*x[,2] + beta_int[3]*x[,5])^2*(1-rho16^2) + sdy^2))
      f_01_2 <- 1-pnorm(-(x2[,2]*beta_t[2] + beta_int[2]*x2[,3]*x2[,4]) - interceptt + log(t),
                      mean = (beta_t[1] + beta_int[1]*x2[,2] + beta_int[3]*x2[,5])*rho16*x2[,6],
                      sd = sqrt((beta_t[1] + beta_int[1]*x2[,2] + beta_int[3]*x2[,5])^2*(1-rho16^2) + sdy^2))
      f_02 <- 1-pnorm(-(x[,1]*beta_t[1] + beta_int[2]*x[,3]*x[,4] + beta_int[3]*x[,1]*x[,5]) - interceptt + log(t),
                      mean = (beta_t[2] + beta_int[1]*x[,1])*rho23*x[,3],
                      sd = sqrt((beta_t[2] + beta_int[1]*x[,1])^2*(1-rho23^2) + sdy^2))
      f_016 <- 1-pnorm(-(x[,2]*beta_t[2] + beta_int[2]*x[,3]*x[,4]) - interceptt + log(t),
                       mean = 0,
                       sd = sqrt((beta_t[1] + beta_int[1]*x[,2] + beta_int[3]*x[,5])^2 + sdy^2))
      f_016_2 <- 1-pnorm(-(x2[,2]*beta_t[2] + beta_int[2]*x2[,3]*x2[,4]) - interceptt + log(t),
                       mean = 0,
                       sd = sqrt((beta_t[1] + beta_int[1]*x2[,2] + beta_int[3]*x2[,5])^2 + sdy^2))
    } else{
      f <- 1-pnorm(-(x %*% beta_t + x[,1]*x[,2]*beta_int[1] + x[,3]*x[,4]*beta_int[2] + x[,1]*x[,5]*beta_int[3]) -interceptt + log(t),
                   mean = 0,
                   sd = sdy)
      f2 <- 1-pnorm(-(x2 %*% beta_t + x2[,1]*x2[,2]*beta_int[1] + x2[,3]*x2[,4]*beta_int[2] + x2[,1]*x2[,5]*beta_int[3]) -interceptt + log(t),
                    mean = 0,
                    sd = sdy)
      f_01 <- 1-pnorm(-(x[,2]*beta_t[2] + beta_int[2]*x[,3]*x[,4]) - interceptt + log(t),
                      mean = 0,
                      sd = sqrt((beta_t[1] + beta_int[1]*x[,2] + beta_int[3]*x[,5])^2 + sdy^2))
      f_01_2 <- 1-pnorm(-(x2[,2]*beta_t[2] + beta_int[2]*x2[,3]*x2[,4]) - interceptt + log(t),
                        mean = 0,
                        sd = sqrt((beta_t[1] + beta_int[1]*x2[,2] + beta_int[3]*x2[,5])^2 + sdy^2))
      f_02 <- 1-pnorm(-(x[,1]*beta_t[1] + beta_int[2]*x[,3]*x[,4] + beta_int[3]*x[,1]*x[,5]) - interceptt + log(t),
                      mean = 0,
                      sd = sqrt((beta_t[2] + beta_int[1]*x[,1])^2* + sdy^2))

      f_016 <- f_01
      f_016_2 <- f_01_2
    }

    # brier
    V_0 <- survML:::estimate_brier(time = dat$y,
                                   event = dat$delta,
                                   approx_times = approx_times,
                                   tau = t,
                                   preds = f,
                                   S_hat = S0,
                                   G_hat = G0)
    V_01 <- survML:::estimate_brier(time = dat$y,
                                    event = dat$delta,
                                    approx_times = approx_times,
                                    tau = t,
                                    preds = f_01,
                                    S_hat = S0,
                                    G_hat = G0)
    V_02 <- survML:::estimate_brier(time = dat$y,
                                    event = dat$delta,
                                    approx_times = approx_times,
                                    tau = t,
                                    preds = f_02,
                                    S_hat = S0,
                                    G_hat = G0)
    V_016 <- survML:::estimate_brier(time = dat$y,
                                    event = dat$delta,
                                    approx_times = approx_times,
                                    tau = t,
                                    preds = f_016,
                                    S_hat = S0,
                                    G_hat = G0)
    V_016_2 <- survML:::estimate_brier(time = dat2$y,
                                     event = dat2$delta,
                                     approx_times = approx_times,
                                     tau = t,
                                     preds = f_016_2,
                                     S_hat = S02,
                                     G_hat = G02)
    V_0_2 <- survML:::estimate_brier(time = dat2$y,
                                   event = dat2$delta,
                                   approx_times = approx_times,
                                   tau = t,
                                   preds = f2,
                                   S_hat = S02,
                                   G_hat = G02)

    var_01 <- mean((V_0$EIF - V_01$EIF)^2)
    var_02 <- mean((V_0$EIF - V_02$EIF)^2)
    brier_01[which(landmark_times == t)] <- var_01
    brier_02[which(landmark_times == t)] <- var_02

    var_01_split <- mean((V_0_2$EIF)^2) + mean((V_01$EIF)^2)
    var_06_split <- mean((V_0$EIF)^2) + mean((V_0_2$EIF)^2)
    var_016_split <- mean((V_016$EIF)^2) + mean((V_016_2$EIF)^2)
    brier_01_split[which(landmark_times == t)] <- var_01_split
    brier_06_split[which(landmark_times == t)] <- var_06_split
    brier_016_split[which(landmark_times == t)] <- var_016_split

    # auc
    V_0 <- survML:::estimate_AUC(time = dat$y,
                                 event = dat$delta,
                                 approx_times = approx_times,
                                 tau = t,
                                 preds = 1-f,
                                 S_hat = S0,
                                 G_hat = G0,
                                 robust = TRUE)
    V_0_2 <- survML:::estimate_AUC(time = dat2$y,
                                 event = dat2$delta,
                                 approx_times = approx_times,
                                 tau = t,
                                 preds = 1-f2,
                                 S_hat = S02,
                                 G_hat = G02,
                                 robust = TRUE)
    V_01 <- survML:::estimate_AUC(time = dat$y,
                                  event = dat$delta,
                                  approx_times = approx_times,
                                  tau = t,
                                  preds = 1-f_01,
                                  S_hat = S0,
                                  G_hat = G0,
                                  robust = TRUE)
    V_02 <- survML:::estimate_AUC(time = dat$y,
                                  event = dat$delta,
                                  approx_times = approx_times,
                                  tau = t,
                                  preds = 1-f_02,
                                  S_hat = S0,
                                  G_hat = G0,
                                  robust = TRUE)
    V_016 <- survML:::estimate_AUC(time = dat$y,
                                  event = dat$delta,
                                  approx_times = approx_times,
                                  tau = t,
                                  preds = 1-f_016,
                                  S_hat = S0,
                                  G_hat = G0,
                                  robust = TRUE)
    V_016_2 <- survML:::estimate_AUC(time = dat2$y,
                                   event = dat2$delta,
                                   approx_times = approx_times,
                                   tau = t,
                                   preds = 1-f_016_2,
                                   S_hat = S02,
                                   G_hat = G02,
                                   robust = TRUE)

    var_01 <- mean((V_0$EIF - V_01$EIF)^2)
    var_02 <- mean((V_0$EIF - V_02$EIF)^2)
    auc_01[which(landmark_times == t)] <- var_01
    auc_02[which(landmark_times == t)] <- var_02
    var_01_split <- mean((V_0_2$EIF)^2) + mean((V_01$EIF)^2)
    var_06_split <- mean((V_0$EIF)^2) + mean((V_0_2$EIF)^2)
    var_016_split <- mean((V_016$EIF)^2) + mean((V_016_2$EIF)^2)
    auc_01_split[which(landmark_times == t)] <- var_01_split
    auc_06_split[which(landmark_times == t)] <- var_06_split
    auc_016_split[which(landmark_times == t)] <- var_016_split
  }

  # cindex
  c_01 <- rep(NA, length(c_landmark_times))
  c_02 <- rep(NA, length(c_landmark_times))
  c_01_split <- rep(NA, length(c_landmark_times))
  c_06_split <- rep(NA, length(c_landmark_times))
  c_016_split <- rep(NA, length(c_landmark_times))

  S0 <- matrix(NA, nrow = nrow(dat), ncol = length(approx_times))
  G0 <- matrix(NA, nrow = nrow(dat), ncol = length(approx_times))
  S02 <- matrix(NA, nrow = nrow(dat), ncol = length(approx_times))
  G02 <- matrix(NA, nrow = nrow(dat), ncol = length(approx_times))
  for (t in approx_times){
    S <- 1-pnorm(-(x %*% beta_t + x[,1]*x[,2]*beta_int[1] + x[,3]*x[,4]*beta_int[2] + x[,1]*x[,5]*beta_int[3]) -interceptt + log(t), sd = sdy)
    G <- 1-pnorm(-(as.matrix(x) %*% beta_c) -interceptc + log(t), sd = sdy)
    S0[,which(approx_times == t)] <- S
    G0[,which(approx_times == t)] <- G
    S2 <- 1-pnorm(-(as.matrix(x2) %*% beta_t + x2[,1]*x2[,2]*beta_int[1] + x2[,3]*x2[,4]*beta_int[2] + x2[,1]*x2[,5]*beta_int[3]) -interceptt + log(t), sd = sdy)
    G2 <- 1-pnorm(-(as.matrix(x2) %*% beta_c) -interceptc + log(t), sd = sdy)
    S02[,which(approx_times == t)] <- S2
    G02[,which(approx_times == t)] <- G2
  }

  for (t in c_landmark_times){
    if (scenario == "4A"){
      f <- -(x %*% beta_t + x[,1]*x[,2]*beta_int[1] + x[,3]*x[,4]*beta_int[2] + x[,1]*x[,5]*beta_int[3])
      f2 <- -(x2 %*% beta_t + x2[,1]*x2[,2]*beta_int[1] + x2[,3]*x2[,4]*beta_int[2] + x2[,1]*x2[,5]*beta_int[3])
      f_01 <- -(x[, 2] * beta_t[2] + beta_int[2]*x[,3]*x[,4] + (beta_t[1] + beta_int[1]*x[,2] + beta_int[3]*x[,5])*rho16*x[,6])
      f_02 <- -(x[, 1] * beta_t[1] + beta_int[2]*x2[,3]*x2[,4] + (beta_t[2] + beta_int[1]*x[,1])*rho23*x[,3])
      f_016 <- -(x[, 2] * beta_t[2] + beta_int[2]*x[,3]*x[,4])
      f_016_2 <- -(x2[, 2] * beta_t[2]+ beta_int[2]*x2[,3]*x2[,4])
    } else{
      f <- -(x %*% beta_t + x[,1]*x[,2]*beta_int[1] + x[,3]*x[,4]*beta_int[2] + x[,1]*x[,5]*beta_int[3])
      f2 <- -(x2 %*% beta_t + x2[,1]*x2[,2]*beta_int[1] + x2[,3]*x2[,4]*beta_int[2] + x2[,1]*x2[,5]*beta_int[3])
      f_01 <- -(x[, 2] * beta_t[2] + beta_int[2]*x[,3]*x[,4])
      f_01_2 <- -(x2[, 2] * beta_t[2] + beta_int[2]*x2[,3]*x2[,4])
      f_02 <- -(x[, 1] * beta_t[1]+ beta_int[2]*x[,3]*x[,4])
      f_016 <- f_01
      f_016_2 <- f_01_2
    }

    V_0 <- survML:::estimate_cindex(time = dat$y,
                                    event = dat$delta,
                                    approx_times = approx_times,
                                    tau = t,
                                    preds = f,
                                    S_hat = S0,
                                    G_hat = G0)
    V_0_2 <- survML:::estimate_cindex(time = dat2$y,
                                      event = dat2$delta,
                                      approx_times = approx_times,
                                      tau = t,
                                      preds = f2,
                                      S_hat = S02,
                                      G_hat = G02)
    V_01 <- survML:::estimate_cindex(time = dat$y,
                                     event = dat$delta,
                                     approx_times = approx_times,
                                     tau = t,
                                     preds = f_01,
                                     S_hat = S0,
                                     G_hat = G0)
    V_02 <- survML:::estimate_cindex(time = dat$y,
                                     event = dat$delta,
                                     approx_times = approx_times,
                                     tau = t,
                                     preds = f_02,
                                     S_hat = S0,
                                     G_hat = G0)
    V_016 <- survML:::estimate_cindex(time = dat$y,
                                      event = dat$delta,
                                      approx_times = approx_times,
                                      tau = t,
                                      preds = f_016,
                                      S_hat = S0,
                                      G_hat = G0)
    V_016_2 <- survML:::estimate_cindex(time = dat2$y,
                                        event = dat2$delta,
                                        approx_times = approx_times,
                                        tau = t,
                                        preds = f_016_2,
                                        S_hat = S02,
                                        G_hat = G02)

    var_01 <- mean((V_0$EIF - V_01$EIF)^2)
    var_02 <- mean((V_0$EIF - V_02$EIF)^2)
    c_01[which(c_landmark_times == t)] <- var_01
    c_02[which(c_landmark_times == t)] <- var_02

    var_01_split <- mean((V_0_2$EIF)^2) + mean((V_01$EIF)^2)
    var_06_split <- mean((V_0$EIF)^2) + mean((V_0_2$EIF)^2)
    var_016_split <- mean((V_016$EIF)^2) + mean((V_016_2$EIF)^2)
    c_01_split[which(c_landmark_times == t)] <- var_01_split
    c_06_split[which(c_landmark_times == t)] <- var_06_split
    c_016_split[which(c_landmark_times == t)] <- var_016_split
  }

  output <- data.frame(vim = rep(c("brier", "AUC"), each = length(landmark_times)),
                       tau = c(landmark_times, landmark_times),
                       vim_1 = c(brier_01, auc_01),
                       vim_2 = c(brier_02, auc_02),
                       vim_1_split = c(brier_01_split, auc_01_split),
                       vim_6_split = c(brier_06_split, auc_06_split),
                       vim_16_split = c(brier_016_split, auc_016_split),
                       correlation = correlation)

  output_2 <- data.frame(vim = "cindex",
                         tau = c_landmark_times,
                         vim_1 = c_01,
                         vim_2 = c_02,
                         vim_1_split = c_01_split,
                         vim_6_split = c_06_split,
                         vim_16_split = c_016_split,
                         correlation = correlation)
  return(bind_rows(output, output_2))
}
