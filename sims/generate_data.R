generate_data <- function(n = 500, scenario = "1", sdy = 1, max_fu = 100){

  if (scenario == "1" | scenario == "3_50"){
    p <- 2
    beta_t <- matrix(c(0.5, -0.3))
    beta_c <- matrix(c(-0.2, 0.2))
    interceptc <- 0
    xnames <- paste0("x", c(1,2))
    Sigma <- diag(1, p)
  } else if (scenario == "1A"){
    p <- 4
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int1 <- 0.2
    beta_int2 <- -0.1
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "5A"){
    p <- 4
    beta_t <- matrix(c(0.5, -0.3, -0.1, 0))
    beta_c <- matrix(rep(0, p))
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
    Sigma[1,4] <- Sigma[4,1] <- 0.9
    # Sigma[2,3] <- Sigma[3,2] <- -0.3
  } else if (scenario == "5B"){
    p <- 4
    beta_t <- matrix(c(0.5, -0.3, -0.1, 0))
    beta_c <- matrix(rep(0, p))
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
    Sigma[1,4] <- Sigma[4,1] <- 0.6
    # Sigma[2,3] <- Sigma[3,2] <- -0.3
  } else if (scenario == "5C"){
    p <- 4
    beta_t <- matrix(c(0.5, -0.3, -0.1, 0))
    beta_c <- matrix(rep(0, p))
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
    Sigma[1,4] <- Sigma[4,1] <- 0.3
    # Sigma[2,3] <- Sigma[3,2] <- -0.3
  } else if (scenario == "5D"){
    p <- 4
    beta_t <- matrix(c(0.5, -0.3, -0.1, 0))
    beta_c <- matrix(rep(0, p))
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
    Sigma[1,4] <- Sigma[4,1] <- 0
    # Sigma[2,3] <- Sigma[3,2] <- -0.3
  } else if (scenario == "2"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, (p-2))))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, (p-2))))
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "2A"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, (p-2))))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, (p-2))))
    beta_int1 <- 0.2
    beta_int2 <- -0.1
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "4A"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, (p-2))))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, (p-2))))
    beta_int1 <- 0.2
    beta_int2 <- -0.1
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
    Sigma[1,5] <- Sigma[5,1] <- 0.7
    Sigma[2,3] <- Sigma[3,2] <- -0.3
  } else if (scenario == "4"){
    p <- 5
    beta_t <- matrix(c(0.5, -0.3, rep(0, (p-2))))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, (p-2))))
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
    Sigma[1,4] <- Sigma[4,1] <- 0.7
    Sigma[2,3] <- Sigma[3,2] <- -0.3
  } else if (scenario == "3_30"){
    p <- 2
    beta_t <- matrix(c(0.5, -0.3))
    beta_c <- matrix(c(-0.2, 0.2))
    interceptc <- 0.75
    xnames <- paste0("x", c(1,2))
    Sigma <- diag(1, p)
  } else if (scenario == "3_40"){
    p <- 2
    beta_t <- matrix(c(0.5, -0.3))
    beta_c <- matrix(c(-0.2, 0.2))
    interceptc <- 0.5
    xnames <- paste0("x", c(1,2))
    Sigma <- diag(1, p)
  } else if (scenario == "3_60"){
    p <- 2
    beta_t <- matrix(c(0.5, -0.3))
    beta_c <- matrix(c(-0.2, 0.2))
    interceptc <- -0.5
    xnames <- paste0("x", c(1,2))
    Sigma <- diag(1, p)
  } else if (scenario == "3_70"){
    p <- 2
    beta_t <- matrix(c(0.5, -0.3))
    beta_c <- matrix(c(-0.2, 0.2))
    interceptc <- -0.75
    xnames <- paste0("x", c(1,2))
    Sigma <- diag(1, p)
  } else if (scenario == "3A_30"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int1 <- 0.2
    beta_int2 <- -0.1
    interceptc <- 0.85
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "3A_40"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int1 <- 0.2
    beta_int2 <- -0.1
    interceptc <- 0.5
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "3A_50"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int1 <- 0.2
    beta_int2 <- -0.1
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "3A_60"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int1 <- 0.2
    beta_int2 <- -0.1
    interceptc <- -0.45
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "3A_70"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int1 <- 0.2
    beta_int2 <- -0.1
    interceptc <- -0.85
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  }

  mu_x <- rep(0, p)

  x <- MASS::mvrnorm(n = n, mu = mu_x, Sigma = Sigma)
  eps <- rnorm(n = n, mean = 0, sd = sdy)
  epsc <- rnorm(n = n, mean = 0, sd = sdy)
  if (scenario != "2A" & scenario != "4A" & scenario != "1A" & scenario != "3A_50" & scenario != "3A_30" & scenario != "3A_40" & scenario != "3A_60" & scenario != "3A_70"){
    logt <- x %*% beta_t + eps
  } else{
    logt <- x %*% beta_t + x[,1]*x[,2]*beta_int1 + x[,3]*x[,4]*beta_int2 + eps
  }

  t <- exp(logt)

  logc <- interceptc + x %*% beta_c + epsc
  c <- exp(logc)
  c <- ifelse(c <= max_fu, c, max_fu)

  y <- ifelse(t <= c, t, c)
  delta <- ifelse(t <= c, 1, 0)

  data <- data.frame(x,
                     y = y,
                     delta = delta,
                     t = t,
                     c = c)
  names(data) <- c(xnames, "y", "delta", "t", "c")
  return(data)
}
