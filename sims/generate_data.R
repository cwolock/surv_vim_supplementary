generate_data <- function(n = 500, scenario = "1", sdy = 1, max_fu = 100){
  
  if (scenario == "1"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, (p-2))))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, (p-2))))
    beta_int <- c(0.1, -0.1, 0.1)
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
    Sigma[1,6] <- Sigma[6,1] <- 0.7
    Sigma[2,3] <- Sigma[3,2] <- -0.3
  } else if (scenario == "2"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, (p-2))))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, (p-2))))
    beta_int <- c(0.1, -0.1, 0.1)
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "3"){
    p <- 5
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int <- c(0.1, -0.1, 0.1)
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "4_30"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int <- c(0.1, -0.1, 0.1)
    interceptc <- 0.85
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "4_40"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int <- c(0.1, -0.1, 0.1)
    interceptc <- 0.5
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "4_50"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int <- c(0.1, -0.1, 0.1)
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "4_60"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int <- c(0.1, -0.1, 0.1)
    interceptc <- -0.45
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "4_70"){
    p <- 25
    beta_t <- matrix(c(0.5, -0.3, rep(0, p-2)))
    beta_c <- matrix(c(-0.2, 0.2, rep(0, p-2)))
    beta_int <- c(0.1, -0.1, 0.1)
    interceptc <- -0.85
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
  } else if (scenario == "permA"){
    p <- 4
    beta_t <- matrix(c(0.5, -0.3, -0.15, 0))
    beta_c <- matrix(rep(0, p))
    beta_int <- c(0,0,0)
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
    Sigma[1,4] <- Sigma[4,1] <- 0
  } else if (scenario == "permB"){
    p <- 4
    beta_t <- matrix(c(0.5, -0.3, -0.15, 0))
    beta_c <- matrix(rep(0, p))
    beta_int <- c(0,0,0)
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
    Sigma[1,4] <- Sigma[4,1] <- 0.3
  } else if (scenario == "permC"){
    p <- 4
    beta_t <- matrix(c(0.5, -0.3, -0.15, 0))
    beta_c <- matrix(rep(0, p))
    beta_int <- c(0,0,0)
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
    Sigma[1,4] <- Sigma[4,1] <- 0.6
  } else if (scenario == "permD"){
    p <- 4
    beta_t <- matrix(c(0.5, -0.3, -0.15, 0))
    beta_c <- matrix(rep(0, p))
    beta_int <- c(0,0,0)
    interceptc <- 0
    xnames <- paste0("x", 1:p)
    Sigma <- diag(1, p)
    Sigma[1,4] <- Sigma[4,1] <- 0.9
  }
  
  mu_x <- rep(0, p)
  
  x <- MASS::mvrnorm(n = n, mu = mu_x, Sigma = Sigma)
  eps <- rnorm(n = n, mean = 0, sd = sdy)
  epsc <- rnorm(n = n, mean = 0, sd = sdy)
  logt <- x %*% beta_t + x[,1]*x[,2]*beta_int[1] + x[,3]*x[,4]*beta_int[2] + x[,1]*x[,5]*beta_int[3] + eps
  
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
