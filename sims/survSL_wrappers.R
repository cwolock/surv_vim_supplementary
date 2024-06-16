# log normal AFT regression for inclusion in survSuperLearner
survSL.AFTreg <- function(time, event, X, newX, new.times, obsWeights, ...) {
  
  if(any(time == 0 & event == 1)) {
    timepos <- as.numeric(time > 0 & event == 1)
    fit.pos <- stats::glm(timepos ~ ., data=cbind(timepos, X)[event == 1,],
                          family='binomial', weights = obsWeights[event == 1])
    pos.pred <- predict(fit.pos, newdata = newX, type = 'response')
  } else {
    fit.pos <- 1
    pos.pred <- rep(1, nrow(newX))
  }
  
  fit.expreg <- survival::survreg(survival::Surv(time[time > 0], event[time > 0]) ~ .,
                                  data = X[time > 0,],
                                  weights = obsWeights[time > 0], dist = "lognormal")
  pred <- predict(fit.expreg, newdata = newX, type = 'quantile', p = seq(0, .999, by=.001))
  pred <- try(t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001),
                                   xout = new.times, method = 'linear', rule = 2)$y)
  })), silent = TRUE)
  if(inherits(pred, "try-error")) stop("Survival regression failed to produce predictions.")
  
  fit <- list(reg.object = fit.expreg, pos.object = fit.pos)
  class(fit) <- c("survSL.AFTreg")
  out <- list(pred = pred, fit = fit)
  return(out)
}

predict.survSL.AFTreg <- function(object, newX, new.times, ...) {
  
  if(inherits(object$pos.object, "glm")) {
    pos.pred <- predict(object$pos.object, newdata = newX, type = 'response')
  } else {
    pos.pred <- rep(1, nrow(newX))
  }
  
  pred <- predict(object$reg.object, newdata = newX, type = 'quantile',
                  p = seq(0, .999, by=.001))
  pred <- t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001),
                                   xout = new.times, method = 'linear', rule = 2)$y)
  }))
  
  return(pred)
}

# log normal AFT regression for inclusion in survSuperLearner
survSL.AFTreg.int <- function(time, event, X, newX, new.times, obsWeights, ...) {
  
  if(any(time == 0 & event == 1)) {
    timepos <- as.numeric(time > 0 & event == 1)
    fit.pos <- stats::glm(timepos ~ ., data=cbind(timepos, X)[event == 1,],
                          family='binomial', weights = obsWeights[event == 1])
    pos.pred <- predict(fit.pos, newdata = newX, type = 'response')
  } else {
    fit.pos <- 1
    pos.pred <- rep(1, nrow(newX))
  }
  
  fit.expreg <- survival::survreg(survival::Surv(time[time > 0], event[time > 0]) ~ .^2,
                                  data = X[time > 0,],
                                  weights = obsWeights[time > 0], dist = "lognormal")
  pred <- predict(fit.expreg, newdata = newX, type = 'quantile', p = seq(0, .999, by=.001))
  pred <- try(t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001),
                                   xout = new.times, method = 'linear', rule = 2)$y)
  })), silent = TRUE)
  if(inherits(pred, "try-error")) stop("Survival regression failed to produce predictions.")
  
  fit <- list(reg.object = fit.expreg, pos.object = fit.pos)
  class(fit) <- c("survSL.AFTreg")
  out <- list(pred = pred, fit = fit)
  return(out)
}

# log logistic AFT regression for inclusion in survSuperLearner
survSL.loglogreg.int <- function(time, event, X, newX, new.times, obsWeights, ...) {
  
  if(any(time == 0 & event == 1)) {
    timepos <- as.numeric(time > 0 & event == 1)
    fit.pos <- stats::glm(timepos ~ ., data=cbind(timepos, X)[event == 1,],
                          family='binomial', weights = obsWeights[event == 1])
    pos.pred <- predict(fit.pos, newdata = newX, type = 'response')
  } else {
    fit.pos <- 1
    pos.pred <- rep(1, nrow(newX))
  }
  
  fit.expreg <- survival::survreg(survival::Surv(time[time > 0], event[time > 0]) ~ .^2,
                                  data = X[time > 0,],
                                  weights = obsWeights[time > 0], dist = "loglogistic")
  pred <- predict(fit.expreg, newdata = newX, type = 'quantile', p = seq(0, .999, by=.001))
  pred <- try(t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001),
                                   xout = new.times, method = 'linear', rule = 2)$y)
  })), silent = TRUE)
  if(inherits(pred, "try-error")) stop("Survival regression failed to produce predictions.")
  
  fit <- list(reg.object = fit.expreg, pos.object = fit.pos)
  class(fit) <- c("survSL.loglogreg")
  out <- list(pred = pred, fit = fit)
  return(out)
}

# exponential AFT regression for inclusion in survSuperLearner
survSL.expreg.int <- function(time, event, X, newX, new.times, obsWeights, ...) {
  
  if(any(time == 0 & event == 1)) {
    timepos <- as.numeric(time > 0 & event == 1)
    fit.pos <- stats::glm(timepos ~ ., data=cbind(timepos, X)[event == 1,],
                          family='binomial', weights = obsWeights[event == 1])
    pos.pred <- predict(fit.pos, newdata = newX, type = 'response')
  } else {
    fit.pos <- 1
    pos.pred <- rep(1, nrow(newX))
  }
  
  fit.expreg <- survival::survreg(survival::Surv(time[time > 0], event[time > 0]) ~ .^2,
                                  data = X[time > 0,],
                                  weights = obsWeights[time > 0], dist = "exponential")
  pred <- predict(fit.expreg, newdata = newX, type = 'quantile', p = seq(0, .999, by=.001))
  pred <- try(t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001),
                                   xout = new.times, method = 'linear', rule = 2)$y)
  })), silent = TRUE)
  if(inherits(pred, "try-error")) stop("Survival regression failed to produce predictions.")
  
  fit <- list(reg.object = fit.expreg, pos.object = fit.pos)
  class(fit) <- c("survSL.expreg")
  out <- list(pred = pred, fit = fit)
  return(out)
}

# Weibull AFT regression for inclusion in survSuperLearner
survSL.weibreg.int <- function(time, event, X, newX, new.times, obsWeights, ...) {
  
  if(any(time == 0 & event == 1)) {
    timepos <- as.numeric(time > 0 & event == 1)
    fit.pos <- stats::glm(timepos ~ ., data=cbind(timepos, X)[event == 1,],
                          family='binomial', weights = obsWeights[event == 1])
    pos.pred <- predict(fit.pos, newdata = newX, type = 'response')
  } else {
    fit.pos <- 1
    pos.pred <- rep(1, nrow(newX))
  }
  
  fit.expreg <- survival::survreg(survival::Surv(time[time > 0], event[time > 0]) ~ .^2,
                                  data = X[time > 0,],
                                  weights = obsWeights[time > 0], dist = "weibull")
  pred <- predict(fit.expreg, newdata = newX, type = 'quantile', p = seq(0, .999, by=.001))
  pred <- try(t(sapply(1:nrow(pred), function(j) {
    pos.pred[j] * (1-stats::approx(pred[j,], seq(0, .999, by=.001),
                                   xout = new.times, method = 'linear', rule = 2)$y)
  })), silent = TRUE)
  if(inherits(pred, "try-error")) stop("Survival regression failed to produce predictions.")
  
  fit <- list(reg.object = fit.expreg, pos.object = fit.pos)
  class(fit) <- c("survSL.weibreg")
  out <- list(pred = pred, fit = fit)
  return(out)
}