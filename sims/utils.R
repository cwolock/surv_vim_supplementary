generate_full_predictions <- function(time,
                                      event,
                                      X,
                                      X_holdout,
                                      landmark_times,
                                      approx_times,
                                      nuisance,
                                      fold_ID){
  if (nuisance == "survSL"){
    event.SL.library <- cens.SL.library <- c("survSL.km",
                                             "survSL.coxph",
                                             "survSL.expreg",
                                             "survSL.weibreg",
                                             "survSL.loglogreg",
                                             "survSL.AFTreg",
                                             "survSL.rfsrc")

    fit <- survSuperLearner::survSuperLearner(time = time,
                                              event = event,
                                              X = X,
                                              newX = rbind(X_holdout, X),
                                              new.times = approx_times,
                                              event.SL.library = event.SL.library,
                                              cens.SL.library = cens.SL.library,
                                              verbose = TRUE,
                                              obsWeights = NULL,
                                              control = list(initWeightAlg = "survSL.rfsrc",
                                                             initWeight = "censoring",
                                                             max.SL.iter = 20,
                                                             event.t.grid = approx_times[approx_times <= max(time)],
                                                             cens.t.grid = approx_times[approx_times <= max(time)]),
                                              cvControl = list(V = 5))
    S_hat_train <- fit$event.SL.predict[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
    G_hat_train <- fit$cens.SL.predict[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
    S_hat <- fit$event.SL.predict[1:nrow(X_holdout),]
    G_hat <- fit$cens.SL.predict[1:nrow(X_holdout),]
    f_hat <- S_hat[,which(approx_times %in% landmark_times),drop=FALSE]
    f_hat_train <- S_hat_train[,which(approx_times %in% landmark_times),drop=FALSE]
  } else if (nuisance == "rfsrc"){
    event_dat <- data.frame(Y = time, Delta = event, X)
    cens_dat <- data.frame(Y = time, Delta = 1 - event, X)
    tune_grid <- expand.grid(mtry = seq(1, ceiling(sqrt(ncol(X))), by = 1),
                             ntree = c(500, 1000), nodesize = c(5, 15, 25))
    event_perfs <- rep(NA, nrow(tune_grid))
    cens_perfs <- rep(NA, nrow(tune_grid))
    for (i in 1:nrow(tune_grid)){
      event_fit_perf <- randomForestSRC::rfsrc(Surv(Y, Delta) ~ .,
                                               data = event_dat,
                                               ntree = tune_grid$ntree[i],
                                               mtry = tune_grid$mtry[i],
                                               nodesize = tune_grid$nodesize[i],
                                               importance = FALSE)$err.rate[tune_grid$ntree[i]]
      cens_fit_perf <- randomForestSRC::rfsrc(Surv(Y, Delta) ~ .,
                                              data = cens_dat,
                                              nodesize = tune_grid$nodesize[i],
                                              ntree = tune_grid$ntree[i],
                                              mtry = tune_grid$mtry[i],
                                              importance = FALSE)$err.rate[tune_grid$ntree[i]]
      event_perfs[i] <- event_fit_perf
      cens_perfs[i] <- cens_fit_perf
    }
    opt_event_index <- which.min(event_perfs)
    opt_cens_index <- which.min(cens_perfs)
    opt_event_mtry <- tune_grid$mtry[opt_event_index]
    opt_event_ntree <- tune_grid$ntree[opt_event_index]
    opt_event_nodesize <- tune_grid$nodesize[opt_event_index]
    opt_cens_mtry <- tune_grid$mtry[opt_cens_index]
    opt_cens_ntree <- tune_grid$ntree[opt_cens_index]
    opt_cens_nodesize <- tune_grid$nodesize[opt_cens_index]
    event_fit <- randomForestSRC::rfsrc(Surv(Y, Delta) ~ .,
                                        data = event_dat,
                                        ntree = opt_event_ntree,
                                        mtry = opt_event_mtry,
                                        importance = FALSE,
                                        perf.type = "none")
    cens_fit <- randomForestSRC::rfsrc(Surv(Y, Delta) ~ .,
                                       data = cens_dat,
                                       ntree = opt_cens_ntree,
                                       mtry = opt_cens_mtry,
                                       importance = FALSE,
                                       perf.type = "none")
    event_q_pred_train <- predict(event_fit, newdata=X, importance='none')$survival
    event_q_pred <- predict(event_fit, newdata=X_holdout, importance='none')$survival
    cens_q_pred <- predict(cens_fit, newdata=X_holdout, importance='none')$survival
    cens_q_pred_train <- predict(cens_fit, newdata=X, importance='none')$survival
    f_hat <- t(sapply(1:nrow(event_q_pred), function(i) {
      stats::approx(c(0,event_fit$time.interest), c(1,event_q_pred[i,]),
                    method='constant', xout = landmark_times, rule = 2)$y
    }))
    f_hat_train <- t(sapply(1:nrow(event_q_pred_train), function(i) {
      stats::approx(c(0,event_fit$time.interest), c(1,event_q_pred_train[i,]),
                    method='constant', xout = landmark_times, rule = 2)$y
    }))
    S_hat <- t(sapply(1:nrow(event_q_pred), function(i) {
      stats::approx(c(0,event_fit$time.interest), c(1,event_q_pred[i,]),
                    method='constant', xout = approx_times, rule = 2)$y
    }))
    S_hat_train <- t(sapply(1:nrow(event_q_pred_train), function(i) {
      stats::approx(c(0,event_fit$time.interest), c(1,event_q_pred_train[i,]),
                    method='constant', xout = approx_times, rule = 2)$y
    }))
    G_hat <- t(sapply(1:nrow(cens_q_pred), function(i) {
      stats::approx(c(0,cens_fit$time.interest), c(1,cens_q_pred[i,]),
                    method='constant', xout = approx_times, rule = 2)$y
    }))
    G_hat_train <- t(sapply(1:nrow(cens_q_pred_train), function(i) {
      stats::approx(c(0,cens_fit$time.interest), c(1,cens_q_pred_train[i,]),
                    method='constant', xout = approx_times, rule = 2)$y
    }))
  } else if (nuisance == "stackG"){
    tune <- list(ntrees = c(250, 500, 1000),
                 max_depth = c(1, 2),
                 minobspernode = 10,
                 shrinkage = 0.01)
    xgb_grid <- create.SL.xgboost(tune = tune)
    SL.library <- c("SL.mean", "SL.glm", "SL.earth", "SL.gam", "SL.ranger", xgb_grid$names)
    surv_out <- survML::stackG(time = time,
                               event = event,
                               X = X,
                               newX = rbind(X_holdout, X),
                               newtimes = approx_times,
                               time_grid_approx = approx_times,
                               bin_size = 0.05,
                               time_basis = "continuous",
                               surv_form = "PI",
                               SL_control = list(SL.library = SL.library,
                                                 V = 5))
    S_hat <- surv_out$S_T_preds[1:nrow(X_holdout),]
    G_hat <- surv_out$S_C_preds[1:nrow(X_holdout),]
    f_hat <- S_hat[,which(approx_times %in% landmark_times),drop=FALSE]
    S_hat_train <- surv_out$S_T_preds[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
    G_hat_train <- surv_out$S_C_preds[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
    f_hat_train <- S_hat_train[,which(approx_times %in% landmark_times),drop=FALSE]
  }
  # weird things happen with a single landmark time for some algorithms
  if (dim(f_hat)[2] != length(landmark_times)){
    f_hat <- t(f_hat)
    f_hat_train = t(f_hat_train)
  }

  return(list(S_hat = S_hat,
              G_hat = G_hat,
              f_hat = f_hat,
              f_hat_train = f_hat_train,
              S_hat_train = S_hat_train,
              G_hat_train = G_hat_train))
}

generate_DR_predictions <- function(time,
                                    event,
                                    X,
                                    X_holdout,
                                    landmark_times,
                                    approx_times,
                                    S_hat,
                                    G_hat){

  tune <- list(ntrees = c(250, 500, 1000),
               max_depth = c(1, 2),
               minobspernode = 10,
               shrinkage = 0.01)
  xgb_grid <- create.SL.xgboost(tune = tune)
  SL.library <- c("SL.mean", "SL.glm", "SL.earth", "SL.gam", "SL.ranger", xgb_grid$names)

  DR_predictions <- matrix(NA, nrow = nrow(X_holdout), ncol = length(landmark_times))
  DR_predictions_train <- matrix(NA, nrow = nrow(X), ncol = length(landmark_times))
  for (i in 1:length(landmark_times)){
    tau <- landmark_times[i]
    Delta_t <- event * (time <= tau) + (time > tau)
    Y_t <- pmin(time, tau)

    # calculate m_0 at a specific time t
    calc_one <- function(t){
      if (approx_times[t] >= tau){
        m <-  rep(1, nrow(S_hat))
      } else{
        m <- S_hat[,which(approx_times == tau)]/S_hat[,t]
      }
      return(m)
    }
    # calculate m_0 over the whole grid
    ms <- matrix(unlist(lapply(1:length(approx_times), FUN = calc_one)),
                 nrow = length(time))

    m_Y <- apply(X = matrix(1:length(Y_t)), MARGIN = 1,
                 FUN = function(x) ms[x,which.min(abs(approx_times - Y_t[x]))])

    G_hat_Y <- apply(X = matrix(1:length(Y_t)), MARGIN = 1,
                     FUN = function(x) G_hat[x,which.min(abs(approx_times - Y_t[x]))])

    term1 <- Delta_t*(Y_t >= tau) / G_hat_Y
    term2 <- (1 - Delta_t) * m_Y / G_hat_Y

    int.vals <- t(sapply(1:length(time), function(j) {
      vals <- diff(1/G_hat[j,])* ms[j,-ncol(ms)]
      if(any(approx_times[-1] > Y_t[j])){
        vals[approx_times[-1] > Y_t[j]] <- 0
      }
      sum(vals)
    }))

    term3 <- t(int.vals)

    DR_pseudo_outcome <- term1 + term2 - term3

    SL_fit <- SuperLearner::SuperLearner(Y = DR_pseudo_outcome,
                                         X = X,
                                         family = stats::gaussian(),
                                         SL.library = SL.library,
                                         method = "method.NNLS",
                                         verbose = FALSE)
    SL_preds <- predict(SL_fit, newdata = X_holdout)$pred
    DR_predictions[,i] <- SL_preds
    SL_preds_train <- predict(SL_fit, newdata = X)$pred
    DR_predictions_train[,i] <- SL_preds_train

  }
  return(list(f_hat = DR_predictions, f_hat_train = DR_predictions_train))
}

generate_reduced_predictions <- function(f_hat,
                                         X_reduced,
                                         X_reduced_holdout,
                                         setting){

  tune <- list(ntrees = c(250, 500, 1000),
               max_depth = c(1, 2),
               minobspernode = 10, shrinkage = 0.01)
  xgb_grid <- create.SL.xgboost(tune = tune)
  SL.library <- c("SL.mean", "SL.glm", "SL.earth",
                  "SL.gam", "SL.ranger", xgb_grid$names)

  long_dat <- data.frame(f_hat = f_hat, X_reduced)
  long_new_dat <- data.frame(X_reduced_holdout)
  reduced_fit <- SuperLearner::SuperLearner(Y = long_dat$f_hat,
                                            X = long_dat[,2:ncol(long_dat),drop=FALSE],
                                            family = stats::gaussian(),
                                            SL.library = SL.library,
                                            method = "method.NNLS",
                                            verbose = FALSE)
  fs_hat <- matrix(predict(reduced_fit, newdata = long_new_dat)$pred,
                   nrow = nrow(X_reduced_holdout),
                   ncol = 1)

  return(list(fs_hat = fs_hat))
}

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

CV_generate_full_predictions_landmark <- function(time,
                                                  event,
                                                  X,
                                                  landmark_times,
                                                  approx_times,
                                                  nuisance,
                                                  folds,
                                                  sample_split){
  .V <- length(unique(folds))
  CV_full_preds_train <- list()
  CV_full_preds <- list()
  CV_S_preds <- list()
  CV_S_preds_train <- list()
  CV_G_preds <- list()

  if (.V == 2 & sample_split){
    time_train <- time
    event_train <- event
    X_train <- X
    X_holdout <- X
    full_preds <- generate_full_predictions(time = time_train,
                                            event = event_train,
                                            X = X_train,
                                            X_holdout = X_holdout,
                                            landmark_times = landmark_times,
                                            approx_times = approx_times,
                                            nuisance = nuisance)
    for (j in 1:.V){
      CV_full_preds_train[[j]] <- full_preds$f_hat_train
      CV_full_preds[[j]] <- full_preds$f_hat[folds == j,]
      CV_S_preds[[j]] <- full_preds$S_hat[folds == j,]
      CV_G_preds[[j]] <- full_preds$G_hat[folds == j,]
    }

  } else{
    for (j in 1:.V){

      if (.V == 1){ # if not actually cross fitting
        time_train <- time[folds == j]
        event_train <- event[folds == j]
        X_train <- X[folds == j,]
      }else{ # if actually cross fitting
        time_train <- time[folds != j]
        event_train <- event[folds != j]
        X_train <- X[folds != j,]
      }
      X_holdout <- X[folds == j,]
      full_preds <- generate_full_predictions(time = time_train,
                                              event = event_train,
                                              X = X_train,
                                              X_holdout = X_holdout,
                                              landmark_times = landmark_times,
                                              approx_times = approx_times,
                                              nuisance = nuisance)
      CV_full_preds_train[[j]] <- full_preds$f_hat_train
      CV_full_preds[[j]] <- full_preds$f_hat
      CV_S_preds[[j]] <- full_preds$S_hat
      CV_G_preds[[j]] <- full_preds$G_hat
    }
  }

  return(list(CV_full_preds_train = CV_full_preds_train,
              CV_full_preds = CV_full_preds,
              CV_S_preds = CV_S_preds,
              CV_G_preds = CV_G_preds))
}

CV_generate_reduced_predictions_landmark <- function(time,
                                                     event,
                                                     X,
                                                     landmark_times,
                                                     folds,
                                                     indx,
                                                     sample_split,
                                                     full_preds_train){
  .V <- length(unique(folds))
  CV_reduced_preds <- list()

  if (.V == 2 & sample_split){
    X_reduced_train <- X[,-indx,drop=FALSE]
    X_reduced_holdout <- X[,-indx,drop=FALSE]
    preds_j <- matrix(NA, nrow = nrow(X_reduced_holdout), ncol = length(landmark_times))
    for (t in landmark_times){
      # in this situation full_preds_train[[1]] and [[2]] are identical (no crossfitting)
      outcomes <- full_preds_train[[1]][,which(landmark_times == t)]
      reduced_preds <- generate_reduced_predictions_regress(f_hat = outcomes,
                                                            X_reduced = X_reduced_train,
                                                            X_reduced_holdout = X_reduced_holdout,
                                                            landmark_times = t,
                                                            setting = "landmark")

      preds_j[,which(landmark_times == t)] <- reduced_preds$fs_hat
    }
    for (j in 1:.V){
      CV_reduced_preds[[j]] <- preds_j[folds == j,]
    }
  } else{
    for (j in 1:.V){
      if (.V == 1){ # if not actually cross fitting
        X_reduced_train <- X[folds == j,-indx,drop=FALSE]
      } else{ # if actually cross fitting
        X_reduced_train <- X[folds != j,-indx,drop=FALSE]
      }

      X_reduced_holdout <- X[folds == j,-indx,drop=FALSE]

      preds_j <- matrix(NA, nrow = nrow(X_reduced_holdout), ncol = length(landmark_times))
      for (t in landmark_times){
        outcomes <- full_preds_train[[j]][,which(landmark_times == t)]
        reduced_preds <- generate_reduced_predictions_regress(f_hat = outcomes,
                                                              X_reduced = X_reduced_train,
                                                              X_reduced_holdout = X_reduced_holdout,
                                                              landmark_times = t,
                                                              setting = "landmark")
        preds_j[,which(landmark_times == t)] <- reduced_preds$fs_hat
      }
      CV_reduced_preds[[j]] <- preds_j
    }
  }
  return(CV_reduced_preds)
}

