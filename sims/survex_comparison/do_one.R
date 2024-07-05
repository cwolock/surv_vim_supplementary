do_one <- function(n_train,
                   method,
                   correlation){

  start <- Sys.time()
  landmark_times <- c(0.5, 0.9)

  # training data
  if (correlation){
    scenario <- "1B"
  } else{
    scenario <- "1C"
  }
  train <- generate_data(n = n_train, scenario = scenario, sdy = 1)
  train <- train %>% select(-c(t, c))

  dimension <- 4
  time <- train$y
  event <- train$delta
  X <- train[,1:dimension]
  approx_times <- sort(unique(c(0, time[event == 1], landmark_times)))
  approx_times <- approx_times[approx_times <= max(landmark_times)]

  if (method == "permutation"){

    event_dat <- data.frame(Y = time, Delta = event, X)
    tune_grid <- expand.grid(mtry = seq(1, ceiling(sqrt(ncol(X))), by = 1),
                             ntree = c(500, 1000), nodesize = c(5, 15, 25))
    event_perfs <- rep(NA, nrow(tune_grid))
    for (i in 1:nrow(tune_grid)){
      event_fit_perf <- randomForestSRC::rfsrc(Surv(Y, Delta) ~ .,
                                               data = event_dat,
                                               ntree = tune_grid$ntree[i],
                                               mtry = tune_grid$mtry[i],
                                               nodesize = tune_grid$nodesize[i],
                                               importance = FALSE)$err.rate[tune_grid$ntree[i]]
      event_perfs[i] <- event_fit_perf
    }
    opt_event_index <- which.min(event_perfs)
    opt_event_mtry <- tune_grid$mtry[opt_event_index]
    opt_event_ntree <- tune_grid$ntree[opt_event_index]
    opt_event_nodesize <- tune_grid$nodesize[opt_event_index]
    rsf <- randomForestSRC::rfsrc(Surv(Y, Delta) ~ .,
                                  data = event_dat,
                                  ntree = opt_event_ntree,
                                  mtry = opt_event_mtry,
                                  nodesize = opt_event_nodesize,
                                  importance = FALSE,
                                  perf.type = "none")
    # rsf <- randomForestSRC::rfsrc(Surv(y, delta)~., data = train)
    rsf_exp <- survex::explain(rsf)

    # model_parts_rsf <- model_parts(rsf_exp)
    model_parts_rsf_auc  <- survex::model_parts(rsf_exp, loss_function=loss_one_minus_cd_auc)

    # plot(model_parts_rsf_auc)

    res <- model_parts_rsf_auc$result %>% filter(`_permutation_` == 0)
    time_indices <- sapply(landmark_times,
                           FUN = function(x) which.min(abs(x - res$`_times_`)))
    pooled_output <- res[time_indices,] %>%  mutate(landmark_time = c(0.5, 0.9)) %>%
      select(landmark_time, x1, x2, x3, x4) %>%
      pivot_longer(cols = -landmark_time, names_to = "indx", values_to = "est") %>%
      mutate(indx = case_when(indx == "x1" ~ "1",
                              indx == "x2" ~ "2",
                              indx == "x3" ~ "3",
                              indx == "x4" ~ "4"))


    # event.SL.library <- cens.SL.library <- lapply(c("survSL.km", "survSL.coxph",
    #                                                 "survSL.expreg.int", "survSL.weibreg.int",
    #                                                 "survSL.loglogreg.int", "survSL.AFTreg.int",
    #                                                 "survSL.rfsrc"), function(alg) {
    #                                                   c(alg, "survscreen.marg")
    #                                                 })
    #
    # event.SL.library <- cens.SL.library <- c("survSL.km", "survSL.coxph",
    #                                                 "survSL.expreg.int", "survSL.weibreg.int",
    #                                                 "survSL.loglogreg.int", "survSL.AFTreg.int",
    #                                                 "survSL.rfsrc")
    #
    # fit <- survSuperLearner::survSuperLearner(time = train$y,
    #                                           event = train$delta,
    #                                           X = train[,1:dimension],
    #                                           newX = train[,1:dimension],
    #                                           new.times = landmark_times,
    #                                           event.SL.library = event.SL.library,
    #                                           cens.SL.library = cens.SL.library,
    #                                           verbose = FALSE,
    #                                           obsWeights = NULL,
    #                                           control = list(initWeightAlg = "survSL.rfsrc",
    #                                                          initWeight = "censoring",
    #                                                          max.SL.iter = 20,
    #                                                          event.t.grid = approx_times[approx_times <= max(time)],
    #                                                          cens.t.grid = approx_times[approx_times <= max(time)]),
    #                                           cvControl = list(V = 5))
    #
    # preds <- predict.survSuperLearner(fit, newdata = X[1:10,], new.times = landmark_times)
    # # S_hat_train <- fit$event.SL.predict[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
    # # G_hat_train <- fit$cens.SL.predict[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
    # # S_hat <- fit$event.SL.predict[1:nrow(X_holdout),]
    #
    #
    # # risk_pred <- function(model, newdata) predict(model, newdata, type = "risk")
    # surv_pred <- function(model, newdata, times) predict.survSuperLearner(model,
    #                                                                       newdata = newdata,
    #                                                                       new.times = times)$event.SL.predict
    # # chf_pred <- function(model, newdata, times) -log(surv_pred(model, newdata, times))
    #
    # manual_cph_explainer <- explain_survival(model = fit,
    #                                          data = train[,1:dimension],
    #                                          y = Surv(train$y, train$delta),
    #                                          # predict_function = risk_pred,
    #                                          predict_survival_function = surv_pred,
    #                                          # predict_cumulative_hazard_function = chf_pred,
    #                                          label="manual survSL")

    # model_parts_survSL <- model_parts(manual_cph_explainer)

    # plot(model_parts_survSL)
  } else if (method == "intrinsic"){

    sample_split <- FALSE
    crossfit <- TRUE
    indxs <- c("1", "2", "3", "4")

    cf_fold_num <- switch((crossfit) + 1, 1, 5)
    ss_fold_num <- 2*cf_fold_num

    V <- switch((sample_split) + 1, cf_fold_num, ss_fold_num)

    folds <- sample(rep(seq_len(V), length = length(time))) # 2V of them

    if (sample_split){
      ss_folds <- c(rep(1, V/2), rep(2, V/2))
    } else{
      ss_folds <- rep(1, V)
    }

    ss_folds <- as.numeric(folds %in% which(ss_folds == 2))

    V0_preds <- CV_generate_full_predictions_landmark(time = time,
                                                      event = event,
                                                      X = X,
                                                      landmark_times = landmark_times,
                                                      approx_times = approx_times,
                                                      nuisance = "rfsrc",
                                                      folds = folds,
                                                      sample_split = sample_split)

    CV_full_preds <- V0_preds$CV_full_preds
    CV_full_preds_train <- V0_preds$CV_full_preds_train
    CV_S_preds <- V0_preds$CV_S_preds
    CV_G_preds <- V0_preds$CV_G_preds

    for (i in 1:length(indxs)){
      char_indx <- as.character(indxs[i])
      indx <- as.numeric(strsplit(char_indx, split = ",")[[1]])

      V0_preds <- CV_generate_reduced_predictions_landmark(time = time,
                                                           event = event,
                                                           X = X,
                                                           landmark_times = landmark_times,
                                                           folds = folds,
                                                           indx = indx,
                                                           sample_split = sample_split,
                                                           full_preds_train = CV_full_preds_train)
      CV_reduced_preds <- V0_preds

      output <- survML::vim_AUC(time = train$y,
                                event = train$delta,
                                approx_times = approx_times,
                                landmark_times = landmark_times,
                                f_hat = lapply(CV_full_preds, function(x) 1-x),
                                fs_hat = lapply(CV_reduced_preds, function(x) 1-x),
                                S_hat = CV_S_preds,
                                G_hat = CV_G_preds,
                                folds = folds,
                                ss_folds = ss_folds,
                                sample_split = sample_split)
      output$indx <- char_indx
      if (!(i == 1)){
        pooled_output <- rbind(pooled_output, output)
      } else{
        pooled_output <- output
      }
    }
  }
  pooled_output <- pooled_output %>% select(landmark_time, indx, est)
  pooled_output$n_train <- n_train
  pooled_output$method <- method
  out1 <- pooled_output %>% filter(landmark_time == 0.5)
  out1$rank <- rank(-out1$est)
  out1$true_rank <- c(2,1,3,4)
  out1$all_right <- all(out1$rank == out1$true_rank)

  out2 <- pooled_output %>% filter(landmark_time == 0.9)
  out2$rank <- rank(-out2$est)
  out2$true_rank <- c(2,1,3,4)
  out2$all_right <- all(out2$rank == out2$true_rank)

  pooled_output <- rbind(out1, out2)

  end <- Sys.time()
  pooled_output$runtime <- difftime(end, start, units = "mins")
  return(pooled_output)
}
