do_one <- function(n_train,
                   nuisance,
                   crossfit){

  start <- Sys.time()

  vims <- c("AUC", "Brier")
  landmark_times <- c(0.5, 0.9)

  # training data
  train <- generate_data(n = n_train, scenario = "3", sdy = 1)

  sample_split <- FALSE
  dimension <- 5
  indxs <- c("1", "2")

  time <- train$y
  event <- train$delta
  X <- train[,1:dimension]
  approx_times <- sort(unique(c(0, time[event == 1], landmark_times)))
  approx_times <- approx_times[approx_times <= max(landmark_times)]

  tune <- list(ntrees = c(250, 500, 1000),
               max_depth = c(1, 2),
               minobspernode = 10,
               shrinkage = 0.01)
  xgb_grid <- create.SL.xgboost(tune = tune)
  SL.library <- c("SL.mean", "SL.glm.interaction", "SL.earth", "SL.gam", "SL.ranger", xgb_grid$names)

  # cf_fold_num <- switch((crossfit) + 1, 1, 5)
  # ss_fold_num <- 2*cf_fold_num

  # V <- switch((sample_split) + 1, cf_fold_num, ss_fold_num)

  # folds <- sample(rep(seq_len(V), length = length(time))) # 2V of them

  # if (sample_split){
  #   ss_folds <- c(rep(1, V/2), rep(2, V/2))
  # } else{
  #   ss_folds <- rep(1, V)
  # }
  #
  # ss_folds <- as.numeric(folds %in% which(ss_folds == 2))
  #
  # V0_preds <- CV_generate_full_predictions_landmark(time = time,
  #                                                   event = event,
  #                                                   X = X,
  #                                                   landmark_times = landmark_times,
  #                                                   approx_times = approx_times,
  #                                                   nuisance = nuisance,
  #                                                   folds = folds,
  #                                                   sample_split = sample_split)
  #
  # CV_full_preds <- V0_preds$CV_full_preds
  # CV_full_preds_train <- V0_preds$CV_full_preds_train
  # CV_S_preds <- V0_preds$CV_S_preds
  # CV_G_preds <- V0_preds$CV_G_preds

  shared_settings <- expand.grid(indx = indxs, vim = vims)

  for (i in 1:length(indxs)){
    char_indx <- as.character(indxs[i])
    indx <- as.numeric(strsplit(char_indx, split = ",")[[1]])

    # V0_preds <- CV_generate_reduced_predictions_landmark(time = time,
    #                                                      event = event,
    #                                                      X = X,
    #                                                      landmark_times = landmark_times,
    #                                                      folds = folds,
    #                                                      indx = indx,
    #                                                      sample_split = sample_split,
    #                                                      full_preds_train = CV_full_preds_train)
    # CV_reduced_preds <- V0_preds

    for (k in 1:length(vims)){
      vim <- vims[k]

      if (i == 1 & k == 1){
        output <- vim(type = vim,
                      time = time,
                      event = event,
                      X = X,
                      landmark_times = landmark_times,
                      approx_times = approx_times,
                      large_feature_vector = 1:ncol(X),
                      small_feature_vector = (1:ncol(X))[-as.numeric(indx)],
                      conditional_surv_generator_control = list(SL.library = SL.library),
                      large_oracle_generator_control = list(SL.library = SL.library),
                      small_oracle_generator_control = list(SL.library = SL.library),
                      cf_fold_num = ifelse(crossfit, 5, 1),
                      sample_split = FALSE,
                      scale_est = FALSE)
        saved_conditional_surv_preds <- output$conditional_surv_preds
        saved_large_oracle_preds <- output$large_oracle_preds
        saved_small_oracle_preds <- output$small_oracle_preds
        saved_folds <- output$folds
        output$result$indx <- char_indx
        pooled_output <- output$result
      } else if (i == 2 & k == 1){
        output <- vim(type = vim,
                      time = time,
                      event = event,
                      X = X,
                      landmark_times = landmark_times,
                      approx_times = approx_times,
                      large_feature_vector = 1:ncol(X),
                      small_feature_vector = (1:ncol(X))[-as.numeric(indx)],
                      conditional_surv_preds = saved_conditional_surv_preds,
                      large_oracle_preds = saved_large_oracle_preds,
                      small_oracle_generator_control = list(SL.library = SL.library),
                      cf_folds = saved_folds$cf_folds,
                      ss_folds = saved_folds$ss_folds,
                      sample_split = FALSE,
                      scale_est = FALSE)
        saved_small_oracle_preds <- output$small_oracle_preds
        output$result$indx <- char_indx
        pooled_output <- rbind(pooled_output, output$result)
      } else{
        output <- vim(type = vim,
                      time = time,
                      event = event,
                      X = X,
                      landmark_times = landmark_times,
                      approx_times = approx_times,
                      large_feature_vector = 1:ncol(X),
                      small_feature_vector = (1:ncol(X))[-as.numeric(indx)],
                      conditional_surv_preds = saved_conditional_surv_preds,
                      large_oracle_preds = saved_large_oracle_preds,
                      small_oracle_preds = saved_small_oracle_preds,
                      cf_folds = saved_folds$cf_folds,
                      ss_folds = saved_folds$ss_folds,
                      sample_split = FALSE,
                      scale_est = FALSE)
        output$result$indx <- char_indx
        pooled_output <- rbind(pooled_output, output$result)
      }
    }
  }
  end <- Sys.time()
  runtime <- difftime(end, start, units = "mins")
  pooled_output$nuisance <- nuisance
  pooled_output$runtime <- runtime
  pooled_output$crossfit <- crossfit
  pooled_output$n_train <- n_train
  return(pooled_output)
}
