do_one <- function(n_train,
                   crossfit){

  start <- Sys.time()

  tau <- 0.9

  # training data
  train <- generate_data(n = n_train, scenario = "3", sdy = 1)

  sample_split <- FALSE
  dimension <- 5
  indxs <- c("1", "2")

  time <- train$y
  event <- train$delta
  X <- train[,1:dimension]
  approx_times <- sort(unique(c(0, time[event == 1], tau)))
  approx_times <- approx_times[approx_times <= tau]

  tune <- list(ntrees = c(250, 500, 1000),
               max_depth = c(1, 2),
               minobspernode = 10,
               shrinkage = 0.01)
  xgb_grid <- create.SL.xgboost(tune = tune)
  SL.library <- c("SL.mean", "SL.glm.interaction", "SL.earth", "SL.gam", "SL.ranger", xgb_grid$names)

  for (i in 1:length(indxs)){
    char_indx <- as.character(indxs[i])
    indx <- as.numeric(strsplit(char_indx, split = ",")[[1]])
    if (i == 1){
      output <- survML::vim(type = "survival_time_MSE",
                            time = time,
                            event = event,
                            X = X,
                            restriction_time = tau,
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
    } else if (i == 2 ){
      output <- survML::vim(type = "survival_time_MSE",
                            time = time,
                            event = event,
                            X = X,
                            restriction_time = tau,
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
    }
  }

  end <- Sys.time()
  runtime <- difftime(end, start, units = "mins")
  pooled_output$runtime <- runtime
  pooled_output$crossfit <- crossfit
  pooled_output$n_train <- n_train
  return(pooled_output)
}
