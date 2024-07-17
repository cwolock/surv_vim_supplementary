do_one <- function(n_train,
                   misspec_type){

  start <- Sys.time()

  nuisance <- "survSL"
  crossfit <- TRUE
  landmark_times <- c(0.5)

  # training data
  train <- generate_data(n = n_train, scenario = "4B", sdy = 1)

  sample_split <- TRUE
  dimension <- 25

  time <- train$y
  event <- train$delta
  X <- train[,1:dimension]
  approx_times <- sort(unique(c(0, time[event == 1], landmark_times)))
  approx_times <- approx_times[approx_times <= 0.9]#max(landmark_times)]

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

  V0_preds <- CV_generate_full_predictions_landmark_misspec(time = time,
                                                            event = event,
                                                            X = X,
                                                            landmark_times = landmark_times,
                                                            approx_times = approx_times,
                                                            nuisance = nuisance,
                                                            folds = folds,
                                                            sample_split = sample_split,
                                                            misspec_type = misspec_type)

  CV_full_preds_robust <- V0_preds$CV_full_preds_robust
  CV_full_preds_train_robust <- V0_preds$CV_full_preds_train_robust
  CV_full_preds_nonrobust <- V0_preds$CV_full_preds_nonrobust
  CV_full_preds_train_nonrobust <- V0_preds$CV_full_preds_train_nonrobust
  CV_S_preds <- V0_preds$CV_S_preds
  CV_G_preds <- V0_preds$CV_G_preds

  V0_preds <- CV_generate_reduced_predictions_landmark(time = time,
                                                       event = event,
                                                       X = X,
                                                       landmark_times = landmark_times,
                                                       folds = folds,
                                                       indx = c(1,6),
                                                       sample_split = sample_split,
                                                       full_preds_train = CV_full_preds_train_robust)

  CV_reduced_preds_robust <- V0_preds

  V0_preds <- CV_generate_reduced_predictions_landmark(time = time,
                                                       event = event,
                                                       X = X,
                                                       landmark_times = landmark_times,
                                                       folds = folds,
                                                       indx = c(1,6),
                                                       sample_split = sample_split,
                                                       full_preds_train = CV_full_preds_train_nonrobust)

  CV_reduced_preds_nonrobust <- V0_preds

  output_robustV_robustf <- survML::vim_AUC(time = train$y,
                            event = train$delta,
                            approx_times = approx_times,
                            landmark_times = landmark_times,
                            f_hat = lapply(CV_full_preds_robust, function(x) 1-x),
                            fs_hat = lapply(CV_reduced_preds_robust, function(x) 1-x),
                            S_hat = CV_S_preds,
                            G_hat = CV_G_preds,
                            folds = folds,
                            ss_folds = ss_folds,
                            sample_split = sample_split,
                            robust = TRUE)
  output_robustV_robustf$robust_V <- TRUE
  output_robustV_robustf$robust_f <- TRUE

  output_nonrobustV_robustf <- survML::vim_AUC(time = train$y,
                                   event = train$delta,
                                   approx_times = approx_times,
                                   landmark_times = landmark_times,
                                   f_hat = lapply(CV_full_preds_robust, function(x) 1-x),
                                   fs_hat = lapply(CV_reduced_preds_robust, function(x) 1-x),
                                   S_hat = CV_S_preds,
                                   G_hat = CV_G_preds,
                                   folds = folds,
                                   ss_folds = ss_folds,
                                   sample_split = sample_split,
                                   robust = FALSE)
  output_nonrobustV_robustf$robust_V <- FALSE
  output_nonrobustV_robustf$robust_f <- TRUE

  output_robustV_nonrobustf <- survML::vim_AUC(time = train$y,
                                            event = train$delta,
                                            approx_times = approx_times,
                                            landmark_times = landmark_times,
                                            f_hat = lapply(CV_full_preds_nonrobust, function(x) 1-x),
                                            fs_hat = lapply(CV_reduced_preds_nonrobust, function(x) 1-x),
                                            S_hat = CV_S_preds,
                                            G_hat = CV_G_preds,
                                            folds = folds,
                                            ss_folds = ss_folds,
                                            sample_split = sample_split,
                                            robust = TRUE)
  output_robustV_nonrobustf$robust_V <- TRUE
  output_robustV_nonrobustf$robust_f <- FALSE

  output_nonrobustV_nonrobustf <- survML::vim_AUC(time = train$y,
                                               event = train$delta,
                                               approx_times = approx_times,
                                               landmark_times = landmark_times,
                                               f_hat = lapply(CV_full_preds_nonrobust, function(x) 1-x),
                                               fs_hat = lapply(CV_reduced_preds_nonrobust, function(x) 1-x),
                                               S_hat = CV_S_preds,
                                               G_hat = CV_G_preds,
                                               folds = folds,
                                               ss_folds = ss_folds,
                                               sample_split = sample_split,
                                               robust = FALSE)
  output_nonrobustV_nonrobustf$robust_V <- FALSE
  output_nonrobustV_nonrobustf$robust_f <- FALSE

  output <- rbind(output_robustV_robustf, output_nonrobustV_robustf,
                  output_robustV_nonrobustf, output_nonrobustV_nonrobustf)

  output$n_train <- n_train
  output$vim <- "AUC"
  output$indx <- "1,6"
  output$misspec_type <- misspec_type
  end <- Sys.time()
  runtime <- difftime(end, start, units = "mins")
  output$runtime <- runtime
  return(output)
}
