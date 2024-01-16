do_one <- function(n_train,
                   misspec_type,
                   robust){

  start <- Sys.time()

  nuisance <- "survSL"
  crossfit <- TRUE
  landmark_times <- c(0.5, 0.9)

  # training data
  train <- generate_data(n = n_train, scenario = "1", sdy = 1)

  sample_split <- FALSE
  dimension <- 2

  time <- train$y
  event <- train$delta
  X <- train[,1:dimension]
  approx_times <- sort(unique(c(0, time[event == 1], landmark_times)))
  approx_times <- approx_times[approx_times <= max(landmark_times)]

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
                                                            misspec_type = misspec_type,
                                                            DR_f0 = robust)

  CV_full_preds <- V0_preds$CV_full_preds
  CV_full_preds_train <- V0_preds$CV_full_preds_train
  CV_S_preds <- V0_preds$CV_S_preds
  CV_G_preds <- V0_preds$CV_G_preds

  V0_preds <- CV_generate_reduced_predictions_landmark(time = time,
                                                       event = event,
                                                       X = X,
                                                       landmark_times = landmark_times,
                                                       folds = folds,
                                                       indx = 1,
                                                       sample_split = sample_split,
                                                       full_preds_train = CV_full_preds_train)

  CV_reduced_preds <- V0_preds

  output <- survML::vim_brier(time = train$y,
                            event = train$delta,
                            approx_times = approx_times,
                            landmark_times = landmark_times,
                            f_hat = CV_full_preds,
                            fs_hat = CV_reduced_preds,
                            S_hat = CV_S_preds,
                            G_hat = CV_G_preds,
                            folds = folds,
                            ss_folds = ss_folds,
                            sample_split = sample_split)

  output$n_train <- n_train
  output$vim <- "brier"
  output$indx <- "1"
  output$misspec_type <- misspec_type
  output$robust <- robust
  end <- Sys.time()
  runtime <- difftime(end, start, units = "mins")
  output$runtime <- runtime
  return(output)
}
