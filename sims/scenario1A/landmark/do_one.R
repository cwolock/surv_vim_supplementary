do_one <- function(n_train,
                   nuisance,
                   crossfit){

  start <- Sys.time()

  vims <- c("AUC", "brier")
  landmark_times <- c(0.5, 0.9)

  # training data
  train <- generate_data(n = n_train, scenario = "1A", sdy = 1)

  sample_split <- FALSE
  dimension <- 2
  indxs <- c("1", "2")

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

  V0_preds <- CV_generate_full_predictions_landmark(time = time,
                                                    event = event,
                                                    X = X,
                                                    landmark_times = landmark_times,
                                                    approx_times = approx_times,
                                                    nuisance = nuisance,
                                                    folds = folds,
                                                    sample_split = sample_split)

  CV_full_preds <- V0_preds$CV_full_preds
  CV_full_preds_train <- V0_preds$CV_full_preds_train
  CV_S_preds <- V0_preds$CV_S_preds
  CV_G_preds <- V0_preds$CV_G_preds

  shared_settings <- expand.grid(indx = indxs, vim = vims)

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

    for (k in 1:length(vims)){
      vim <- vims[k]
      if (vim == "brier"){
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
      } else if (vim == "AUC"){
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
      } else if (vim == "rsquared"){
        output <- survML::vim_rsquared(time = train$y,
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
      }
      output$vim <- vim
      output$indx <- char_indx
      if (!(i == 1 & k == 1)){
        pooled_output <- rbind(pooled_output, output)
      } else{
        pooled_output <- output
      }
    }
  }
  end <- Sys.time()
  runtime <- difftime(end, start, units = "mins")
  pooled_output$runtime <- runtime
  pooled_output$crossfit <- crossfit
  pooled_output$n_train <- n_train
  pooled_output$nuisance <- nuisance
  return(pooled_output)
}
