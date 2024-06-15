options(expressions= 20000) # this is necessary for the `combinations` function to work with large n
do_one <- function(n_train,
                   nuisance,
                   crossfit){
  start <- Sys.time()

  # training data
  train <- generate_data(n = n_train, scenario = "4A", sdy = 1)

  sample_split <- TRUE
  dimension <- 25
  indxs <- c("1", "4", "1,4")

  time <- train$y
  event <- train$delta
  X <- train[,1:dimension]

  tau <- 0.9
  approx_times <- sort(unique(c(0, time[time <= tau & event == 1])))

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

  nuisance_preds <- CV_generate_full_predictions_landmark(time = time,
                                                          event = event,
                                                          X = X,
                                                          landmark_times = tau,
                                                          approx_times = approx_times,
                                                          nuisance = nuisance,
                                                          folds = folds,
                                                          sample_split = sample_split)
  CV_S_preds <- nuisance_preds$CV_S_preds
  CV_S_preds_train <- nuisance_preds$CV_S_preds_train
  CV_G_preds <- nuisance_preds$CV_G_preds
  
  V0_preds <- CV_generate_predictions_cindex(time = time,
                                             event = event,
                                             X = X,
                                             approx_times = approx_times,
                                             folds = folds,
                                             sample_split = sample_split,
                                             CV_S_preds_train =  CV_S_preds_train,
                                             CV_S_preds = CV_S_preds,
                                             indx = NULL,
                                             subsample_n = 500,
                                             params =  list(
                                               mstop = c(100, 250, 500),
                                               nu = c(0.1),
                                               sigma = c(0.01, 0.05),
                                               learner = c("glm")))
  
  CV_full_preds <- V0_preds

  for (i in 1:length(indxs)){
    char_indx <- as.character(indxs[i])
    indx <- as.numeric(strsplit(char_indx, split = ",")[[1]])

    V0_preds <- CV_generate_predictions_cindex(time = time,
                                               event = event,
                                               X = X,
                                               approx_times = approx_times,
                                               folds = folds,
                                               sample_split = sample_split,
                                               indx = indx,
                                               CV_S_preds_train =  CV_S_preds_train,
                                               CV_S_preds = CV_S_preds,
                                               subsample_n = 500,
                                               params =  list(
                                                 mstop = c(100, 250, 500),
                                                 nu = c(0.1),
                                                 sigma = c(0.01, 0.05),
                                                 learner = c("glm")))
    
    CV_reduced_preds <- V0_preds

    output <- survML::vim_cindex(time = train$y,
                                  event = train$delta,
                                  approx_times = approx_times,
                                  f_hat = CV_full_preds,
                                  fs_hat = CV_reduced_preds,
                                  S_hat = CV_S_preds,
                                  G_hat = CV_G_preds,
                                  folds = folds,
                                  sample_split = sample_split,
                                  ss_folds = ss_folds,
                                  tau = tau)

    output$indx <- rep(char_indx, nrow(output))

    if (i != 1){
      pooled_output <- rbind(pooled_output, output)
    } else{
      pooled_output <- output
    }
  }

  end <- Sys.time()
  runtime <- difftime(end, start, units = "mins")
  pooled_output$runtime <- runtime
  pooled_output$vim <- "cindex"
  pooled_output$crossfit <- crossfit
  pooled_output$n_train <- n_train
  pooled_output$nuisance <- nuisance
  return(pooled_output)
}
