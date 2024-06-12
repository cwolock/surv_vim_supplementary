do_one <- function(seed,
                   global_seed,
                   approach){
  
  sample_split <- TRUE
  nfolds <- 5
  crossfit <- TRUE
  nuisance <- "stackG"
  
  dat <- readRDS("/home/cwolock/surv_vim_supplementary/data_analysis/cleaned_RS_data_combined36.rds")
  # dat <- readRDS("/Users/cwolock/Dropbox/UW/DISSERTATION/surv_vim_supplementary/data_analysis/cleaned_RS_data_combined36.rds")
  landmark_times <- c(545, 730, 912)
  
  cf_fold_num <- switch((crossfit) + 1, 1, nfolds)
  ss_fold_num <- 2*cf_fold_num
  V <- switch((sample_split) + 1, cf_fold_num, ss_fold_num)
  
  # do stratified cross-fitting/sample-splitting
  male_indices <- which(dat$female == 0)
  female_indices <- which(dat$female == 1)
  
  dat <- dat[c(male_indices, female_indices),]
  
  male_folds <- sample(rep(seq_len(V), length = length(male_indices))) # 2V of them
  female_folds <- sample(rep(seq_len(V), length = length(female_indices))) # 2V of them
  folds <- c(male_folds, female_folds)
  if (sample_split){
    ss_folds <- c(rep(1, V/2), rep(2, V/2))
  } else{
    ss_folds <- rep(1, V)
  }
  male_ss_folds <- as.numeric(male_folds %in% which(ss_folds == 2))
  female_ss_folds <- as.numeric(female_folds %in% which(ss_folds == 2))
  ss_folds <- c(male_ss_folds, female_ss_folds)
  
  
  
  
  # 
  # folds <- sample(rep(seq_len(V), length = nrow(dat))) # 2V of them
  # if (sample_split){
  #   ss_folds <- c(rep(1, V/2), rep(2, V/2))
  # } else{
  #   ss_folds <- rep(1, V)
  # }
  # ss_folds <- as.numeric(folds %in% which(ss_folds == 2))
  
  time <- dat$HIV36fu
  event <- dat$HIV36
  X <- dat %>% select(-c(HIV36fu, HIV36))
  
  X_names <- names(X)
  vacc_index <- which(X_names == "VACC")
  age_index <- which(X_names == "AGE")
  geo_index <- which(X_names == "REG" | X_names == "REG1" | X_names == "REG2")
  bmi_index <- which(X_names == "BMIGE25")
  sexhealth_index <- which(X_names %in% c("ANYSTI", "GENSOR", "GENDIS"))
  sexbehave_index <- which(X_names %in% c("HETERO", "CONDOM", "ANALSX",
                                          "SXHIVP", "EXCHSX", "USXALC",
                                          "MAINPRT", "USXHIVPyes",
                                          "LIVWPRTyes", "OTHPRTyes"))
  social_index <- which(X_names %in% c("URBAN", "FMLDWELL", "HOMESRV"))
  brs_index <- which(X_names != "female" & X_names != "AGE" & X_names != "REG" &
                       X_names != "BMIGE25" & X_names != "VACC" & X_names != "REG1" &
                       X_names != "REG2")
  all_index <- 1:length(X_names)
  
  age_index <- c(age_index, vacc_index)
  geo_index <- c(geo_index, vacc_index)
  bmi_index <- c(bmi_index, vacc_index)
  sexhealth_index <- c(sexhealth_index, vacc_index)
  sexbehave_index <- c(sexbehave_index, vacc_index)
  social_index <- c(social_index, vacc_index)
  brs_index <- c(brs_index, vacc_index)
  
  all_but_geo_index <- unique(c(brs_index, age_index, bmi_index, vacc_index))
  
  vacc_index_text <- paste(vacc_index, collapse = ",")
  age_index_text <- paste(age_index, collapse = ",")
  geo_index_text <- paste(geo_index, collapse = ",")
  bmi_index_text <- paste(bmi_index, collapse = ",")
  sexhealth_index_text <- paste0(sexhealth_index, collapse = ",")
  sexbehave_index_text <- paste0(sexbehave_index, collapse = ",")
  social_index_text <- paste0(social_index, collapse = ",")
  brs_variables_text <- as.character(brs_index)
  brs_index_text <- paste(brs_index, collapse = ",")
  all_but_geo_text <- paste(all_but_geo_index, collapse = ",")
  
  if (approach == "marginal"){
    
    all_index_text <- c(age_index_text,
                        bmi_index_text,
                        sexhealth_index_text, 
                        sexbehave_index_text,
                        social_index_text)
    all_index_names <- c("age",
                         "bmi", "BRS_sexhealth",
                         "BRS_sexbehave", "BRS_social")
  } else if (approach == "conditional"){
    
    all_index_text <- c(geo_index_text,
                        age_index_text,
                        bmi_index_text,
                        sexhealth_index_text, 
                        sexbehave_index_text,
                        social_index_text)
    all_index_names <- c("geo", "age",
                         "bmi", "BRS_sexhealth",
                         "BRS_sexbehave", "BRS_social")
  }
  
  approx_times <- sort(unique(c(time[event == 1], landmark_times)))
  approx_times <- approx_times[approx_times <= max(landmark_times)]
  
  V0_preds <- CV_generate_full_predictions_landmark(time = time,
                                                    event = event,
                                                    X = X,
                                                    landmark_times = landmark_times,
                                                    approx_times = approx_times,
                                                    nuisance = nuisance,
                                                    folds = folds,
                                                    sample_split = sample_split)
  
  CV_full_preds_landmark_train <- V0_preds$CV_full_preds_train
  CV_S_preds <- V0_preds$CV_S_preds
  CV_S_preds_train <- V0_preds$CV_S_preds_train
  CV_G_preds <- V0_preds$CV_G_preds
  
  # now that nuisances have been estimated, drop females
  dat_male <- dat %>% filter(female == 0) %>% select(-c(female))
  time_male <- dat_male$HIV36fu
  event_male <- dat_male$HIV36
  
  dat_female <- dat %>% filter(female == 1) %>% select(-c(female))
  time_female <- dat_female$HIV36fu
  event_female <- dat_female$HIV36
  
  if (approach == "marginal"){
    all_but_geo <- as.numeric(strsplit(all_but_geo_text, split = ",")[[1]])
    # these can be considered "baseline" models - they only use geography
    V0_preds <- CV_generate_reduced_predictions_landmark(time = time,
                                                         event = event,
                                                         X = X,
                                                         landmark_times = landmark_times,
                                                         folds = folds,
                                                         sample_split = sample_split,
                                                         indx = all_but_geo,
                                                         full_preds_train = CV_full_preds_landmark_train)
    
    CV_reduced_preds_landmark <- V0_preds
    
    V0_preds <- CV_generate_predictions_cindex(time = time,
                                               event = event,
                                               X = X,
                                               approx_times = approx_times,
                                               folds = folds,
                                               sample_split = sample_split,
                                               CV_S_preds_train =  CV_S_preds_train,
                                               CV_S_preds = CV_S_preds,
                                               indx = all_but_geo,
                                               subsample_n = 1000,
                                               params =  list(
                                                 mstop = c(100, 250, 500, 1000),
                                                 nu = c(0.1),
                                                 sigma = c(0.01, 0.05),
                                                 learner = c("glm")))
    
    CV_reduced_preds_cindex <- V0_preds
    
  } else if (approach == "conditional"){
    V0_preds <- CV_generate_reduced_predictions_landmark(time = time,
                                                         event = event,
                                                         X = X,
                                                         landmark_times = landmark_times,
                                                         folds = folds,
                                                         sample_split = sample_split,
                                                         indx = vacc_index,
                                                         full_preds_train = CV_full_preds_landmark_train)
    
    CV_full_preds_landmark <- V0_preds
    
    V0_preds <- CV_generate_predictions_cindex(time = time,
                                               event = event,
                                               X = X,
                                               approx_times = approx_times,
                                               folds = folds,
                                               sample_split = sample_split,
                                               CV_S_preds_train =  CV_S_preds_train,
                                               CV_S_preds = CV_S_preds,
                                               indx = vacc_index,
                                               subsample_n = 1000,
                                               params =  list(
                                                 mstop = c(100, 250, 500, 1000),
                                                 nu = c(0.1),
                                                 sigma = c(0.01, 0.05),
                                                 learner = c("glm")))
    
    CV_full_preds_cindex <- V0_preds
  }
  
  for (i in 1:length(all_index_text)){
    char_indx <- as.character(all_index_text[i])
    char_indx_name <- all_index_names[i]
    indx <- as.numeric(strsplit(char_indx, split = ",")[[1]])
    if (approach == "marginal"){
      
      indx <- all_but_geo_index[-which(all_but_geo_index %in% indx)]
      V0_preds <- CV_generate_reduced_predictions_landmark(time = time,
                                                           event = event,
                                                           X = X,
                                                           landmark_times = landmark_times,
                                                           folds = folds,
                                                           sample_split = sample_split,
                                                           indx = indx,
                                                           full_preds_train = CV_full_preds_landmark_train)
      CV_full_preds_landmark <- V0_preds
      
      V0_preds <- CV_generate_predictions_cindex(time = time,
                                                 event = event,
                                                 X = X,
                                                 approx_times = approx_times,
                                                 folds = folds,
                                                 sample_split = sample_split,
                                                 CV_S_preds_train =  CV_S_preds_train,
                                                 CV_S_preds = CV_S_preds,
                                                 indx = indx,
                                                 subsample_n = 1000,
                                                 params =  list(
                                                   mstop = c(100, 250, 500, 1000),
                                                   nu = c(0.1),
                                                   sigma = c(0.01, 0.05),
                                                   learner = c("glm")))
      
      CV_full_preds_cindex <- V0_preds
      
    } else if (approach == "conditional"){
      V0_preds <- CV_generate_reduced_predictions_landmark(time = time,
                                                           event = event,
                                                           X = X,
                                                           landmark_times = landmark_times,
                                                           folds = folds,
                                                           sample_split = sample_split,
                                                           indx = indx,
                                                           full_preds_train = CV_full_preds_landmark_train)
      CV_reduced_preds_landmark <- V0_preds
      
      V0_preds <- CV_generate_predictions_cindex(time = time,
                                                 event = event,
                                                 X = X,
                                                 approx_times = approx_times,
                                                 folds = folds,
                                                 sample_split = sample_split,
                                                 CV_S_preds_train =  CV_S_preds_train,
                                                 CV_S_preds = CV_S_preds,
                                                 indx = indx,
                                                 subsample_n = 1000,
                                                 params =  list(
                                                   mstop = c(100, 250, 500, 1000),
                                                   nu = c(0.1),
                                                   sigma = c(0.01, 0.05),
                                                   learner = c("glm")))
      
      CV_reduced_preds_cindex <- V0_preds
    }
    
    ### MALE
    CV_full_preds_landmark_male <- lapply(1:V,
                                          FUN = function(x){
                                            dat_j <- dat[folds == x,]
                                            males_j <- which(dat_j$female == 0)
                                            return(CV_full_preds_landmark[[x]][males_j,])
                                          })
    CV_reduced_preds_landmark_male <- lapply(1:V,
                                             FUN = function(x){
                                               dat_j <- dat[folds == x,]
                                               males_j <- which(dat_j$female == 0)
                                               return(CV_reduced_preds_landmark[[x]][males_j,])
                                             })
    CV_full_preds_cindex_male <- lapply(1:V,
                                          FUN = function(x){
                                            dat_j <- dat[folds == x,]
                                            males_j <- which(dat_j$female == 0)
                                            return(CV_full_preds_cindex[[x]][males_j])
                                          })
    CV_reduced_preds_cindex_male <- lapply(1:V,
                                             FUN = function(x){
                                               dat_j <- dat[folds == x,]
                                               males_j <- which(dat_j$female == 0)
                                               return(CV_reduced_preds_cindex[[x]][males_j])
                                             })
    CV_S_preds_male <- lapply(1:V,
                              FUN = function(x){
                                dat_j <- dat[folds == x,]
                                males_j <- which(dat_j$female == 0)
                                return(CV_S_preds[[x]][males_j,])
                              })
    CV_G_preds_male <- lapply(1:V,
                              FUN = function(x){
                                dat_j <- dat[folds == x,]
                                males_j <- which(dat_j$female == 0)
                                return(CV_G_preds[[x]][males_j,])
                              })
    
    output_auc_male <- survML::vim_AUC(time = time_male,
                                  event = event_male,
                                  approx_times = approx_times,
                                  landmark_times = landmark_times,
                                  f_hat = lapply(CV_full_preds_landmark_male, function(x) 1-x),
                                  fs_hat = lapply(CV_reduced_preds_landmark_male, function(x) 1-x),
                                  S_hat = CV_S_preds_male,
                                  G_hat = CV_G_preds_male,
                                  folds = male_folds,
                                  ss_folds = male_ss_folds,
                                  sample_split = sample_split,
                                  scale_est = TRUE)
    
    output_auc_male$vim <- "AUC"
    output_auc_male$SAB <- "male"
    output_auc_male <- output_auc_male %>% mutate(tau = landmark_time) %>%
      select(-landmark_time)
    output_cindex_male <- survML::vim_cindex(time = time_male,
                                        event = event_male,
                                        approx_times = approx_times,
                                        tau = max(approx_times),
                                        f_hat = CV_full_preds_cindex_male,
                                        fs_hat = CV_reduced_preds_cindex_male,
                                        S_hat = CV_S_preds_male,
                                        G_hat = CV_G_preds_male,
                                        folds = male_folds,
                                        ss_folds = male_ss_folds,
                                        sample_split = sample_split,
                                        scale_est = TRUE)
    output_cindex_male$vim <- "cindex"
    output_cindex_male$SAB <- "male"
    output_cindex_male <- output_cindex_male %>% mutate(tau = restriction_time) %>%
      select(-restriction_time)
    
    
    ### FEMALE
    CV_full_preds_landmark_female <- lapply(1:V,
                                          FUN = function(x){
                                            dat_j <- dat[folds == x,]
                                            females_j <- which(dat_j$female == 1)
                                            return(CV_full_preds_landmark[[x]][females_j,])
                                          })
    CV_reduced_preds_landmark_female <- lapply(1:V,
                                             FUN = function(x){
                                               dat_j <- dat[folds == x,]
                                               females_j <- which(dat_j$female == 1)
                                               return(CV_reduced_preds_landmark[[x]][females_j,])
                                             })
    CV_full_preds_cindex_female <- lapply(1:V,
                                        FUN = function(x){
                                          dat_j <- dat[folds == x,]
                                          females_j <- which(dat_j$female == 1)
                                          return(CV_full_preds_cindex[[x]][females_j])
                                        })
    CV_reduced_preds_cindex_female <- lapply(1:V,
                                           FUN = function(x){
                                             dat_j <- dat[folds == x,]
                                             females_j <- which(dat_j$female == 1)
                                             return(CV_reduced_preds_cindex[[x]][females_j])
                                           })
    CV_S_preds_female <- lapply(1:V,
                              FUN = function(x){
                                dat_j <- dat[folds == x,]
                                females_j <- which(dat_j$female == 1)
                                return(CV_S_preds[[x]][females_j,])
                              })
    CV_G_preds_female <- lapply(1:V,
                              FUN = function(x){
                                dat_j <- dat[folds == x,]
                                females_j <- which(dat_j$female == 1)
                                return(CV_G_preds[[x]][females_j,])
                              })
    
    output_auc_female <- survML::vim_AUC(time = time_female,
                                       event = event_female,
                                       approx_times = approx_times,
                                       landmark_times = landmark_times,
                                       f_hat = lapply(CV_full_preds_landmark_female, function(x) 1-x),
                                       fs_hat = lapply(CV_reduced_preds_landmark_female, function(x) 1-x),
                                       S_hat = CV_S_preds_female,
                                       G_hat = CV_G_preds_female,
                                       folds = female_folds,
                                       ss_folds = female_ss_folds,
                                       sample_split = sample_split,
                                       scale_est = TRUE)
    
    output_auc_female$vim <- "AUC"
    output_auc_female$SAB <- "female"
    output_auc_female <- output_auc_female %>% mutate(tau = landmark_time) %>%
      select(-landmark_time)
    output_cindex_female <- survML::vim_cindex(time = time_female,
                                             event = event_female,
                                             approx_times = approx_times,
                                             tau = max(approx_times),
                                             f_hat = CV_full_preds_cindex_female,
                                             fs_hat = CV_reduced_preds_cindex_female,
                                             S_hat = CV_S_preds_female,
                                             G_hat = CV_G_preds_female,
                                             folds = female_folds,
                                             ss_folds = female_ss_folds,
                                             sample_split = sample_split,
                                             scale_est = TRUE)
    output_cindex_female$vim <- "cindex"
    output_cindex_female$SAB <- "female"
    output_cindex_female <- output_cindex_female %>% mutate(tau = restriction_time) %>%
      select(-restriction_time)
    
    
    output <- rbind(output_auc_male, output_cindex_male, output_auc_female, output_cindex_female)
    output$indx <- rep(char_indx, nrow(output))
    output$indx_name <- rep(char_indx_name, nrow(output))
    if (!(i == 1)){
      pooled_output <- rbind(pooled_output, output)
    } else{
      pooled_output <- output
    }
  }
  
  dat <- pooled_output
  dat$approach <- approach
  dat$seed <- seed
  dat$global_seed <- global_seed
  return(dat)
}

