boost_c_index <- function(time, # follow up times
                          event, # event indicators
                          X, # feature matrix
                          S_hat, # conditional survival function predictions
                          V = 5, # cross validation fold number
                          approx_times, # times for approximating integrals
                          indx = NULL, # index of feature to be removed
                          tuning = "CV", # whether to do CV or simply fit a model
                          produce_fit = FALSE, # whether to produce a fit after CV
                          subsample_n = 500,
                          params = list(mstop = c(1000), # tuning parameter list
                                        nu = c(0.3),
                                        sigma = c(0.1),
                                        learner = "glm")){
  
  if (subsample_n < length(time)){
    subsample_inds <- sample(1:length(time), size = subsample_n, replace = FALSE)
    time <- time[subsample_inds]
    event <- event[subsample_inds]
    X <- X[subsample_inds,,drop=FALSE]
    S_hat <- S_hat[subsample_inds,,drop=FALSE]
  }
  
  print(length(time))
  print(length(event))
  print(dim(X))
  print(dim(S_hat))
  
  reverse_sorted_mstop <- params$mstop[order(params$mstop, decreasing = TRUE)]
  tau <- max(approx_times)
  
  k <- length(approx_times)
  n <- length(time)
  folds <- sample(rep(seq_len(V), length = n))
  
  param_grid <- expand.grid(mstop = max(params$mstop),
                            nu = params$nu,
                            sigma = params$sigma,
                            learner = params$learner)
  param_grid_full <- expand.grid(mstop = reverse_sorted_mstop,
                                 nu = params$nu,
                                 sigma = params$sigma,
                                 learner = params$learner)
  param_grid_full$CV_risk <- NA
  mod_list <- vector("list", nrow(param_grid))
  
  K <- length(params$mstop)
  for (i in 1:nrow(param_grid)){
    
    mstop_curr = param_grid[i,1]
    nu_curr <- param_grid[i,2]
    sigma_curr <- param_grid[i,3]
    learner_curr <- param_grid[i,4]
    
    Sweights_train <- function(i1, i2){
      # prob that i2 exceeds i1
      exceed_prob <- -sum(S_hat_train[i2,2:k] * diff(S_hat_train[i1,1:k]))
      return(exceed_prob)
    }
    
    # the i,j entry of this matrix is P(i fails before j)
    my_Cindex <- function (sigma = 0.1) {
      approxGrad <- function(x) { ## sigmoid function for gradient
        exp(x/sigma) / (sigma * (1 + exp(x/sigma))^2)
      }
      approxLoss <- function(x) { ## sigmoid function for loss
        1 / (1 + exp(x / sigma))
      }
      ngradient = function(y, f, w = 1) { ## negative gradient
        if (!all(w %in% c(0,1)))
          stop(sQuote("weights"), " must be either 0 or 1 for family ",
               sQuote("UnoC"))
        survtime <- y[,1]
        event <- y[,2]
        if (length(w) == 1) w <- rep(1, length(event))
        if (length(f) == 1) {
          f <- rep(f, length(survtime))
        }
        n <- length(survtime)
        etaj <- matrix(f, nrow = n, ncol = n, byrow = TRUE)
        etak <- matrix(f, nrow = n, ncol = n)
        etaMat <- etak - etaj
        rm(etaj); rm(etak);
        weights_out <- wweights#compute_weights(y, w)
        M1 <- approxGrad(etaMat) * weights_out
        ng <- colSums(M1) - rowSums(M1)
        return(ng)
      }
      risk = function(y, f, w = 1) { ## empirical risk
        survtime <- y[,1]
        event <- y[,2]
        if (length(f) == 1) {
          f <- rep(f, length(y))
        }
        n <- length(survtime)
        etaj <- matrix(f, nrow = n, ncol = n, byrow = TRUE)
        etak <- matrix(f, nrow = n, ncol = n)
        etaMat <- (etak - etaj)
        rm(etaj); rm(etak);
        weights_out <- wweights#compute_weights(y, w)
        M1 <- approxLoss(etaMat) * weights_out
        return(- sum(M1))
      }
      mboost::Family( ## build the family object
        ngradient = ngradient,
        risk = risk,
        weights = "zeroone",
        offset = function(y, w = 1) {0},
        check_y = function(y) {
          if (!inherits(y, "Surv"))
            stop("response is not an object of class ", sQuote("Surv"),
                 " but ", sQuote("family = UnoC()"))
          y},
        rclass = function(f){},
        name = paste("Efron-type c-index boosting")
      )
    }
    
    if (tuning == "CV"){
      CV_risks <- matrix(NA, nrow = V, ncol = K)
      for (j in 1:V){
        print(paste0("CV fold ", j))
        
        time_train <- time[folds != j]
        event_train <- event[folds != j]
        X_train <- X[folds != j,,drop=FALSE]
        S_hat_train <- S_hat[folds != j,]
        X_test <- X[folds == j,,drop=FALSE]
        S_hat_test <- S_hat[folds == j,]
        time_test <- time[folds == j]
        event_test <- event[folds == j]
        
        w <- rep(1, length(time_train))
        index_grid <- gtools::combinations(n = length(time_train),
                                           r = 2,
                                           v = 1:length(time_train),
                                           repeats.allowed = TRUE)
        weight_vec <- mapply(FUN = Sweights_train, index_grid[,1], index_grid[,2])
        wweights <- matrix(NA, length(time_train), length(time_train))
        wweights[lower.tri(wweights, diag=TRUE)] <- weight_vec
        wweights <- t(wweights)
        wweights[lower.tri(wweights)] <- 1-t(wweights)[lower.tri(wweights)]
        Wmat <- w %o% w
        wweights <- wweights * Wmat
        diag(wweights) <- 0.5
        wweights <- wweights / sum(wweights)
        
        dtrain <- data.frame(time = time_train,
                             event = event_train,
                             X_train)
        dtest <- data.frame(X_test)
        
        feature_names <- names(X_train)
        
        feature_names <- apply(as.matrix(feature_names),
                               MARGIN = 1,
                               FUN = function(x) paste0("", x, ""))
        feature_form <- paste(feature_names, collapse = "+")
        formula_text <- paste0("Surv(time, event) ~ ", feature_form)
        if (learner_curr == "tree"){
          mod <- mboost::blackboost(Surv(time, event) ~ .,
                                    family = my_Cindex(sigma = sigma_curr),
                                    control = mboost::boost_control(mstop = mstop_curr, 
                                                                    trace = TRUE, 
                                                                    nu = nu_curr),
                                    tree_controls = partykit::ctree_control(maxdepth = 2),
                                    data = dtrain)
          print(mod)
        } else if (learner_curr == "glm"){
          mod <- mboost::glmboost(Surv(time, event) ~ .^2,
                                  family = my_Cindex(sigma = sigma_curr),
                                  control = mboost::boost_control(mstop = mstop_curr,
                                                                  trace = TRUE,
                                                                  nu = nu_curr),
                                  data = dtrain)
        } else if (learner_curr == "gam"){
          mod <- mboost::gamboost(Surv(time, event) ~ .,
                                  family = my_Cindex(sigma = sigma_curr),
                                  control = mboost::boost_control(mstop = mstop_curr, trace = FALSE, nu = nu_curr),
                                  data = dtrain)
        }
        
        for (l in 1:K){
          sub_mstop = reverse_sorted_mstop[l]
          mboost::mstop(mod) <- sub_mstop
          preds <- -predict(mod, newdata = dtest)[,1]
          
          risk_j <- -survML:::estimate_cindex(time = time_test,
                                              event = event_test,
                                              approx_times = approx_times,
                                              preds = preds,
                                              S_hat = S_hat_test,
                                              G_hat = matrix(0.5, nrow = nrow(S_hat_test),
                                                             ncol = ncol(S_hat_test)),
                                              tau = tau)$plug_in
          CV_risks[j,l] <- risk_j
        }
      }
      param_grid_full$CV_risk[(((i-1)*K)+1):(i*K)] <- colMeans(CV_risks)
      
    } else if (tuning == "none"){
      
      time_train <- time
      event_train <- event
      X_train <- X
      S_hat_train <- S_hat
      dtrain <- data.frame(time = time_train,
                           event = event_train,
                           X_train)
      
      #Survival function version
      w <- rep(1, length(time_train))
      index_grid <- gtools::combinations(n = length(time_train),
                                         r = 2,
                                         v = 1:length(time_train),
                                         repeats.allowed = TRUE)
      weight_vec <- mapply(FUN = Sweights_train, index_grid[,1], index_grid[,2])
      wweights <- matrix(NA, length(time_train), length(time_train))
      wweights[lower.tri(wweights, diag=TRUE)] <- weight_vec
      wweights <- t(wweights)
      wweights[lower.tri(wweights)] <- 1-t(wweights)[lower.tri(wweights)]
      Wmat <- w %o% w
      wweights <- wweights * Wmat
      diag(wweights) <- 0.5 # diagonal should be 1/2, I think?
      wweights <- wweights / sum(wweights)
      
      feature_names <- names(X_train)
      feature_names <- apply(as.matrix(feature_names),
                             MARGIN = 1,
                             FUN = function(x) paste0("", x, ""))
      feature_form <- paste(feature_names, collapse = "+")
      formula_text <- paste0("Surv(time, event) ~ ", feature_form)
      
      if (learner_curr == "tree"){
        mod <- mboost::blackboost(Surv(time, event) ~ .,
                                  family = my_Cindex(sigma = sigma_curr),
                                  control = mboost::boost_control(mstop = mstop_curr, trace = TRUE, nu = nu_curr),
                                  tree_controls = partykit::ctree_control(maxdepth = 2),
                                  data = dtrain)
      } else if (learner_curr == "glm"){
        mod <- mboost::glmboost(Surv(time, event) ~ .^2,
                                family = my_Cindex(sigma = sigma_curr),
                                control = mboost::boost_control(mstop = mstop_curr, trace = TRUE, nu = nu_curr),
                                data = dtrain)
      } else if (learner_curr == "gam"){
        mod <- mboost::gamboost(Surv(time, event) ~ .,
                                family = my_Cindex(sigma = sigma_curr),
                                control = mboost::boost_control(mstop = mstop_curr, trace = FALSE, nu = nu_curr),
                                data = dtrain)
      }
      
      mod_list[[i]] <- mod
      param_grid$CV_risk[i] <- 999
    }
  }
  
  if (tuning == "none"){
    param_grid_full <- param_grid
  }
  
  param_grid_full <- param_grid_full[order(param_grid_full$mstop), ]
  opt_index <- which.min(round(param_grid_full$CV_risk, digits = 3))
  
  if (tuning == "CV" & produce_fit){
    mstop_curr = param_grid_full[opt_index,1]
    nu_curr <- param_grid_full[opt_index,2]
    sigma_curr <- param_grid_full[opt_index,3]
    learner_curr <- param_grid_full[opt_index,4]
    time_train <- time
    event_train <- event
    X_train <- X
    S_hat_train <- S_hat
    dtrain <- data.frame(time = time_train,
                         event = event_train,
                         X_train)
    
    #Survival function version
    w <- rep(1, length(time_train))
    index_grid <- gtools::combinations(n = length(time_train),
                                       r = 2,
                                       v = 1:length(time_train),
                                       repeats.allowed = TRUE)
    weight_vec <- mapply(FUN = Sweights_train, index_grid[,1], index_grid[,2])
    wweights <- matrix(NA, length(time_train), length(time_train))
    wweights[lower.tri(wweights, diag=TRUE)] <- weight_vec
    wweights <- t(wweights)
    wweights[lower.tri(wweights)] <- 1-t(wweights)[lower.tri(wweights)]
    Wmat <- w %o% w
    wweights <- wweights * Wmat
    diag(wweights) <- 0.5 # diagonal should be 1/2, I think?
    wweights <- wweights / sum(wweights)
    
    feature_names <- names(X_train)
    feature_names <- apply(as.matrix(feature_names),
                           MARGIN = 1,
                           FUN = function(x) paste0("btree(", x, ")"))
    feature_form <- paste(feature_names, collapse = "+")
    formula_text <- paste0("Surv(time, event) ~ ", feature_form)
    
    if (learner_curr == "tree"){
      mod <- mboost::blackboost(Surv(time, event) ~ .,
                                family = my_Cindex(sigma = sigma_curr),
                                control = mboost::boost_control(mstop = mstop_curr, trace = TRUE, nu = nu_curr),
                                tree_controls = partykit::ctree_control(maxdepth = 2),
                                data = dtrain)
    } else if (learner_curr == "glm"){
      mod <- mboost::glmboost(Surv(time, event) ~ .^2,
                              family = my_Cindex(sigma = sigma_curr),
                              control = mboost::boost_control(mstop = mstop_curr, trace = TRUE, nu = nu_curr),
                              data = dtrain)
    } else if (learner_curr == "gam"){
      mod <- mboost::gamboost(Surv(time, event) ~ .,
                              family = my_Cindex(sigma = sigma_curr),
                              control = mboost::boost_control(mstop = mstop_curr, trace = FALSE, nu = nu_curr),
                              data = dtrain)
    }
    opt_model <- mod
  } else if (tuning == "none" & produce_fit){
    opt_model <- mod_list[[opt_index]]
  } else{
    opt_model <- NULL
  }
  return(list(param_grid = param_grid_full,
              opt_index = opt_index,
              opt_model = opt_model))
}