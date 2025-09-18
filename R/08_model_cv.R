# R/08_model_cv.R

misclass_rate <- function(y, yhat) mean(y != yhat)
rmse_log      <- function(t, p) sqrt(mean((t - p)^2, na.rm = TRUE))
pearson_r     <- function(t, p) stats::cor(t, p, use = "complete.obs")

run_kfold <- function(X, Y, is_on, k = 10, seed = 123,
                      chosen_combo = "RF_Clf_RF_Reg", chosen_iteration = 1) {
  
  set.seed(seed)
  n <- nrow(X); folds <- sample(rep(1:k, length.out = n))
  auc_mat <- matrix(0, nrow = k, ncol = 3, dimnames = list(NULL, c("Logistic","RF_Clf","SVM_Clf")))
  mis_mat <- matrix(0, nrow = k, ncol = 3, dimnames = list(NULL, c("Logistic","RF_Clf","SVM_Clf")))
  
  combos <- c("Logistic_Lasso","Logistic_RF_Reg","Logistic_M","Logistic_Sv",
              "RF_Clf_Lasso","RF_Clf_RF_Reg","RF_Clf_M","RF_Clf_Sv",
              "SVM_Clf_Lasso","SVM_Clf_RF_Reg","SVM_Clf_M","SVM_Clf_Sv")
  rmse_mat <- matrix(0, nrow = k, ncol = length(combos), dimnames = list(NULL, combos))
  r_mat    <- matrix(0, nrow = k, ncol = length(combos), dimnames = list(NULL, combos))
  
  rf_roc_to_plot <- NULL; store_actual <- store_pred <- NULL
  
  for (i in 1:k) {
    te <- which(folds == i); tr <- setdiff(seq_len(n), te)
    Xtr <- X[tr, ]; Ytr <- Y[tr]; ONtr <- is_on[tr]
    Xte <- X[te, ]; Yte <- Y[te]; ONte <- is_on[te]
    
    dtr <- data.frame(Xtr, is_on = factor(ONtr, levels = c(FALSE, TRUE)))
    dte <- data.frame(Xte, is_on = factor(ONte, levels = c(FALSE, TRUE)))
    
    # classifiers
    log_m <- glm(is_on ~ ., data = dtr, family = "binomial")
    log_p <- predict(log_m, newdata = dte, type = "response")
    log_c <- log_p > 0.5
    
    rf_c  <- randomForest::randomForest(is_on ~ ., data = dtr)
    rf_p  <- predict(rf_c, newdata = dte, type = "prob")[, 2]
    rf_cy <- rf_p > 0.5
    rf_roc <- pROC::roc(ONte, rf_p); if (i == 1) rf_roc_to_plot <- rf_roc
    
    svm_c <- e1071::svm(is_on ~ ., data = dtr, probability = TRUE, kernel = "radial")
    svm_p <- attr(predict(svm_c, newdata = dte, probability = TRUE), "probabilities")[, 2]
    svm_cy <- svm_p > 0.5
    
    auc_mat[i, "Logistic"] <- pROC::auc(pROC::roc(ONte, log_p))
    auc_mat[i, "RF_Clf" ]  <- pROC::auc(rf_roc)
    auc_mat[i, "SVM_Clf"]  <- pROC::auc(pROC::roc(ONte, svm_p))
    mis_mat[i, ] <- c(misclass_rate(ONte, log_c), misclass_rate(ONte, rf_cy), misclass_rate(ONte, svm_cy))
    
    # regressors on ON genes only
    on_tr  <- which(ONtr); on_te <- which(ONte)
    Xtr_on <- Xtr[on_tr, ]; Ytr_on <- Ytr[on_tr]
    Xte_on <- Xte[on_te, ]; Yte_on <- Yte[on_te]
    
    las   <- glmnet::cv.glmnet(Xtr_on, Ytr_on, alpha = 1, family = "gaussian")
    y_las <- as.numeric(predict(las, newx = Xte_on, s = las$lambda.min))
    
    rf_r  <- randomForest::randomForest(x = Xtr_on, y = Ytr_on)
    y_rfr <- predict(rf_r, newdata = Xte_on)
    
    mar   <- earth::earth(Ytr_on ~ ., data = data.frame(Xtr_on, Ytr_on = Ytr_on), degree = 1)
    y_mar <- predict(mar, newdata = data.frame(Xte_on))
    
    svm_r <- e1071::svm(Ytr_on ~ ., data = data.frame(Xtr_on, Ytr_on = Ytr_on), kernel = "radial")
    y_svm <- predict(svm_r, newdata = data.frame(Xte_on))
    
    combo_list <- list(
      Logistic_Lasso  = list(on_idx = which(log_c), pred = y_las),
      Logistic_RF_Reg = list(on_idx = which(log_c), pred = y_rfr),
      Logistic_M      = list(on_idx = which(log_c), pred = y_mar),
      Logistic_Sv     = list(on_idx = which(log_c), pred = y_svm),
      
      RF_Clf_Lasso    = list(on_idx = which(rf_cy), pred = y_las),
      RF_Clf_RF_Reg   = list(on_idx = which(rf_cy), pred = y_rfr),
      RF_Clf_M        = list(on_idx = which(rf_cy), pred = y_mar),
      RF_Clf_Sv       = list(on_idx = which(rf_cy), pred = y_svm),
      
      SVM_Clf_Lasso   = list(on_idx = which(svm_cy), pred = y_las),
      SVM_Clf_RF_Reg  = list(on_idx = which(svm_cy), pred = y_rfr),
      SVM_Clf_M       = list(on_idx = which(svm_cy), pred = y_mar),
      SVM_Clf_Sv      = list(on_idx = which(svm_cy), pred = y_svm)
    )
    
    for (nm in names(combo_list)) {
      mod <- combo_list[[nm]]
      on_pred_on <- intersect(on_te, mod$on_idx)
      rel <- match(on_pred_on, on_te); rel <- rel[!is.na(rel)]
      if (length(rel) > 0) {
        rmse_mat[i, nm] <- rmse_log(Yte_on[rel], mod$pred[rel])
        r_mat[i, nm]    <- pearson_r(Yte_on[rel], mod$pred[rel])
        if (nm == chosen_combo && i == chosen_iteration) {
          store_actual <- Yte_on[rel]; store_pred <- mod$pred[rel]
        }
      } else {
        rmse_mat[i, nm] <- NA; r_mat[i, nm] <- NA
      }
    }
  }
  
  list(
    auc = colMeans(auc_mat),
    misclass = colMeans(mis_mat, na.rm = TRUE),
    rmse = colMeans(rmse_mat, na.rm = TRUE),
    r = colMeans(r_mat, na.rm = TRUE),
    rf_roc = rf_roc_to_plot,
    store_actual = store_actual,
    store_pred = store_pred
  )
}
