# R/09_var_importance.R

rf_importances_multi <- function(X, Y, is_on, seeds = c(111,222,333,444,555,666,777,888,999,1010)) {
  clf_list <- list(); reg_list <- list()
  on_idx <- which(is_on); X_on <- X[on_idx, ]; Y_on <- Y[on_idx]
  
  for (sd in seeds) {
    set.seed(sd)
    d <- data.frame(X, is_on = factor(is_on, levels = c(FALSE, TRUE)))
    rf_c <- randomForest::randomForest(is_on ~ ., data = d)
    clf_list[[length(clf_list) + 1]] <- randomForest::importance(rf_c)[, "MeanDecreaseGini"]
    
    rf_r <- randomForest::randomForest(x = X_on, y = Y_on)
    reg_list[[length(reg_list) + 1]] <- randomForest::importance(rf_r)[, "IncNodePurity"]
  }
  
  clf_mat <- do.call(cbind, clf_list); reg_mat <- do.call(cbind, reg_list)
  list(
    clf_df = data.frame(Feature = rownames(clf_mat), MeanDecreaseGini = rowMeans(clf_mat, na.rm = TRUE)),
    reg_df = data.frame(Feature = rownames(reg_mat), IncNodePurity   = rowMeans(reg_mat, na.rm = TRUE)),
    clf_mat = clf_mat, reg_mat = reg_mat
  )
}
