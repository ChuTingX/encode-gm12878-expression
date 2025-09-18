# R/07_features.R

#' Build feature matrix X using best bin per track
make_feature_matrix <- function(histone_files, bins_list, pc_and_bin, hist_cov_list) {
  L <- lapply(seq_along(histone_files), function(i) {
    cov  <- hist_cov_list[[i]]
    mats <- lapply(bins_list, mean_bin_signal, cov = cov)
    M    <- do.call(rbind, mats)
    pc   <- pc_and_bin[[ histone_files[i] ]]$pseudocount
    bb   <- pc_and_bin[[ histone_files[i] ]]$best_bin
    log2(M + pc)[ , bb, drop = FALSE]
  })
  X <- do.call(cbind, L); colnames(X) <- histone_files; X
}

#' Build response Y (log2 RPM) and ON mask
make_labels <- function(cage_signal, thresh_log2 = 0.1) {
  Y <- rep(NA_real_, length(cage_signal))
  nz <- which(cage_signal > 0)
  Y[nz] <- log2(cage_signal[nz])
  Y[is.na(Y)] <- 0
  list(Y = Y, is_on = Y > thresh_log2)
}
