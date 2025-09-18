# R/06_best_bin.R

#' Pick best bin with fixed pseudocount by maximizing |cor| to expression
find_best_bin_fixed_pc <- function(hist_cov, bins_list, expression, pc = 0.1) {
  mats <- lapply(bins_list, mean_bin_signal, cov = hist_cov)
  if (length(mats) == 0) return(list(pseudocount = pc, best_bin = 1))
  M <- do.call(rbind, mats)
  if (is.null(M) || nrow(M) == 0 || all(is.na(M))) return(list(pseudocount = pc, best_bin = 1))
  cors <- apply(log2(M + pc), 2, function(x) stats::cor(x, expression, use = "pairwise.complete.obs"))
  if (all(is.na(cors))) return(list(pseudocount = pc, best_bin = 1))
  i <- which.max(abs(cors)); if (length(i) == 0 || is.na(i)) i <- 1
  list(pseudocount = pc, best_bin = i)
}

#' Loop over all histone files to get best bin specs
best_bins_all <- function(hist_cov_list, histone_files, D1_bins_list, D1_signal, pc) {
  out <- vector("list", length(histone_files)); names(out) <- histone_files
  for (i in seq_along(histone_files)) {
    out[[i]] <- find_best_bin_fixed_pc(hist_cov_list[[i]], D1_bins_list, D1_signal, pc = pc)
  }
  out
}
