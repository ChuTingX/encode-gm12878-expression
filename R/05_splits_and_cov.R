# R/05_splits_and_cov.R

#' Train/test gene split at gene level
split_genes <- function(n_genes, prop_train = 1/3, seed = 123) {
  set.seed(seed)
  D1_size <- round(n_genes * prop_train)
  D1_idx  <- sample(seq_len(n_genes), size = D1_size)
  list(D1 = D1_idx, D2 = setdiff(seq_len(n_genes), D1_idx))
}

#' Mean coverage in bins given a coverage RleList
mean_bin_signal <- function(bins, cov) {
  bins_chr <- split(bins, GenomicRanges::seqnames(bins))
  res <- numeric(length(bins)); start_idx <- 1
  for (chr in names(bins_chr)) {
    chr_bins <- bins_chr[[chr]]
    if (!chr %in% names(cov)) {
      end_idx <- start_idx + length(chr_bins) - 1
      res[start_idx:end_idx] <- NA; start_idx <- end_idx + 1; next
    }
    v   <- IRanges::Views(cov[[chr]], start = GenomicRanges::start(chr_bins),
                          end   = GenomicRanges::end(chr_bins))
    m   <- IRanges::viewMeans(v)
    end_idx <- start_idx + length(chr_bins) - 1
    res[start_idx:end_idx] <- if (length(m) == 0) rep(NA, length(chr_bins)) else m
    start_idx <- end_idx + 1
  }
  res
}

#' Import bigWigs and compute coverage per file (filtered to chrom)
coverage_for_files <- function(histone_files, chrom) {
  BiocParallel::bplapply(histone_files, function(hf) {
    x <- rtracklayer::import(hf, format = "bigWig")
    x <- x[seqnames(x) == chrom]
    GenomicRanges::coverage(x, weight = "score")
  })
}
