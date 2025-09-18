# R/03_cage_signal.R

#' Compute normalized CAGE signal (RPM) at +/- window around TSS
compute_cage_rpm <- function(cage_bw_file, tss_gr, chrom, window_size) {
  cage <- rtracklayer::import(cage_bw_file, format = "bigWig")
  cage <- cage[seqnames(cage) == chrom]
  win  <- promoters_around(tss_gr, upstream = window_size, downstream = window_size + 1)
  
  ov   <- GenomicRanges::findOverlaps(win, cage)
  ssum <- tapply(cage$score[BiocGenerics::subjectHits(ov)],
                 BiocGenerics::queryHits(ov), sum, default = 0)
  cage_signal <- numeric(length(win))
  cage_signal[as.integer(names(ssum))] <- ssum
  
  total <- sum(cage$score)
  cage_signal / (total / 1e6)
}

#' Select most expressed TSS per gene by CAGE
select_max_tss_per_gene <- function(tss_gr, cage_signal) {
  tss_by_gene  <- GenomicRanges::split(tss_gr, tss_gr$gene_id)
  sig_by_gene  <- split(cage_signal, tss_gr$gene_id)
  idx_max      <- sapply(sig_by_gene, which.max)
  picked_tss   <- mapply(function(gr, i) gr[i], tss_by_gene, idx_max, SIMPLIFY = FALSE)
  selected_tss <- unlist(GRangesList(picked_tss))
  selected_sig <- mapply(function(v, i) v[i], sig_by_gene, idx_max)
  list(tss = selected_tss, signal = as.numeric(selected_sig))
}
