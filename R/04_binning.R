# R/04_binning.R

#' Create bins around one TSS
create_bins_one <- function(tss, upstream, downstream, nbins_up, nbins_down) {
  s <- as.character(GenomicRanges::strand(tss))
  pos <- GenomicRanges::start(tss)
  if (s == "-") {
    region_start <- pos - downstream; region_end <- pos + upstream
  } else {
    region_start <- pos - upstream;   region_end <- pos + downstream
  }
  region <- GRanges::GRanges(GenomicRanges::seqnames(tss),
                             IRanges::IRanges(min(region_start, region_end), max(region_start, region_end)),
                             strand = s)
  total_bins <- nbins_up + nbins_down + 1
  brk <- seq(GenomicRanges::start(region), GenomicRanges::end(region), length.out = total_bins + 1)
  bins <- GRanges::GRanges(
    GenomicRanges::seqnames(region),
    IRanges::IRanges(start = brk[-(total_bins + 1)], end = brk[-1]),
    strand = s
  )
  GenomicRanges::trim(bins)
}

#' Vectorized binning for many TSS
make_all_bins <- function(selected_tss, upstream, downstream, nbins_up, nbins_down) {
  lapply(seq_along(selected_tss), \(i)
         create_bins_one(selected_tss[i], upstream, downstream, nbins_up, nbins_down))
}