# R/02_gtf_tss.R

#' Import transcripts on a chromosome and length-filter
import_and_filter_transcripts <- function(gtf_file, chrom, min_length) {
  g <- rtracklayer::import(gtf_file, format = "gtf")
  g <- g[seqnames(g) == chrom]
  tx <- g[g$type == "transcript"]
  tx[width(tx) >= min_length]
}

#' Derive 1-bp TSS GRanges with gene_id
make_tss <- function(transcripts) {
  tss_pos <- ifelse(as.character(GenomicRanges::strand(transcripts)) == "+",
                    GenomicRanges::start(transcripts),
                    GenomicRanges::end(transcripts))
  GRanges::GRanges(
    seqnames = GenomicRanges::seqnames(transcripts),
    ranges   = IRanges::IRanges(tss_pos, width = 1),
    strand   = GenomicRanges::strand(transcripts),
    gene_id  = transcripts$gene_id
  ) |> unique()
}

#' Promoter window around TSS (strand-aware)
promoters_around <- function(x, upstream = 2000, downstream = 200) {
  out <- x
  plus  <- which(GenomicRanges::strand(x) == "+")
  minus <- which(GenomicRanges::strand(x) == "-")
  GenomicRanges::start(out)[plus]  <- GenomicRanges::start(x)[plus] - upstream
  GenomicRanges::end(out)[plus]    <- GenomicRanges::start(x)[plus] + downstream - 1
  GenomicRanges::start(out)[minus] <- GenomicRanges::end(x)[minus] - downstream + 1
  GenomicRanges::end(out)[minus]   <- GenomicRanges::end(x)[minus] + upstream
  GenomicRanges::trim(out)
}
