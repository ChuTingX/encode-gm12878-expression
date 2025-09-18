# R/01_setup_parallel.R

#' Register BiocParallel workers cross-platform
#' @param workers integer
register_workers <- function(workers = parallel::detectCores()) {
  if (.Platform$OS.type == "windows") {
    BiocParallel::register(BiocParallel::SnowParam(workers))
  } else {
    BiocParallel::register(BiocParallel::MulticoreParam(workers))
  }
  invisible(TRUE)
}