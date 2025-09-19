#!/usr/bin/env Rscript

# Reads conf/dataspec.yml, discovers appropriate files from the web directories,
# downloads them into data/raw/, decompresses GTF if needed,
# and writes conf/paths.yml with the local paths.

suppressPackageStartupMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml")
  if (!requireNamespace("digest", quietly = TRUE)) install.packages("digest")
  if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")
})

# ---- helpers ----
ensure_dir <- function(path) dir.create(path, showWarnings = FALSE, recursive = TRUE)

read_text <- function(url) {
  con <- url(url); on.exit(close(con), add = TRUE)
  readLines(con, warn = FALSE)
}

extract_hrefs <- function(lines) {
  # very lightweight: extract href="..."
  href_lines <- grep('href="', lines, value = TRUE, ignore.case = TRUE)
  hrefs <- gsub('.*href="([^"]+)".*', '\\1', href_lines)
  unique(hrefs)
}

# keep only files (not ../) and optionally only given extensions
filter_files <- function(names, allowed_exts = NULL) {
  names <- names[!grepl("^\\?dir=", names)]                     # sometimes indexes add query params
  names <- names[!grepl("^\\.{1,2}/?$", names)]                 # drop . and ..
  names <- names[!grepl("/$", names)]                           # drop directories
  if (!is.null(allowed_exts) && length(allowed_exts) > 0) {
    pat <- paste0("(", paste(gsub("\\.", "\\\\.", allowed_exts), collapse = "|"), ")$")
    names <- names[grepl(pat, names, ignore.case = TRUE)]
  }
  names
}

all_patterns_match <- function(x, patterns) {
  if (is.null(patterns) || length(patterns) == 0) return(TRUE)
  all(vapply(patterns, function(p) grepl(p, x, ignore.case = TRUE), logical(1)))
}

any_pattern_match <- function(x, patterns) {
  if (is.null(patterns) || length(patterns) == 0) return(TRUE)
  any(vapply(patterns, function(p) grepl(p, x, ignore.case = TRUE), logical(1)))
}

download_if_needed <- function(url, dest) {
  ensure_dir(dirname(dest))
  if (!file.exists(dest)) {
    message("Downloading: ", url)
    utils::download.file(url, dest, mode = "wb", quiet = FALSE)
  } else {
    message("Exists: ", dest)
  }
  dest
}

gunzip_if_needed <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    out <- sub("\\.gz$", "", path, ignore.case = TRUE)
    if (!file.exists(out)) {
      message("Decompressing: ", basename(path), " -> ", basename(out))
      R.utils::gunzip(path, destname = out, overwrite = TRUE, remove = FALSE)
    }
    return(out)
  }
  path
}

write_paths_yml <- function(gtf_file, cage_file, histone_files) {
  y <- list(
    gtf_file = gtf_file,
    cage_bw_file = cage_file,
    histone_files = as.list(histone_files)
  )
  yaml::write_yaml(y, "conf/paths.yml")
  message("\nWrote conf/paths.yml with ", length(histone_files), " predictor files.")
}

# ---- main ----
spec <- yaml::read_yaml("conf/dataspec.yml")
allowed_exts <- spec$allowed_exts %||% c(".bw", ".bigWig")
`%||%` <- function(a, b) if (!is.null(a)) a else b

ensure_dir("data/raw")

# --- GTF ---
gtf_mode <- spec$gtf$mode
gtf_url <- NULL
if (gtf_mode == "explicit") {
  gtf_url <- spec$gtf$explicit_url
  if (is.null(gtf_url)) stop("gtf.explicit_url is NULL")
} else if (gtf_mode == "latest_for_release") {
  stop("gtf.mode 'latest_for_release' is not implemented in this script. Use 'explicit'.")
}
gtf_dest <- file.path("data/raw", basename(gtf_url))
download_if_needed(gtf_url, gtf_dest)
gtf_local <- gunzip_if_needed(gtf_dest)

# --- CAGE ---
cage_file_local <- NULL
if (spec$cage$mode == "explicit") {
  cu <- spec$cage$explicit_url
  if (is.null(cu)) stop("cage.explicit_url is NULL while mode='explicit'")
  cage_dest <- file.path("data/raw", basename(cu))
  download_if_needed(cu, cage_dest)
  cage_file_local <- cage_dest
} else {
  base <- spec$cage$base_url
  inc <- spec$cage$include_patterns
  exc <- spec$cage$exclude_patterns %||% character(0)
  idx <- read_text(base)
  hrefs <- extract_hrefs(idx)
  files <- filter_files(hrefs, allowed_exts)
  cand <- files[vapply(files, function(fn) all_patterns_match(fn, inc) && !any_pattern_match(fn, exc), logical(1))]
  if (length(cand) < 1) stop("No CAGE file matched include_patterns at: ", base)
  # pick the first match
  chosen <- cand[1]
  cage_url <- paste0(base, chosen)
  cage_dest <- file.path("data/raw", basename(chosen))
  download_if_needed(cage_url, cage_dest)
  cage_file_local <- cage_dest
}

# --- Predictors ---
pred_mode <- spec$predictors$mode
pred_base <- spec$predictors$base_url
histone_locals <- character(0)

if (pred_mode == "construct") {
  pfx  <- spec$predictors$construct$prefix
  cell <- spec$predictors$construct$cell_line
  marks <- spec$predictors$construct$marks
  sfx  <- spec$predictors$construct$suffix
  fnames <- paste0(pfx, cell, marks, sfx)
  urls   <- paste0(pred_base, fnames)
  for (i in seq_along(urls)) {
    dest <- file.path("data/raw", basename(urls[i]))
    try(download_if_needed(urls[i], dest), silent = FALSE)
    if (file.exists(dest)) histone_locals <- c(histone_locals, dest)
  }
} else if (pred_mode == "pattern") {
  inc <- spec$predictors$pattern$include_patterns
  exc <- spec$predictors$pattern$exclude_patterns %||% character(0)
  idx <- read_text(pred_base)
  hrefs <- extract_hrefs(idx)
  files <- filter_files(hrefs, allowed_exts)
  cand <- files[vapply(files, function(fn) all_patterns_match(fn, inc) && !any_pattern_match(fn, exc), logical(1))]
  if (length(cand) < 1) stop("No predictor files matched include_patterns at: ", pred_base)
  for (fn in cand) {
    url  <- paste0(pred_base, fn)
    dest <- file.path("data/raw", basename(fn))
    try(download_if_needed(url, dest), silent = FALSE)
    if (file.exists(dest)) histone_locals <- c(histone_locals, dest)
  }
} else if (pred_mode == "explicit_list") {
  urls <- spec$predictors$explicit_urls
  if (length(urls) < 1) stop("predictors.explicit_urls is empty while mode='explicit_list'")
  for (u in urls) {
    dest <- file.path("data/raw", basename(u))
    try(download_if_needed(u, dest), silent = FALSE)
    if (file.exists(dest)) histone_locals <- c(histone_locals, dest)
  }
} else {
  stop("Unknown predictors.mode: ", pred_mode)
}

if (length(histone_locals) < 1) stop("No predictor files were downloaded or found.")

# --- Write conf/paths.yml for the pipeline ---
write_paths_yml(gtf_local, cage_file_local, histone_locals)

message("\nDiscovery & download completed.\n- GTF: ", gtf_local, "\n- CAGE: ", cage_file_local,
        "\n- Predictors (", length(histone_locals), "): first => ", histone_locals[1])
