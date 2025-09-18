# R/00_utils_io.R

#' Load YAML config into a list
#' @param path path to YAML file
#' @return named list
read_yaml_list <- function(path) {
  yaml::read_yaml(path)
}

#' Safe saveRDS with auto dir create
save_rds <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, path)
  message("Saved: ", path)
}

#' Safe readRDS check
  if (!file.exists(path)) stop("File not found: ", path)
  readRDS(path)
}