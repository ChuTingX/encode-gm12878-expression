# R/99_sessioninfo.R

write_session_info <- function(path = "data/processed/sessionInfo.txt") {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  capture.output(sessionInfo(), file = path)
  message("Session info saved to: ", path)
}
