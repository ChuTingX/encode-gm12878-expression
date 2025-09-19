# scripts/render_report.R
if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown")
dir.create("reports", showWarnings = FALSE, recursive = TRUE)

cand <- c("reports/analysis.Rmd", "reports/paper.Rmd")
for (f in cand) {
  if (file.exists(f)) {
    message("Rendering: ", f)
    rmarkdown::render(f, output_format = "html_document", output_dir = "reports")
  } else {
    message("Skipping (not found): ", f)
  }
}
