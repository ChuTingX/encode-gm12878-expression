# scripts/00_bootstrap_env.R
# One-time environment setup with renv + Bioconductor.

message(">>> Bootstrapping renv + Bioconductor environment...")

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!file.exists("renv.lock")) {
  renv::init(bare = TRUE)
} else {
  renv::restore()
}

# Core packages (Bioc first, then CRAN)
bioc <- c("rtracklayer","GenomicRanges","BiocParallel")
cran <- c("glmnet","randomForest","earth","pROC","e1071","ggplot2",
          "yaml","optparse","targets","tarchetypes","testthat","rmarkdown",
          "digest","R.utils")

BiocManager::install(c(bioc, cran), ask = FALSE, update = TRUE)
renv::snapshot(prompt = FALSE)

message(">>> Environment ready. Use renv::restore() on new machines.")
