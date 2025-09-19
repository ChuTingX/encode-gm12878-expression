#!/usr/bin/env Rscript
# scripts/run_all.R

suppressPackageStartupMessages({
  library(optparse); library(yaml)
  library(rtracklayer); library(GenomicRanges); library(BiocParallel)
})

# ---- CLI ----
option_list <- list(
  make_option(c("-p","--paths"),  type="character", default="conf/paths.yml"),
  make_option(c("-a","--params"), type="character", default="conf/params.yml")
)
opt <- parse_args(OptionParser(option_list = option_list))
paths  <- yaml::read_yaml(opt$paths)
params <- yaml::read_yaml(opt$params)

stop_if_empty <- function(x, msg) if (is.null(x) || (is.list(x) && length(x) == 0) || (is.character(x) && length(x) == 0)) stop(msg)

stop_if_empty(paths$gtf_file,  "paths.yml: gtf_file is not set. Run scripts/discover_and_download.R")
stop_if_empty(paths$cage_bw_file, "paths.yml: cage_bw_file is not set. Run scripts/discover_and_download.R")
stop_if_empty(paths$histone_files, "paths.yml: histone_files is empty. Run scripts/discover_and_download.R")

# ---- source modules ----
src <- function(f) source(file.path("R", f))
src("00_utils_io.R"); src("01_setup_parallel.R")
src("02_gtf_tss.R");  src("03_cage_signal.R")
src("04_binning.R");  src("05_splits_and_cov.R")
src("06_best_bin.R"); src("07_features.R")
src("08_model_cv.R"); src("09_var_importance.R")
src("99_sessioninfo.R")

# ---- run ----
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
register_workers(params$num_workers)

message("1) Import GTF, filter transcripts, build TSS...")
tx  <- import_and_filter_transcripts(paths$gtf_file, params$chrom, params$min_length)
tss <- make_tss(tx)

message("2) Compute CAGE RPM, select per-gene max TSS...")
cage_rpm <- compute_cage_rpm(paths$cage_bw_file, tss, params$chrom, params$window_size)
sel <- select_max_tss_per_gene(tss, cage_rpm)
selected_tss <- sel$tss; selected_signal <- sel$signal
save_rds(list(tss=selected_tss, signal=selected_signal), "data/processed/selected_tss_signal.rds")

message("3) Create bins around TSS...")
bins <- make_all_bins(selected_tss, params$upstream, params$downstream, params$nbins_up, params$nbins_down)
save_rds(bins, "data/processed/bins.rds")

message("4) Split genes into D1/D2...")
sp <- split_genes(length(selected_tss), prop_train = 1/3, seed = params$seed)
D1_bins <- bins[sp$D1]; D2_bins <- bins[sp$D2]
D1_sig  <- selected_signal[sp$D1]; D2_sig <- selected_signal[sp$D2]

message("5) Compute coverage on D1 and choose best bins (fixed pc)...")
cov_d1 <- coverage_for_files(paths$histone_files, params$chrom)
pcbin  <- best_bins_all(cov_d1, paths$histone_files, D1_bins, D1_sig, pc = params$fixed_pseudocount)
save_rds(pcbin, "data/processed/best_bins.rds")

message("6) Coverage on D2, build features + labels...")
cov_d2 <- coverage_for_files(paths$histone_files, params$chrom)
X      <- make_feature_matrix(paths$histone_files, D2_bins, pcbin, cov_d2)
lab    <- make_labels(D2_sig, thresh_log2 = 0.1)
Y <- lab$Y; is_on <- lab$is_on
save_rds(list(X=X,Y=Y,is_on=is_on), "data/processed/dataset_D2.rds")

message("7) 10-fold CV (classifiers + regressors, two-step combos)...")
cvres <- run_kfold(X, Y, is_on, k = params$kfolds,
                   seed = params$seed,
                   chosen_combo = params$chosen_combo,
                   chosen_iteration = params$chosen_iteration)
save_rds(cvres, "data/processed/cv_results.rds")

message("8) Variable importance across seeds...")
imps <- rf_importances_multi(X, Y, is_on)
save_rds(imps, "data/processed/importances.rds")

message("9) Write session info...")
write_session_info()

cat("\n=== Summary (means over folds) ===\n")
cat("Misclass:\n"); print(cvres$misclass)
cat("AUC:\n");      print(cvres$auc)
cat("RMSE:\n");     print(cvres$rmse)
cat("Pearson r:\n");print(cvres$r)
cat("\nDone. Cached outputs in data/processed/ .\n")
