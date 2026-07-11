library(testthat)
library(tidyverse)
library(patchwork)
library(multcompView)
library(ggsignif)
library(cluster)

# Load only function definitions from 06_clustering.R (combine_panels,
# run_clustering, run_validation, compute_roc), stopping before the driver
# code that loads real ~250k-pixel data. Real-data loads at the top
# (sn_boundary, index_config) are skipped -- these tests use their own
# synthetic sn_boundary/cluster_cols/clust_labels instead.
exprs    <- parse("../src/06_clustering.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl("st_read\\(|^results <-|^n_idx|^p_sweep|^ggsave|^p_map|^p_scatter|^p_ts|^val_results|^p_val|^roc_results|^roc_combined|^roc_colors|^auc_label|^p_roc|^for \\(idx in INDICES\\)", src_text)
stop_at <- which(grepl("^for \\(idx in INDICES\\)", src_text))[1]
stopifnot(!is.na(stop_at))
for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]], envir = .GlobalEnv)
}

# These are normally set from the real data load / top-level params;
# defined here directly since the real loads above were skipped.
K <- 3
cluster_cols <- setNames(c("#2166AC", "#92C5DE", "#D73027"), c("1", "2", "3"))
clust_labels <- c("1 - no CSD", "2 - partial", "3 - full CSD")

# ----------------------------------------------------------------------------
# combine_panels() -- direct regression guard for the bug this session
# found: chart code used to hardcode results$ndvi/results$ndwi, crashing
# whenever INDICES had only one entry.
# ----------------------------------------------------------------------------
test_that("combine_panels returns a single panel unchanged when given one", {
  p <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) + geom_point()
  result <- combine_panels(list(p))
  expect_identical(result, p)
})

test_that("combine_panels combines two panels into a patchwork object", {
  p1 <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) + geom_point()
  p2 <- ggplot(data.frame(x = 2, y = 2), aes(x, y)) + geom_point()
  result <- combine_panels(list(p1, p2))
  expect_s3_class(result, "patchwork")
})

test_that("combine_panels with combine_op = `/` stacks rather than errors", {
  p1 <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) + geom_point()
  p2 <- ggplot(data.frame(x = 2, y = 2), aes(x, y)) + geom_point()
  expect_no_error(combine_panels(list(p1, p2), combine_op = `/`))
})

# ----------------------------------------------------------------------------
# run_clustering() -- synthetic tau data with 3 obviously-separated groups;
# verifies stress-rank relabeling puts the lowest-composite group in
# cluster 1 and the highest in cluster K, regardless of kmeans' raw
# (arbitrary) cluster numbering.
# ----------------------------------------------------------------------------
test_that("run_clustering relabels clusters by ascending composite stress", {
  set.seed(1)
  n_per_group <- 20
  d <- data.frame(
    pixel_id = 1:60,
    lon = -120, lat = 38,
    tau_r1   = c(rnorm(n_per_group, -0.8, 0.05), rnorm(n_per_group, 0, 0.05), rnorm(n_per_group, 0.8, 0.05)),
    tau_var  = c(rnorm(n_per_group, -0.8, 0.05), rnorm(n_per_group, 0, 0.05), rnorm(n_per_group, 0.8, 0.05)),
    tau_skew = rnorm(60, 0, 0.1),
    Actual_Drop_15_16 = rnorm(60, 0.1, 0.02)
  )
  tau_path <- tempfile(fileext = ".csv")
  write_csv(d, tau_path)

  cfg <- list(tau_path = tau_path, drop_col = "Actual_Drop_15_16",
             clusters_out = tempfile(fileext = ".csv"), label = "TEST")

  result <- run_clustering(cfg, K = 3, SEED = 1)

  # The 20 lowest-tau_r1/tau_var pixels (rows 1:20) must all land in
  # cluster "1" (least stressed); the 20 highest (rows 41:60) in cluster
  # "3" (most stressed) -- regardless of kmeans' internal cluster numbering.
  low_group_clusters  <- result$d$cluster[result$d$pixel_id %in% 1:20]
  high_group_clusters <- result$d$cluster[result$d$pixel_id %in% 41:60]
  expect_true(all(low_group_clusters == "1"))
  expect_true(all(high_group_clusters == "3"))
})

# ----------------------------------------------------------------------------
# compute_roc() -- synthetic data with a hand-computable confusion matrix.
# ----------------------------------------------------------------------------
test_that("compute_roc's confusion matrix and AUC match hand-calculated values", {
  # 10 pixels: top 20% (2 pixels) by Actual_Drop_15_16 = "severe".
  # Cluster 3 membership perfectly predicts severity here by construction.
  d <- data.frame(
    pixel_id = 1:10,
    tau_r1  = c(5, 4, 3, 2, 1, -1, -2, -3, -4, -5),
    tau_var = 0,
    cluster = factor(c(3, 3, 2, 2, 2, 1, 1, 1, 1, 1), levels = 1:3),
    Actual_Drop_15_16 = c(10, 9, 3, 3, 3, 1, 1, 1, 1, 1)
  )
  cfg <- list(drop_col = "Actual_Drop_15_16", label = "TEST")
  res <- list(d = d, cfg = cfg)

  result <- compute_roc(res, K = 3)

  # Top 20% of 10 = top 2 rows (Actual_Drop_15_16 = 10, 9) -> both cluster 3
  # -> perfect separation: TP=2, FP=0, TN=8, FN=0.
  expect_equal(result$auc, 1, tolerance = 1e-6)
})

test_that("compute_roc gives AUC near 0.5 when the predictor carries no signal", {
  set.seed(2)
  n <- 200
  d <- data.frame(
    pixel_id = 1:n,
    tau_r1  = rnorm(n),
    tau_var = rnorm(n),
    cluster = factor(sample(1:3, n, replace = TRUE), levels = 1:3),
    Actual_Drop_15_16 = rnorm(n)  # unrelated to cluster/tau
  )
  cfg <- list(drop_col = "Actual_Drop_15_16", label = "TEST")
  res <- list(d = d, cfg = cfg)

  result <- compute_roc(res, K = 3)
  expect_true(result$auc > 0.3 && result$auc < 0.7)
})
