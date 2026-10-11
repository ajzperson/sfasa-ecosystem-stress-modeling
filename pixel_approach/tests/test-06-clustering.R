library(testthat)
library(tidyverse)
library(patchwork)
library(cluster)

# Tests for 06_clustering.R. Helpers are sourced by relative path; the script's own source() is skipped.
source("../src/clustering_helpers.R")

# Globals present before this file; anything newer is removed at the end
.pre_existing_globals <- ls(.GlobalEnv)

# Load only the function definitions, skipping real-data loads
exprs    <- parse("../src/06_clustering.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl("^source\\(|st_read\\(|^fire_flags <-|^pixel_meta <-|dir\\.create\\(|^results <-|^map_panels <-|^n_idx|^p_sweep|^ggsave|^p_map|^p_scatter|^p_ts|^val_results|^p_val|^roc_results|^roc_combined|^roc_colors|^auc_label|^p_roc|^for \\(idx in INDICES\\)", src_text)
stop_at <- which(grepl("^for \\(idx in INDICES\\)", src_text))[1]
stopifnot(!is.na(stop_at))
for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]], envir = .GlobalEnv)
}

K <- 3
cluster_cols <- setNames(c("#2166AC", "#92C5DE", "#D73027"), c("1", "2", "3"))
clust_labels <- c("1 - no CSD", "2 - partial", "3 - full CSD")

# Assigned into .GlobalEnv, where run_clustering() looks them up
assign("fire_flags", tibble(pixel_id = 1:60, burned_baseline_00_14 = FALSE),
       envir = .GlobalEnv)
assign("pixel_meta", tibble(pixel_id = 1:60), envir = .GlobalEnv)

# combine_panels()
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

# run_clustering(): three well-separated groups must map to clusters 1..3 in order
test_that("run_clustering relabels clusters by ascending composite stress", {
  set.seed(1)
  n_per_group <- 20
  d <- data.frame(
    pixel_id = 1:60,
    lon = -120, lat = 38,
    tau_r1   = c(rnorm(n_per_group, -0.8, 0.05), rnorm(n_per_group, 0, 0.05), rnorm(n_per_group, 0.8, 0.05)),
    tau_var  = rnorm(60, 0, 0.1),
    tau_skew = c(rnorm(n_per_group, -0.8, 0.05), rnorm(n_per_group, 0, 0.05), rnorm(n_per_group, 0.8, 0.05)),
    Actual_Drop_15_16 = rnorm(60, 0.1, 0.02)
  )
  tau_path <- tempfile(fileext = ".csv")
  write_csv(d, tau_path)

  # Fake index matrix; tiny values so residualizing barely changes tau
  mat_path <- tempfile(fileext = ".rds")
  mat <- matrix(rnorm(60 * 3, 0, 0.001), nrow = 60,
                dimnames = list(NULL, c("2010-01-01", "2010-02-01", "2010-03-01")))
  saveRDS(mat, mat_path)

  cfg <- list(tau_path = tau_path, mat_path = mat_path, drop_col = "Actual_Drop_15_16",
             clusters_out = tempfile(fileext = ".csv"), label = "TEST")

  result <- run_clustering(cfg, K = 3, SEED = 1)

  low_group_clusters  <- result$d$cluster[result$d$pixel_id %in% 1:20]
  high_group_clusters <- result$d$cluster[result$d$pixel_id %in% 41:60]
  expect_true(all(low_group_clusters == "1"))
  expect_true(all(high_group_clusters == "3"))
})

# Remove globals this file added so they don't leak into later test files
rm(list = setdiff(ls(.GlobalEnv), .pre_existing_globals), envir = .GlobalEnv)
