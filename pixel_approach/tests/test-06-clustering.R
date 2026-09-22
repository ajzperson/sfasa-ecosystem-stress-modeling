library(testthat)
library(tidyverse)
library(patchwork)
library(multcompView)
library(ggsignif)
library(cluster)

# run_clustering() calls compute_base_level/residualize_tau/run_kmeans_sweep,
# defined in clustering_helpers.R -- sourced directly here (correct relative
# path from this test's own location) rather than via 06_clustering.R's own
# source("pixel_approach/src/clustering_helpers.R") call, which is
# repo-root-relative and breaks under testthat's working directory. That
# call is skipped below when 06_clustering.R's expressions are evaluated.
source("../src/clustering_helpers.R")

# Snapshot of .GlobalEnv before this file writes anything into it (see
# cleanup at the bottom of this file for why).
.pre_existing_globals <- ls(.GlobalEnv)

# Load only function definitions from 06_clustering.R (combine_panels,
# run_clustering -- its only two real, uncommented functions), stopping
# before the driver code that loads real ~250k-pixel data. Real-data loads
# at the top (sn_boundary, index_config) are skipped -- these tests use
# their own synthetic sn_boundary/cluster_cols/clust_labels instead.
# run_validation/compute_roc (the ANOVA/Tukey/ROC block) are entirely
# commented out in the source, so parsing it never actually defines them --
# not tested here since there's nothing live to test.
exprs    <- parse("../src/06_clustering.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl("^source\\(|st_read\\(|^fire_flags <-|^pixel_meta <-|dir\\.create\\(|^results <-|^map_panels <-|^n_idx|^p_sweep|^ggsave|^p_map|^p_scatter|^p_ts|^val_results|^p_val|^roc_results|^roc_combined|^roc_colors|^auc_label|^p_roc|^for \\(idx in INDICES\\)", src_text)
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

# fire_flags/pixel_meta are normally loaded from real ~250k-pixel data
# (pixel_fire.csv, pixel_metadata.csv); pixel_id 1:60 covers every pixel
# used by the synthetic tau data in the tests below, with no burned pixels
# so none get excluded by run_clustering's fire filter. Assigned into
# .GlobalEnv explicitly, not with plain `<-` -- run_clustering() (above)
# was itself eval'd into .GlobalEnv, so its free-variable lookups resolve
# there, not in whatever environment testthat sources this file into.
assign("fire_flags", tibble(pixel_id = 1:60, burned_baseline_00_14 = FALSE),
       envir = .GlobalEnv)
assign("pixel_meta", tibble(pixel_id = 1:60), envir = .GlobalEnv)

# ----------------------------------------------------------------------------
# combine_panels() -- regression guard: must return a single panel
# unchanged when INDICES has only one entry, not assume both NDVI and NDWI
# always ran.
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
    # tau_r1 and tau_skew are the two features run_clustering() actually
    # clusters on (feat <- transmute(resid_tau_r1, resid_tau_skew) |>
    # scale()) -- both are group-separated so K-means has a clean signal
    # on the full 2D feature space. scale() z-normalizes each column to
    # unit variance independently, so an unstructured column here would be
    # just as influential as a group-separated one post-scaling -- tau_var
    # is left as pure noise since it's genuinely not a clustering feature
    # (excluded by design; see clustering_helpers.R).
    tau_r1   = c(rnorm(n_per_group, -0.8, 0.05), rnorm(n_per_group, 0, 0.05), rnorm(n_per_group, 0.8, 0.05)),
    tau_var  = rnorm(60, 0, 0.1),
    tau_skew = c(rnorm(n_per_group, -0.8, 0.05), rnorm(n_per_group, 0, 0.05), rnorm(n_per_group, 0.8, 0.05)),
    Actual_Drop_15_16 = rnorm(60, 0.1, 0.02)
  )
  tau_path <- tempfile(fileext = ".csv")
  write_csv(d, tau_path)

  # compute_base_level() reads this like a real ndvi_matrix.rds/ndwi_matrix.rds
  # -- rows in pixel_meta$pixel_id order, columns are dates, pre-2015 columns
  # get averaged into each pixel's base_level. Amplitude kept tiny (0.001)
  # so residualize_tau()'s regression against it barely perturbs tau_r1/
  # tau_skew -- large enough would blur the deliberately clean group
  # separation constructed in `d` above.
  mat_path <- tempfile(fileext = ".rds")
  mat <- matrix(rnorm(60 * 3, 0, 0.001), nrow = 60,
                dimnames = list(NULL, c("2010-01-01", "2010-02-01", "2010-03-01")))
  saveRDS(mat, mat_path)

  cfg <- list(tau_path = tau_path, mat_path = mat_path, drop_col = "Actual_Drop_15_16",
             clusters_out = tempfile(fileext = ".csv"), label = "TEST")

  result <- run_clustering(cfg, K = 3, SEED = 1)

  # The 20 lowest-tau_r1/tau_skew pixels (rows 1:20) must all land in
  # cluster "1" (least stressed); the 20 highest (rows 41:60) in cluster
  # "3" (most stressed) -- regardless of kmeans' internal cluster numbering.
  low_group_clusters  <- result$d$cluster[result$d$pixel_id %in% 1:20]
  high_group_clusters <- result$d$cluster[result$d$pixel_id %in% 41:60]
  expect_true(all(low_group_clusters == "1"))
  expect_true(all(high_group_clusters == "3"))
})

# run_clustering()'s closure was forced into .GlobalEnv above, so its
# free-variable lookups (fire_flags, compute_base_level, etc.) would
# resolve there -- but .GlobalEnv persists across test files within one
# testthat::test_dir() run. Without this cleanup, leftover bindings like
# K/SEED/index_config would leak into later test files that guard their
# own defaults with `if (!exists("K")) K <- ...`, silently overriding
# them.
rm(list = setdiff(ls(.GlobalEnv), .pre_existing_globals), envir = .GlobalEnv)
