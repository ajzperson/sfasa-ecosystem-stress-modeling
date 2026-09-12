library(testthat)
library(tidyverse)
library(patchwork)
library(zoo)

# Loads cluster_cols/clust_labels/WIN/build_level_trajectory/build_panel
# from 12_csd_rolling_timeline.R, stopping before the driver code that
# loads real ~250k-pixel data (index_config, results <- lapply(...)).
# build_cluster_trajectories() itself isn't covered -- it reads
# cfg$clusters_path/r1_path/skew_path (parameterizable) but also calls
# build_level_trajectory(), which reads vpd_matrix.rds/ppt_matrix.rds from
# hardcoded repo-root-relative paths, not exposed via cfg; those two are
# faked below as real files at that literal (harmless, nested, cleaned-up)
# path instead, specifically so build_level_trajectory() -- the most
# complex, matrix-arithmetic-heavy function in this file -- gets covered.
exprs    <- parse("../src/12_csd_rolling_timeline.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl("^index_config <-", src_text)
stop_at <- which(grepl("^results <-", src_text))[1]
stopifnot(!is.na(stop_at))
for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]])
}

test_that("build_panel filters to the requested statistic and colors by cluster", {
  traj <- tibble(
    date = as.Date(c("2010-01-01", "2010-02-01", "2010-01-01", "2010-02-01")),
    cluster = factor(c(1, 1, 2, 2), levels = 1:3),
    statistic = c("r1", "r1", "skew", "skew"),
    value = c(0.1, 0.2, -0.1, -0.2)
  )
  p <- build_panel(traj, "r1", "Rolling r1", "TEST LABEL")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "TEST LABEL")
  # Only the 2 "r1" rows should reach the plot -- "skew" rows filtered out.
  expect_equal(nrow(p$data), 2)
})

test_that("build_level_trajectory reconstructs clim_resid and rolls it correctly", {
  # 3 pixels (one per cluster), 72 months (2009-01 to 2014-12) -- exactly
  # 54 growing-season (Mar-Nov) months, i.e. exactly WIN, so there's
  # exactly one rolling window per cluster with a hand-computable value.
  # seas_* = 0 and VPD/PPT = 0 make predicted = b0 and y_deseas = mat_sub
  # unchanged, so clim_resid = mat_sub - b0 = mat_sub (b0 = 0) -- a
  # constant series per pixel, so its rolling mean is exactly that
  # constant: cluster k's level = pixel k's constant value (1, 2, 3).
  col_dates <- format(seq(as.Date("2009-01-01"), as.Date("2014-12-01"), by = "month"))
  n_months  <- length(col_dates)
  stopifnot(n_months == 72)

  mat <- matrix(c(rep(1, n_months), rep(2, n_months), rep(3, n_months)),
                nrow = 3, byrow = TRUE, dimnames = list(NULL, col_dates))
  mat_path <- tempfile(fileext = ".rds")
  saveRDS(mat, mat_path)

  zero_mat <- matrix(0, nrow = 3, ncol = n_months, dimnames = list(NULL, col_dates))

  # vpd_matrix.rds/ppt_matrix.rds are read from this literal hardcoded
  # path inside build_level_trajectory -- created here as real files
  # (nested under this test's own cwd, harmless) and removed afterward.
  fixture_dir <- "pixel_approach/data/processed"
  dir_existed <- dir.exists(fixture_dir)
  dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(zero_mat, file.path(fixture_dir, "vpd_matrix.rds"))
  saveRDS(zero_mat, file.path(fixture_dir, "ppt_matrix.rds"))
  on.exit({
    if (!dir_existed) unlink("pixel_approach", recursive = TRUE)
  }, add = TRUE)

  seas_cols <- setNames(as.list(rep(0, 12)), sprintf("seas_%02d", 1:12))
  clusters <- as_tibble(c(list(pixel_id = 1:3, cluster = 1:3, b0 = 0, b1 = 0, b2 = 0),
                          seas_cols))

  cfg <- list(mat_path = mat_path)
  result <- build_level_trajectory(cfg, clusters)

  expect_equal(nrow(result), 3)  # exactly 1 window x 3 clusters
  expect_equal(result$value[order(result$cluster)], c(1, 2, 3))
  expect_true(all(result$statistic == "level"))
  expect_equal(unique(result$date), as.Date("2014-11-01"))  # last growing-season month
})
