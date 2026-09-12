library(testthat)
library(tidyverse)
library(patchwork)

source("../src/clustering_helpers.R")

# Snapshot of .GlobalEnv before this file writes anything into it (see
# cleanup at the bottom -- same reasoning as test-06-clustering.R).
.pre_existing_globals <- ls(.GlobalEnv)

# Load only VALIDATION_YEARS/year_label and run_beetle_validation() from
# 08_beetle_validation.R, stopping before the driver code that loads real
# ~250k-pixel data (index_config, beetle_flags). run_beetle_validation()
# references year_label as a free variable and calls build_cluster_style()
# (from clustering_helpers.R, sourced above) -- both must be visible via
# its closure, so its definition is eval'd into .GlobalEnv, matching
# test-06-clustering.R's run_clustering() pattern.
exprs    <- parse("../src/08_beetle_validation.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|^OUT_CHART|^index_config <-|^beetle_flags <-|^win_cols <-|^stopifnot\\(|^beetle_flags\\$",
  src_text
)
stop_at <- which(grepl("^results <-", src_text))[1]
stopifnot(!is.na(stop_at))
for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]], envir = .GlobalEnv)
}

# beetle_flags, as run_beetle_validation() receives it, already has `tpa`
# (pmax across the validation-window tpa_<year> columns) and `any_mortality`
# (tpa > 0) -- that reduction happens in the driver code above, which is
# skipped, so fixtures here supply those two columns directly.

test_that("run_beetle_validation computes mortality rate and relative risk correctly", {
  # 3 clusters x 10 pixels. Cluster 1: 2/10 mortality (tpa=5 when present).
  # Cluster 2: 5/10 mortality. Cluster 3: 8/10 mortality. Hand-computable
  # relative risk (Cluster 3 vs 1): 80% / 20% = 4.00x.
  pixel_ids <- 1:30
  clusters  <- factor(rep(1:3, each = 10), levels = 1:3)
  mortality <- c(rep(TRUE, 2), rep(FALSE, 8),   # cluster 1: 2/10
                 rep(TRUE, 5), rep(FALSE, 5),   # cluster 2: 5/10
                 rep(TRUE, 8), rep(FALSE, 2))   # cluster 3: 8/10
  tpa <- ifelse(mortality, 5, 0)

  pixel_path <- tempfile(fileext = ".csv")
  write_csv(tibble(pixel_id = pixel_ids, cluster = clusters), pixel_path)

  beetle_flags <- tibble(pixel_id = pixel_ids, tpa = tpa, any_mortality = mortality)

  cfg <- list(pixel_path = pixel_path, label = "TEST",
              out_csv = tempfile(fileext = ".csv"))

  result <- run_beetle_validation(cfg, beetle_flags)

  expect_equal(result$zonal$pct_with_mortality, c(20, 50, 80))
  expect_equal(result$zonal$n_with_mortality, c(2, 5, 8))
  # rr_c3_vs_c1 isn't returned directly, but ratio_to_c1 (in zonal) is the
  # same computation (pct_with_mortality / cluster-1's pct_with_mortality).
  expect_equal(result$zonal$ratio_to_c1, c(1, 2.5, 4), tolerance = 1e-6)
})

test_that("run_beetle_validation's ratio_to_c1 is exactly 1.00x for Cluster 1 by construction", {
  pixel_ids <- 1:20
  clusters  <- factor(rep(1:2, each = 10), levels = 1:2)
  mortality <- c(rep(TRUE, 3), rep(FALSE, 7), rep(TRUE, 6), rep(FALSE, 4))
  tpa <- ifelse(mortality, 2, 0)

  pixel_path <- tempfile(fileext = ".csv")
  write_csv(tibble(pixel_id = pixel_ids, cluster = clusters), pixel_path)
  beetle_flags <- tibble(pixel_id = pixel_ids, tpa = tpa, any_mortality = mortality)
  cfg <- list(pixel_path = pixel_path, label = "TEST",
              out_csv = tempfile(fileext = ".csv"))

  result <- run_beetle_validation(cfg, beetle_flags)
  expect_equal(result$zonal$ratio_to_c1[result$zonal$cluster == 1], 1)
})

test_that("run_beetle_validation's mean_tpa_affected averages only mortality pixels, not all", {
  # Cluster with 2 pixels: one has mortality (tpa=10), one has none (tpa=0).
  # mean_tpa_affected should be 10 (mortality pixels only), not 5 (all pixels).
  pixel_ids <- 1:2
  pixel_path <- tempfile(fileext = ".csv")
  write_csv(tibble(pixel_id = pixel_ids, cluster = factor(c(1, 1), levels = 1)),
            pixel_path)
  beetle_flags <- tibble(pixel_id = pixel_ids, tpa = c(10, 0),
                          any_mortality = c(TRUE, FALSE))
  cfg <- list(pixel_path = pixel_path, label = "TEST",
              out_csv = tempfile(fileext = ".csv"))

  result <- run_beetle_validation(cfg, beetle_flags)
  expect_equal(result$zonal$mean_tpa_affected, 10)
  expect_equal(result$zonal$mean_tpa_all, 5)
})

rm(list = setdiff(ls(.GlobalEnv), .pre_existing_globals), envir = .GlobalEnv)
