library(testthat)
library(tidyverse)
library(patchwork)

source("../src/clustering_helpers.R")

# Snapshot of .GlobalEnv before this file writes anything into it (see
# cleanup at the bottom -- same reasoning as test-06-clustering.R).
.pre_existing_globals <- ls(.GlobalEnv)

# Load only VALIDATION_YEARS/year_label and run_fire_validation() from
# 09_fire_validation.R, stopping before the driver code that loads real
# ~250k-pixel data. run_fire_validation() references year_label as a free
# variable and calls build_cluster_style() (from clustering_helpers.R,
# sourced above) -- both must be visible via its closure, so its
# definition is eval'd into .GlobalEnv (same pattern as
# test-08-beetle-validation.R's run_beetle_validation()).
exprs    <- parse("../src/09_fire_validation.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|^OUT_CHART|^index_config <-|^fire_flags <-|^year_cols <-|^stopifnot\\(|^fire_flags\\$",
  src_text
)
stop_at <- which(grepl("^results <-", src_text))[1]
stopifnot(!is.na(stop_at))
for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]], envir = .GlobalEnv)
}

# fire_flags, as run_fire_validation() receives it, already has
# burned_baseline_00_14 (loaded as-is from pixel_fire.csv) and
# burned_validation_15_16 (rowSums OR-reduction across the validation-
# window burned_<year> columns) -- that reduction happens in the driver
# code above, which is skipped, so fixtures here supply both columns
# directly.

test_that("run_fire_validation computes burn rate and relative risk correctly", {
  # 3 clusters x 10 pixels. Cluster 1: 1/10 burned. Cluster 2: 3/10.
  # Cluster 3: 5/10. Hand-computable relative risk (Cluster 3 vs 1):
  # 50% / 10% = 5.00x.
  pixel_ids <- 1:30
  clusters  <- factor(rep(1:3, each = 10), levels = 1:3)
  burned    <- c(rep(TRUE, 1), rep(FALSE, 9),
                 rep(TRUE, 3), rep(FALSE, 7),
                 rep(TRUE, 5), rep(FALSE, 5))

  pixel_path <- tempfile(fileext = ".csv")
  write_csv(tibble(pixel_id = pixel_ids, lon = -120, lat = 38, cluster = clusters),
            pixel_path)

  fire_flags <- tibble(pixel_id = pixel_ids,
                        burned_baseline_00_14 = FALSE,
                        burned_validation_15_16 = burned)

  cfg <- list(pixel_path = pixel_path, label = "TEST",
              out_csv = tempfile(fileext = ".csv"))

  result <- run_fire_validation(cfg, fire_flags)

  expect_equal(result$zonal$pct_burned_validation, c(10, 30, 50))
  expect_equal(result$zonal$ratio_to_c1, c(1, 3, 5), tolerance = 1e-6)
})

test_that("run_fire_validation reports the baseline window separately from the validation window", {
  # All 10 pixels burned in the baseline window (pre-clustering, should be
  # excluded upstream by 06_clustering.R -- this is just a sanity readout,
  # per the file's own header comment), but only 4 in the validation
  # window. The two rates must not be conflated.
  pixel_ids <- 1:10
  pixel_path <- tempfile(fileext = ".csv")
  write_csv(tibble(pixel_id = pixel_ids, lon = -120, lat = 38,
                    cluster = factor(rep(1, 10), levels = 1)),
            pixel_path)
  fire_flags <- tibble(pixel_id = pixel_ids,
                        burned_baseline_00_14 = TRUE,
                        burned_validation_15_16 = c(rep(TRUE, 4), rep(FALSE, 6)))
  cfg <- list(pixel_path = pixel_path, label = "TEST",
              out_csv = tempfile(fileext = ".csv"))

  result <- run_fire_validation(cfg, fire_flags)
  expect_equal(result$zonal$pct_burned_baseline_00_14, 100)
  expect_equal(result$zonal$pct_burned_validation, 40)
})

rm(list = setdiff(ls(.GlobalEnv), .pre_existing_globals), envir = .GlobalEnv)
