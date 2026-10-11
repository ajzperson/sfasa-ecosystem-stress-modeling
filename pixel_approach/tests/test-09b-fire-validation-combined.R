library(testthat)
library(tidyverse)

source("../src/clustering_helpers.R")

# Tests for 09b_fire_validation_combined.R, run on a synthetic pixels table.
exprs    <- parse("../src/09b_fire_validation_combined.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|^fire_flags <-|^year_cols <-|^stopifnot\\(|^fire_flags\\$|^pixels <-",
  src_text
)
stop_at <- which(grepl("^write_csv\\(", src_text))[1]
stopifnot(!is.na(stop_at))

# Burned 10/30/50/70% -> RR (4 vs 1) = 7.00; numeric cluster as read_csv() returns it
n_burned <- c(1, 3, 5, 7)
burned <- unlist(lapply(n_burned, function(n) c(rep(TRUE, n), rep(FALSE, 10 - n))))
pixels <- tibble(
  pixel_id = 1:40,
  lon = -120, lat = 38,
  cluster = rep(1:4, each = 10),
  burned_baseline_00_14 = FALSE,
  burned_validation_15_16 = burned
)

for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]])
}

test_that("combined fire validation computes zonal burn rate per cluster", {
  expect_equal(zonal$pct_burned_validation, c(10, 30, 50, 70))
})

test_that("combined fire validation's relative risk compares worst (K) vs best (1) cluster", {
  expect_equal(K, 4)
  expect_equal(rr_worst_vs_best, 7, tolerance = 1e-6)
})

test_that("combined fire validation's ratio_to_c1 is 1.00x for Cluster 1 and matches relative risk for Cluster K", {
  expect_equal(zonal$ratio_to_c1, c(1, 3, 5, 7), tolerance = 1e-6)
})

test_that("combined fire validation keeps baseline and validation window rates separate", {
  expect_true(all(zonal$pct_burned_baseline_00_14 == 0))
})
