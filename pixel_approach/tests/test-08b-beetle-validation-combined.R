library(testthat)
library(tidyverse)

source("../src/clustering_helpers.R")

# Tests for 08b_beetle_validation_combined.R, run on a synthetic pixels table.
exprs    <- parse("../src/08b_beetle_validation_combined.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|^beetle_flags <-|^win_cols <-|^stopifnot\\(|^beetle_flags\\$|^pixels <-",
  src_text
)
stop_at <- which(grepl("^write_csv\\(", src_text))[1]
stopifnot(!is.na(stop_at))

# Mortality 10/30/60/90% -> RR (4 vs 1) = 9.00
n_mortality <- c(1, 3, 6, 9)
mortality <- unlist(lapply(n_mortality, function(n) c(rep(TRUE, n), rep(FALSE, 10 - n))))
pixels <- tibble(
  pixel_id = 1:40,
  lon = -120, lat = 38,
  # Numeric, as read_csv() returns it
  cluster = rep(1:4, each = 10),
  tpa = ifelse(mortality, 3, 0),
  any_mortality = mortality
)

for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]])
}

test_that("combined beetle validation computes zonal mortality rate per cluster", {
  expect_equal(zonal$pct_with_mortality, c(10, 30, 60, 90))
  expect_equal(zonal$n_with_mortality, n_mortality)
})

test_that("combined beetle validation's relative risk compares worst (K) vs best (1) cluster", {
  expect_equal(K, 4)
  expect_equal(rr_worst_vs_best, 9, tolerance = 1e-6)
})

test_that("combined beetle validation's ratio_to_c1 is 1.00x for Cluster 1 and matches relative risk for Cluster K", {
  expect_equal(zonal$ratio_to_c1, c(1, 3, 6, 9), tolerance = 1e-6)
})

test_that("combined beetle validation's mean_tpa_affected averages only mortality pixels", {
  expect_true(all(zonal$mean_tpa_affected == 3))
})
