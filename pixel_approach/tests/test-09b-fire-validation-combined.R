library(testthat)
library(tidyverse)

source("../src/clustering_helpers.R")

# 09b_fire_validation_combined.R has no functions of its own -- flat
# top-level pipeline, structurally identical to
# 08b_beetle_validation_combined.R (join pixel clusters + fire exposure ->
# zonal stats -> relative risk -> ratio_to_c1). Real ~250k-pixel loads
# (fire_flags, pixels) are skipped and replaced with a small synthetic
# `pixels` fixture below (already shaped like the post-join result);
# chart/CSV-writing lines are skipped too. `eval()` uses the default envir
# since nothing here is wrapped in a function.
exprs    <- parse("../src/09b_fire_validation_combined.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|^fire_flags <-|^year_cols <-|^stopifnot\\(|^fire_flags\\$|^pixels <-",
  src_text
)
stop_at <- which(grepl("^write_csv\\(", src_text))[1]
stopifnot(!is.na(stop_at))

# 4 clusters x 10 pixels, burn rate stepping 10/30/50/70% -- hand-
# computable relative risk (Cluster 4 vs 1): 70% / 10% = 7.00x, and
# ratio_to_c1 = c(1, 3, 5, 7). Numeric (not factor) cluster column --
# matches what read_csv() infers for a "1"/"2"/"3"/"4" CSV column (the
# real load this replaces); K <- max(zonal$cluster) errors on an
# unordered factor but works fine on numeric (see
# test-08b-beetle-validation-combined.R for the same issue).
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
