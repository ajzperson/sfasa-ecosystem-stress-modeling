library(testthat)
library(tidyverse)

source("../src/clustering_helpers.R")

# 08b_beetle_validation_combined.R has no functions of its own -- flat
# top-level pipeline (join pixel clusters + beetle mortality -> zonal
# stats -> relative risk -> ratio_to_c1). Real ~250k-pixel loads
# (beetle_flags, pixels) are skipped and replaced with a small synthetic
# `pixels` fixture below (already shaped like the post-join result, since
# the join itself is skipped); chart/CSV-writing lines are skipped too.
# `eval()` uses the default envir (this file's own environment) since
# nothing here is wrapped in a function -- same reasoning as
# test-06b-clustering-combined.R.
exprs    <- parse("../src/08b_beetle_validation_combined.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|^beetle_flags <-|^win_cols <-|^stopifnot\\(|^beetle_flags\\$|^pixels <-",
  src_text
)
stop_at <- which(grepl("^write_csv\\(", src_text))[1]
stopifnot(!is.na(stop_at))

# 4 clusters x 10 pixels, mortality rate stepping 10/30/60/90% -- hand-
# computable relative risk (Cluster 4 vs 1): 90% / 10% = 9.00x, and
# ratio_to_c1 = c(1, 3, 6, 9).
n_mortality <- c(1, 3, 6, 9)
mortality <- unlist(lapply(n_mortality, function(n) c(rep(TRUE, n), rep(FALSE, 10 - n))))
pixels <- tibble(
  pixel_id = 1:40,
  lon = -120, lat = 38,
  # Plain numeric, not factor -- matches what read_csv() actually infers
  # for a "1"/"2"/"3"/"4" column round-tripped through write_csv() (the
  # real load this replaces); K <- max(zonal$cluster) below errors on an
  # unordered factor but works fine on numeric.
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
  # Every mortality pixel has tpa=3 by construction; mean_tpa_affected must
  # be exactly 3 regardless of cluster (not diluted by the tpa=0 pixels).
  expect_true(all(zonal$mean_tpa_affected == 3))
})
