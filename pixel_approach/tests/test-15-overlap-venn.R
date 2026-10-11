library(testthat)
library(tidyverse)
library(sf)

# Tests for the helpers and group filters in 15_overlap_venn.R (geometry and chart skipped).
exprs    <- parse("../src/15_overlap_venn.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl("^d <-", src_text)
stop_at <- which(grepl("^poly_A <-", src_text))[1]
stopifnot(!is.na(stop_at))

# 4 pixels each: Both, O1-only, O2-only, Neither
d <- tibble(
  pixel_id = 1:16,
  ndvi_cluster = c(rep(3, 12), rep(1, 4)),
  ndwi_cluster = c(rep(3, 12), rep(1, 4)),
  in_agreement_both = c(rep(TRUE, 4), rep(TRUE, 4), rep(FALSE, 4), rep(FALSE, 4)),
  in_combined_top   = c(rep(TRUE, 4), rep(FALSE, 4), rep(TRUE, 4), rep(FALSE, 4)),
  any_mortality = c(rep(TRUE, 2), rep(FALSE, 2),      # Both: 2/4 = 50%
                     rep(TRUE, 1), rep(FALSE, 3),      # O1-only: 1/4 = 25%
                     rep(FALSE, 4),                    # O2-only: 0/4 = 0%
                     rep(FALSE, 4)),                   # Neither
  burned_validation = c(rep(TRUE, 3), rep(FALSE, 1),   # Both: 3/4 = 75%
                         rep(TRUE, 1), rep(FALSE, 3),   # O1-only: 1/4 = 25%
                         rep(TRUE, 2), rep(FALSE, 2),   # O2-only: 2/4 = 50%
                         rep(FALSE, 4))                 # Neither
)

for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]])
}

test_that("O1/O2/both group sizes match the agreement/combined-top overlap definitions", {
  expect_equal(nrow(O1), 8)
  expect_equal(nrow(O2), 8)
  expect_equal(nrow(both), 4)
  expect_equal(n_ndvi3, 12)
  expect_equal(n_ndwi3, 12)
})

test_that("O1/O2/both rates match hand-computed values (O1 and O2 both include Both)", {
  expect_equal(unname(both_n["pct"]), 50)
  expect_equal(unname(both_fire["pct"]), 75)
  expect_equal(unname(o1_n["pct"]), round(100 * 3 / 8, 1))    # (2+1)/8
  expect_equal(unname(o2_n["pct"]), round(100 * 2 / 8, 1))    # (2+0)/8
})

test_that("deg2rad converts degrees to radians", {
  expect_equal(deg2rad(180), pi)
  expect_equal(deg2rad(90), pi / 2)
  expect_equal(deg2rad(0), 0)
})

test_that("rate and circle_pts (duplicated from 14b_agreement_venn.R) behave the same in this file", {
  expect_equal(unname(rate(tibble(x = c(TRUE, FALSE)), "x")["pct"]), 50)
  pts <- circle_pts(0, 0, 2, n = 5)
  expect_equal(nrow(pts), 5)
  expect_equal(pts$x[1], 2)
})
