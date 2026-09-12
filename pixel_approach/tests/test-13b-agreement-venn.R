library(testthat)
library(tidyverse)

# Loads only rate()/circle_pts() from 13b_agreement_venn.R -- both pure
# functions with no free-variable dependencies. Real ~250k-pixel data load
# (d <- read_csv(...)) and everything downstream of it (group filters,
# chart) are skipped; the group-filter logic itself (ndvi_cluster == 3 &
# ndwi_cluster != 3, etc.) is the same case-when-style classification
# already covered by test-13-cluster-agreement-validation.R, so isn't
# duplicated here.
exprs    <- parse("../src/13b_agreement_venn.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^d <-|^ndvi_only <-|^ndwi_only <-|^both <-|^ndvi_n <-|^ndvi_fire <-|^ndwi_n <-|^ndwi_fire <-|^both_n <-|^both_fire <-|^cat\\(",
  src_text
)
stop_at <- which(grepl("^p <- ggplot\\(\\)", src_text))[1]
stopifnot(!is.na(stop_at))
for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]])
}

test_that("rate computes n and pct TRUE for a non-empty subset", {
  subset <- tibble(any_mortality = c(TRUE, TRUE, FALSE, FALSE, FALSE))
  result <- rate(subset, "any_mortality")
  expect_equal(unname(result["n"]), 5)
  expect_equal(unname(result["pct"]), 40)
})

test_that("rate returns n=0, pct=0 for an empty subset instead of NaN", {
  subset <- tibble(any_mortality = logical(0))
  result <- rate(subset, "any_mortality")
  expect_equal(unname(result["n"]), 0)
  expect_equal(unname(result["pct"]), 0)
})

test_that("rate ignores NA values via na.rm", {
  subset <- tibble(any_mortality = c(TRUE, NA, FALSE, FALSE))
  result <- rate(subset, "any_mortality")
  expect_equal(unname(result["n"]), 4)  # n is nrow(), NAs still counted
  expect_equal(unname(result["pct"]), round(100 * 1 / 3, 1))  # pct excludes the NA
})

test_that("circle_pts generates n points on a circle of radius r centered at (cx, cy)", {
  # seq(0, 2*pi, length.out=n) spans n-1 steps (both endpoints included,
  # and 0/2*pi coincide) -- with n=5 that's pi/2 per step, so theta hits
  # exactly pi (the far side of the circle) at index 3.
  pts <- circle_pts(cx = 2, cy = -1, r = 3, n = 5)
  expect_equal(nrow(pts), 5)
  expect_equal(pts$x[1], 2 + 3)                          # theta = 0
  expect_equal(pts$y[1], -1, tolerance = 1e-10)
  expect_equal(pts$x[3], 2 - 3, tolerance = 1e-10)        # theta = pi
  expect_equal(pts$y[3], -1, tolerance = 1e-10)
})

test_that("circle_pts defaults to 200 points", {
  pts <- circle_pts(0, 0, 1)
  expect_equal(nrow(pts), 200)
})
