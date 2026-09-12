library(testthat)
library(tidyverse)
library(patchwork)
library(cluster)

source("../src/clustering_helpers.R")

# 06b_clustering_combined.R has no functions of its own -- it's a flat
# top-level pipeline (inner-join NDVI+NDWI tau -> fire exclusion ->
# residualize -> 4-feature K-means -> stress-rank relabel -> profile).
# Real ~250k-pixel loads (fire_flags, pixel_meta, ndvi_tau/ndwi_tau,
# base_level_ndvi/ndwi) are skipped and replaced with small synthetic
# fixtures below; chart/CSV-writing lines are skipped too. `eval()` here
# uses the default envir (this file's own environment, not .GlobalEnv --
# unlike test-06-clustering.R, nothing here is wrapped in a function with
# a forced-different closure), so the resulting `d`/`profile` are directly
# visible by name in the test_that() blocks below.
exprs    <- parse("../src/06b_clustering_combined.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|st_read\\(|^fire_flags <-|^pixel_meta <-|^ndvi_tau <-|^ndwi_tau <-|^base_level_ndvi <-|^base_level_ndwi <-|^dir\\.create\\(|^write_csv\\(|^p_sweep <-|^ggsave\\(|^p_map <-",
  src_text
)

# 4 pixel groups (10 pixels each, pixel_id blocks 1-10/11-20/21-30/31-40),
# all 4 features (tau_r1/tau_skew x NDVI/NDWI) group-separated the same
# way -- these are exactly the features run_kmeans_sweep()'s 4D `feat`
# clusters on, so K-means recovers these blocks with a clean, deterministic
# stress-rank ordering. Two extra NDVI-only pixels (41, 42) test that the
# inner_join drops pixels missing from either index. Pixel 5 is flagged
# burned to test the fire-exclusion filter.
set.seed(1)
group_means <- c(-1.2, -0.4, 0.4, 1.2)
make_group_col <- function() unlist(lapply(group_means, function(m) rnorm(10, m, 0.05)))

ndvi_tau <- tibble(
  pixel_id  = 1:42,
  lon = -120, lat = 38,
  tau_r1    = c(make_group_col(), rnorm(2, 0, 0.05)),
  tau_skew  = c(make_group_col(), rnorm(2, 0, 0.05))
)
ndwi_tau <- tibble(
  pixel_id  = 1:40,
  tau_r1    = make_group_col(),
  tau_skew  = make_group_col()
)

fire_flags <- tibble(pixel_id = 1:40, burned_baseline_00_14 = (pixel_id == 5))

# Tiny amplitude so residualize_tau()'s regression barely perturbs the
# deliberately clean group separation above (see test-06-clustering.R for
# the same reasoning).
base_level_ndvi <- tibble(pixel_id = 1:40, base_level_ndvi = rnorm(40, 0, 0.001))
base_level_ndwi <- tibble(pixel_id = 1:40, base_level_ndwi = rnorm(40, 0, 0.001))

for (i in seq_along(exprs)) {
  if (!skip_lines[i]) eval(exprs[[i]])
}

test_that("combined clustering keeps only pixels present in both NDVI and NDWI", {
  expect_false(any(c(41, 42) %in% d$pixel_id))
  expect_equal(sort(d$pixel_id), setdiff(1:40, 5))  # pixel 5 excluded (burned)
})

test_that("combined clustering excludes burned-baseline pixels", {
  expect_false(5 %in% d$pixel_id)
  expect_equal(nrow(d), 39)
})

test_that("combined clustering relabels the 4 synthetic groups by ascending stress", {
  expect_true(all(d$cluster[d$pixel_id %in% setdiff(1:10, 5)] == "1"))
  expect_true(all(d$cluster[d$pixel_id %in% 11:20] == "2"))
  expect_true(all(d$cluster[d$pixel_id %in% 21:30] == "3"))
  expect_true(all(d$cluster[d$pixel_id %in% 31:40] == "4"))
})

test_that("cluster profile summarizes every cluster with correct group sizes", {
  expect_equal(nrow(profile), 4)
  expect_equal(profile$n, c(9, 10, 10, 10))
})
