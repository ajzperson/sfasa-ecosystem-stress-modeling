library(testthat)
library(tidyverse)
library(sf)
library(patchwork)

# 13_cluster_agreement_validation.R has no functions of its own (besides
# the chart-only build_group_panel at the very end, not covered) -- flat
# top-level pipeline (inner-join NDVI+NDWI clusters -> classify agreement
# group -> join real-world outcomes -> zonal rates -> relative risk). Real
# ~250k-pixel loads (ndvi, ndwi, fire, beetle) are skipped and replaced
# with small synthetic fixtures below; chart-building lines are skipped
# too. `eval()` uses the default envir since nothing kept here is wrapped
# in a function.
exprs    <- parse("../src/13_cluster_agreement_validation.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "st_read\\(|^ndvi <-|^ndwi <-|^fire <-|^year_cols <-|^fire\\$|^beetle <-",
  src_text
)
stop_at <- which(grepl("^plot_df <-", src_text))[1]
stopifnot(!is.na(stop_at))

# 3 agreement groups x 10 pixels, plus a "Neither" group (should be
# dropped) and one NDVI-only pixel absent from ndwi (should be dropped by
# the inner_join). Hand-computable rates:
#   Both (1-10):         fire 6/10=60%, beetle 3/10=30%
#   NDVI-3 only (11-20):  fire 3/10=30%, beetle 3/10=30%
#   NDWI-3 only (21-30):  fire 2/10=20%, beetle 1/10=10%
#   Neither (31-33):     excluded entirely -- group == "Neither"
# Relative risk (Both vs NDWI-3 only): Fire 60/20=3.00x, Beetle 30/10=3.00x.
ndvi <- tibble(
  pixel_id = c(1:30, 31:33, 999),
  lon = -120, lat = 38,
  cluster_ndvi = c(rep(3, 10), rep(3, 10), rep(1, 10), rep(1, 3), 3)
)
ndwi <- tibble(
  pixel_id = 1:30,  # 999 (NDVI-only) and 31:33 stay out of ndwi on purpose
  cluster_ndwi = c(rep(3, 10), rep(1, 10), rep(3, 10))
)
fire <- tibble(
  pixel_id = 1:33,
  burned_validation = c(rep(TRUE, 6), rep(FALSE, 4),    # Both: 6/10
                         rep(TRUE, 3), rep(FALSE, 7),    # NDVI-3 only: 3/10
                         rep(TRUE, 2), rep(FALSE, 8),    # NDWI-3 only: 2/10
                         rep(FALSE, 3))                  # Neither
)
beetle <- tibble(
  pixel_id = 1:33,
  any_mortality = c(rep(TRUE, 3), rep(FALSE, 7),         # Both: 3/10
                     rep(TRUE, 3), rep(FALSE, 7),         # NDVI-3 only: 3/10
                     rep(TRUE, 1), rep(FALSE, 9),         # NDWI-3 only: 1/10
                     rep(FALSE, 3))                       # Neither
)

for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]])
}

test_that("agreement classification correctly assigns Both/NDVI-3-only/NDWI-3-only and drops Neither", {
  expect_equal(nrow(joined), 30)  # 31:33 (Neither) and 999 (inner_join-only) excluded
  expect_setequal(joined$pixel_id[joined$group == "Both"], 1:10)
  expect_setequal(joined$pixel_id[joined$group == "NDVI-3 only"], 11:20)
  expect_setequal(joined$pixel_id[joined$group == "NDWI-3 only"], 21:30)
})

test_that("inner_join drops the NDVI-only pixel absent from ndwi", {
  expect_false(999 %in% joined$pixel_id)
})

test_that("summary_tbl computes the correct fire/beetle rate per group", {
  s <- summary_tbl
  expect_equal(s$Fire[s$group == "Both"], 60)
  expect_equal(s$Fire[s$group == "NDVI-3 only"], 30)
  expect_equal(s$Fire[s$group == "NDWI-3 only"], 20)
  expect_equal(s$Beetle[s$group == "Both"], 30)
  expect_equal(s$Beetle[s$group == "NDWI-3 only"], 10)
})

test_that("relative risk (Both vs NDWI-3 only) matches hand-computed ratios", {
  expect_equal(rr_fire, 3, tolerance = 1e-6)
  expect_equal(rr_beetle, 3, tolerance = 1e-6)
})
