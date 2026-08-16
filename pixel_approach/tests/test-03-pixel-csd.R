library(testthat)
library(tidyverse)
library(zoo)
library(future.apply)
library(progressr)
options(progressr.enable = TRUE)

# Load only the function definitions from 03_pixel_stl_csd.R (roll_r1,
# roll_var, roll_skew, edge_na_run, run_pixel_csd), stopping before its
# driver block (Sys.setenv/plan/for-loop over real ~1GB matrices). The
# script's own top-level loads (pixel_meta/vpd_mat_full/ppt_mat_full) use
# project-root-relative paths, so wd is temporarily switched there — the
# tests below don't use those real files anyway (synthetic inputs instead).
exprs    <- parse("../src/03_pixel_stl_csd.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
stop_at  <- which(grepl("^Sys\\.setenv\\(VECLIB_MAXIMUM_THREADS", src_text))
stopifnot(length(stop_at) == 1)
old_wd <- getwd()
setwd("../..")
for (e in exprs[seq_len(stop_at - 1)]) eval(e, envir = .GlobalEnv)
setwd(old_wd)

# ----------------------------------------------------------------------------
# Pure helper functions
# ----------------------------------------------------------------------------

test_that("roll_r1 detects strong positive autocorrelation", {
  # A perfectly linear series doesn't give a value near 1 under this
  # formula even in principle -- dividing by length(x) rather than the
  # actual pair count, plus endpoints only contributing to one product
  # each, attenuates it. 0.7 still clearly demonstrates strong positive
  # detection without assuming an unrealistic exact value.
  x <- 1:20
  d <- seq(as.Date("2020-01-01"), by = "month", length.out = 20)
  expect_gt(roll_r1(x, d), 0.7)
})

test_that("roll_r1 excludes Nov->Mar cross-year gaps from the calculation", {
  # Two versions of the same values, differing only in whether the
  # Nov->Mar pair is a genuine 1-month gap or a 4-month (Nov->Mar) gap.
  # If the Nov->Mar pair were wrongly included, the two results would
  # differ; since it's correctly excluded, they must be identical.
  x <- c(1, 5, 2, 8, 3, 9)
  d_within_season  <- as.Date(c("2020-07-01", "2020-08-01", "2020-09-01",
                                "2020-10-01", "2020-11-01", "2020-12-01"))
  d_cross_year_gap <- as.Date(c("2020-07-01", "2020-08-01", "2020-09-01",
                                "2020-10-01", "2020-11-01", "2021-03-01"))
  # Both have the same diff-in-months pattern except the last pair (1 vs 4
  # months apart) -- roll_r1 should treat that last pair as excluded in the
  # second case, giving a DIFFERENT (not equal) result than if it were
  # naively treated as adjacent.
  r_within <- roll_r1(x, d_within_season)
  r_cross  <- roll_r1(x, d_cross_year_gap)
  expect_false(isTRUE(all.equal(r_within, r_cross)))
})

test_that("roll_var equals plain variance", {
  x <- c(2, 4, 4, 4, 5, 5, 7, 9)
  expect_equal(roll_var(x), var(x))
})

test_that("roll_skew detects right-skew as positive (real, un-negated skewness)", {
  # c(1,1,1,1,10) has positive raw (Fisher-Pearson) skewness (long right
  # tail) -- roll_skew returns the real, un-negated value.
  x <- c(1, 1, 1, 1, 10)
  expect_gt(roll_skew(x), 0)
})

test_that("roll_skew is ~0 for a symmetric distribution", {
  x <- c(1, 2, 3, 4, 5, 6, 7)  # symmetric around 4
  expect_equal(roll_skew(x), 0, tolerance = 1e-10)
})

test_that("edge_na_run returns 0 when there's no leading NA", {
  expect_equal(edge_na_run(c(1, 2, NA, 3)), 0L)
})

test_that("edge_na_run counts the leading NA run length", {
  expect_equal(edge_na_run(c(NA, NA, NA, 1, 2)), 3L)
})

test_that("edge_na_run on rev(x) counts the trailing NA run", {
  x <- c(1, 2, 3, NA, NA)
  expect_equal(edge_na_run(rev(x)), 2L)
})

# ----------------------------------------------------------------------------
# run_pixel_csd() integration test on small synthetic data (not the real
# ~250k-pixel x 1GB matrices -- a deterministic 4-year synthetic series
# with a clear seasonal + trend pattern, small enough for STL to converge
# fast while still exercising the full interpolate -> STL -> deseasonalize
# -> climate-correct -> rolling-CSD -> Kendall's tau pipeline end to end.
# ----------------------------------------------------------------------------

test_that("run_pixel_csd processes a normal pixel end-to-end and correctly drops a high-NA pixel", {
  n_months <- 48  # Jan 2000 - Dec 2003
  dates    <- seq(as.Date("2000-01-01"), by = "month", length.out = n_months)
  month_of <- as.integer(format(dates, "%m"))
  t_idx    <- seq_len(n_months)

  # Deterministic seasonal + trend + small secondary-frequency wobble (so
  # STL has real structure to decompose without needing actual randomness).
  seasonal_pattern <- 0.10 * sin(2 * pi * (month_of - 1) / 12)
  wobble           <- 0.01 * sin(2 * pi * t_idx / 5)
  trend            <- 0.001 * t_idx

  pixel_normal  <- 0.5 + seasonal_pattern + trend + wobble
  pixel_high_na <- pixel_normal
  pixel_high_na[sample(n_months, 20)] <- NA  # >30% of 48 months

  mat <- rbind(pixel_normal, pixel_high_na)
  colnames(mat) <- as.character(dates)
  rownames(mat) <- NULL

  vpd_mat <- matrix(0.5 + 0.05 * sin(2 * pi * (month_of - 1) / 12 + 1) + 0.0005 * t_idx,
                    nrow = 2, ncol = n_months, byrow = TRUE)
  ppt_mat <- matrix(50 + 10 * sin(2 * pi * (month_of - 1) / 12 + 2) - 0.05 * t_idx,
                    nrow = 2, ncol = n_months, byrow = TRUE)
  colnames(vpd_mat) <- colnames(mat)
  colnames(ppt_mat) <- colnames(mat)

  pixel_meta <- data.frame(pixel_id = 1:2, lon = c(-120, -120), lat = c(38, 38))

  out_file   <- tempfile(fileext = ".csv")
  r1_file    <- tempfile(fileext = ".rds")
  var_file   <- tempfile(fileext = ".rds")
  skew_file  <- tempfile(fileext = ".rds")
  future::plan(future::multisession, workers = 2)
  run_pixel_csd(mat, "SYNTH-TEST", out_file, r1_file, var_file, skew_file,
                vpd_mat, ppt_mat, pixel_meta,
                WIN = 12, NA_THRESHOLD = 0.30, EDGE_NA_LIMIT = 6, N_WORKERS = 2)
  future::plan(future::sequential)

  result <- read_csv(out_file, show_col_types = FALSE)

  expect_equal(nrow(result), 2)
  expect_false(is.na(result$tau_r1[1]))   # normal pixel processed successfully
  expect_true(is.na(result$tau_r1[2]))    # high-NA pixel dropped
  expect_true(all(result$tau_r1 >= -1 & result$tau_r1 <= 1, na.rm = TRUE))
  expect_true(all(result$tau_var >= -1 & result$tau_var <= 1, na.rm = TRUE))

  # Rolling r1/var/skew series: 2 rows (pixels) x (48 growing-season-ish
  # months - WIN=12 + 1) columns; row 1 (normal pixel) real values, row 2
  # (dropped pixel) all-NA, same shape as row 1 either way.
  r1_matrix   <- readRDS(r1_file)
  var_matrix  <- readRDS(var_file)
  skew_matrix <- readRDS(skew_file)

  expect_equal(nrow(r1_matrix), 2)
  expect_equal(nrow(var_matrix), 2)
  expect_equal(nrow(skew_matrix), 2)
  expect_true(ncol(r1_matrix) > 0)
  expect_equal(ncol(r1_matrix), ncol(var_matrix))
  expect_equal(ncol(r1_matrix), ncol(skew_matrix))
  expect_false(any(is.na(r1_matrix[1, ])))   # normal pixel: fully populated
  expect_true(all(is.na(r1_matrix[2, ])))    # dropped pixel: all-NA placeholder
  expect_true(all(r1_matrix[1, ] >= -1 & r1_matrix[1, ] <= 1))
  # Column names are the window end dates, parseable back to real dates
  expect_false(anyNA(as.Date(colnames(r1_matrix))))
})
