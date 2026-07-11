library(testthat)
library(tidyverse)

# Load only the function definitions from 04_actual_drop.R (safe_min,
# add_actual_drop), stopping before its driver loop over the real
# ndvi_pixel_tau.csv / ndwi_pixel_tau.csv files.
exprs    <- parse("../src/04_actual_drop.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
stop_at  <- which(grepl("^for \\(idx in INDICES\\)", src_text))
stopifnot(length(stop_at) == 1)
for (e in exprs[seq_len(stop_at - 1)]) eval(e, envir = .GlobalEnv)

test_that("safe_min returns NA for an all-NA row instead of Inf", {
  expect_equal(safe_min(c(NA_real_, NA_real_)), NA_real_)
})

test_that("safe_min returns the real minimum when some values are present", {
  expect_equal(safe_min(c(3, NA, 1, NA, 2)), 1)
})

test_that("add_actual_drop computes baseline mean minus 2015-2016 min, growing season only", {
  # 3 years covering baseline (2010), a non-growing-season month to be
  # excluded (Dec), and the collapse window (2015-2016).
  dates <- c("2010-06-01", "2010-07-01",   # baseline growing season
            "2010-12-01",                  # baseline, but NOT growing season -> excluded
            "2015-06-01", "2016-07-01")    # collapse growing season
  mat <- matrix(c(
    0.5, 0.6, 99,   0.3, 0.2,   # pixel 1: baseline mean = mean(0.5,0.6) = 0.55; min(0.3,0.2) = 0.2 -> drop = 0.35
    0.4, 0.4, -99,  0.4, 0.4    # pixel 2: baseline mean = 0.4; min = 0.4 -> drop = 0
  ), nrow = 2, byrow = TRUE)
  colnames(mat) <- dates

  mat_path <- tempfile(fileext = ".rds")
  saveRDS(mat, mat_path)

  tau_path <- tempfile(fileext = ".csv")
  write_csv(data.frame(pixel_id = 1:2, tau_r1 = c(0.1, 0.2)), tau_path)

  add_actual_drop(tau_path, mat_path, "Actual_Drop_15_16", "TEST")

  result <- read_csv(tau_path, show_col_types = FALSE)
  expect_equal(result$Actual_Drop_15_16, c(0.35, 0.00), tolerance = 1e-6)
  # The Dec value (99 / -99) must not have leaked into either baseline
  # mean or collapse min -- if it had, pixel 1's drop would be very
  # different (baseline mean would be ~33 instead of 0.55).
})

test_that("add_actual_drop's alignment check passes despite pixel_id being double vs seq_len's integer", {
  # Regression test: the alignment stopifnot must use == (value equality),
  # not identical(), because pixel_id round-trips through CSV as double
  # while seq_len() returns integer -- identical() would false-positive
  # reject correctly-aligned data purely due to type.
  dates <- c("2010-06-01", "2015-06-01")
  mat <- matrix(c(0.5, 0.3, 0.4, 0.4), nrow = 2, byrow = TRUE)
  colnames(mat) <- dates
  mat_path <- tempfile(fileext = ".rds")
  saveRDS(mat, mat_path)

  tau_path <- tempfile(fileext = ".csv")
  write_csv(data.frame(pixel_id = 1:2, tau_r1 = c(0.1, 0.2)), tau_path)  # pixel_id becomes double after CSV round-trip

  expect_no_error(add_actual_drop(tau_path, mat_path, "Actual_Drop_15_16", "TEST"))
})

test_that("add_actual_drop's alignment check catches a genuine pixel_id misalignment", {
  dates <- c("2010-06-01", "2015-06-01")
  mat <- matrix(c(0.5, 0.3, 0.4, 0.4), nrow = 2, byrow = TRUE)
  colnames(mat) <- dates
  mat_path <- tempfile(fileext = ".rds")
  saveRDS(mat, mat_path)

  tau_path <- tempfile(fileext = ".csv")
  write_csv(data.frame(pixel_id = c(2, 1), tau_r1 = c(0.1, 0.2)), tau_path)  # reversed order

  expect_error(add_actual_drop(tau_path, mat_path, "Actual_Drop_15_16", "TEST"))
})

test_that("add_actual_drop is idempotent -- re-running overwrites rather than duplicating the column", {
  dates <- c("2010-06-01", "2015-06-01")
  mat <- matrix(c(0.5, 0.3), nrow = 1)
  colnames(mat) <- dates
  mat_path <- tempfile(fileext = ".rds")
  saveRDS(mat, mat_path)

  tau_path <- tempfile(fileext = ".csv")
  write_csv(data.frame(pixel_id = 1, tau_r1 = 0.1), tau_path)

  add_actual_drop(tau_path, mat_path, "Actual_Drop_15_16", "TEST")
  add_actual_drop(tau_path, mat_path, "Actual_Drop_15_16", "TEST")

  result <- read_csv(tau_path, show_col_types = FALSE)
  expect_equal(sum(names(result) == "Actual_Drop_15_16"), 1)
  expect_equal(result$Actual_Drop_15_16, 0.2, tolerance = 1e-6)
})
