library(testthat)
library(tidyverse)
library(patchwork)

source("../src/clustering_helpers.R")

# Tests for 09_fire_validation.R.

# Globals present before this file; anything newer is removed at the end
.pre_existing_globals <- ls(.GlobalEnv)

# Load the function into .GlobalEnv, skipping real-data loads
exprs    <- parse("../src/09_fire_validation.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|^OUT_CHART|^index_config <-|^fire_flags <-|^year_cols <-|^stopifnot\\(|^fire_flags\\$",
  src_text
)
stop_at <- which(grepl("^results <-", src_text))[1]
stopifnot(!is.na(stop_at))
for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]], envir = .GlobalEnv)
}

# Fixtures supply both burned columns directly

test_that("run_fire_validation computes burn rate and relative risk correctly", {
  # Burned 1/10, 3/10, 5/10 -> RR (3 vs 1) = 5.00
  pixel_ids <- 1:30
  clusters  <- factor(rep(1:3, each = 10), levels = 1:3)
  burned    <- c(rep(TRUE, 1), rep(FALSE, 9),
                 rep(TRUE, 3), rep(FALSE, 7),
                 rep(TRUE, 5), rep(FALSE, 5))

  pixel_path <- tempfile(fileext = ".csv")
  write_csv(tibble(pixel_id = pixel_ids, lon = -120, lat = 38, cluster = clusters),
            pixel_path)

  fire_flags <- tibble(pixel_id = pixel_ids,
                        burned_baseline_00_14 = FALSE,
                        burned_validation_15_16 = burned)

  cfg <- list(pixel_path = pixel_path, label = "TEST",
              out_csv = tempfile(fileext = ".csv"))

  result <- run_fire_validation(cfg, fire_flags)

  expect_equal(result$zonal$pct_burned_validation, c(10, 30, 50))
  expect_equal(result$zonal$ratio_to_c1, c(1, 3, 5), tolerance = 1e-6)
})

test_that("run_fire_validation reports the baseline window separately from the validation window", {
  # Baseline and validation rates must stay separate
  pixel_ids <- 1:10
  pixel_path <- tempfile(fileext = ".csv")
  write_csv(tibble(pixel_id = pixel_ids, lon = -120, lat = 38,
                    cluster = factor(rep(1, 10), levels = 1)),
            pixel_path)
  fire_flags <- tibble(pixel_id = pixel_ids,
                        burned_baseline_00_14 = TRUE,
                        burned_validation_15_16 = c(rep(TRUE, 4), rep(FALSE, 6)))
  cfg <- list(pixel_path = pixel_path, label = "TEST",
              out_csv = tempfile(fileext = ".csv"))

  result <- run_fire_validation(cfg, fire_flags)
  expect_equal(result$zonal$pct_burned_baseline_00_14, 100)
  expect_equal(result$zonal$pct_burned_validation, 40)
})

rm(list = setdiff(ls(.GlobalEnv), .pre_existing_globals), envir = .GlobalEnv)
