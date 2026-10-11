library(testthat)
library(tidyverse)
library(patchwork)

source("../src/clustering_helpers.R")

# Tests for 08_beetle_validation.R.

# Globals present before this file; anything newer is removed at the end
.pre_existing_globals <- ls(.GlobalEnv)

# Load the function into .GlobalEnv, skipping real-data loads
exprs    <- parse("../src/08_beetle_validation.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|^OUT_CHART|^index_config <-|^beetle_flags <-|^win_cols <-|^stopifnot\\(|^beetle_flags\\$",
  src_text
)
stop_at <- which(grepl("^results <-", src_text))[1]
stopifnot(!is.na(stop_at))
for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]], envir = .GlobalEnv)
}

# Fixtures supply tpa and any_mortality directly

test_that("run_beetle_validation computes mortality rate and relative risk correctly", {
  # Mortality 2/10, 5/10, 8/10 -> RR (3 vs 1) = 4.00
  pixel_ids <- 1:30
  clusters  <- factor(rep(1:3, each = 10), levels = 1:3)
  mortality <- c(rep(TRUE, 2), rep(FALSE, 8),   # cluster 1: 2/10
                 rep(TRUE, 5), rep(FALSE, 5),   # cluster 2: 5/10
                 rep(TRUE, 8), rep(FALSE, 2))   # cluster 3: 8/10
  tpa <- ifelse(mortality, 5, 0)

  pixel_path <- tempfile(fileext = ".csv")
  write_csv(tibble(pixel_id = pixel_ids, cluster = clusters), pixel_path)

  beetle_flags <- tibble(pixel_id = pixel_ids, tpa = tpa, any_mortality = mortality)

  cfg <- list(pixel_path = pixel_path, label = "TEST",
              out_csv = tempfile(fileext = ".csv"))

  result <- run_beetle_validation(cfg, beetle_flags)

  expect_equal(result$zonal$pct_with_mortality, c(20, 50, 80))
  expect_equal(result$zonal$n_with_mortality, c(2, 5, 8))
  expect_equal(result$zonal$ratio_to_c1, c(1, 2.5, 4), tolerance = 1e-6)
})

test_that("run_beetle_validation's ratio_to_c1 is exactly 1.00x for Cluster 1 by construction", {
  pixel_ids <- 1:20
  clusters  <- factor(rep(1:2, each = 10), levels = 1:2)
  mortality <- c(rep(TRUE, 3), rep(FALSE, 7), rep(TRUE, 6), rep(FALSE, 4))
  tpa <- ifelse(mortality, 2, 0)

  pixel_path <- tempfile(fileext = ".csv")
  write_csv(tibble(pixel_id = pixel_ids, cluster = clusters), pixel_path)
  beetle_flags <- tibble(pixel_id = pixel_ids, tpa = tpa, any_mortality = mortality)
  cfg <- list(pixel_path = pixel_path, label = "TEST",
              out_csv = tempfile(fileext = ".csv"))

  result <- run_beetle_validation(cfg, beetle_flags)
  expect_equal(result$zonal$ratio_to_c1[result$zonal$cluster == 1], 1)
})

test_that("run_beetle_validation's mean_tpa_affected averages only mortality pixels, not all", {
  # mean_tpa_affected averages only pixels with mortality
  pixel_ids <- 1:2
  pixel_path <- tempfile(fileext = ".csv")
  write_csv(tibble(pixel_id = pixel_ids, cluster = factor(c(1, 1), levels = 1)),
            pixel_path)
  beetle_flags <- tibble(pixel_id = pixel_ids, tpa = c(10, 0),
                          any_mortality = c(TRUE, FALSE))
  cfg <- list(pixel_path = pixel_path, label = "TEST",
              out_csv = tempfile(fileext = ".csv"))

  result <- run_beetle_validation(cfg, beetle_flags)
  expect_equal(result$zonal$mean_tpa_affected, 10)
  expect_equal(result$zonal$mean_tpa_all, 5)
})

rm(list = setdiff(ls(.GlobalEnv), .pre_existing_globals), envir = .GlobalEnv)
