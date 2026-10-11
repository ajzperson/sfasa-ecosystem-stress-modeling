library(testthat)
library(terra)
library(sf)
library(lubridate)

source("../src/extract_pixel_matrix.R")

# Tests for extract_pixel_matrix() and the 01 scripts' QC masks, on small synthetic rasters.

make_square_boundary <- function(xmin, ymin, xmax, ymax) {
  st_as_sf(st_sfc(
    st_polygon(list(rbind(c(xmin, ymin), c(xmax, ymin),
                          c(xmax, ymax), c(xmin, ymax), c(xmin, ymin)))),
    crs = 4326
  ))
}

# extract_pixel_matrix()

test_that("extract_pixel_matrix produces correct dimensions and monthly means", {
  dates <- as.Date(c("2020-01-15", "2020-01-20", "2020-02-10"))
  r <- rast(nrows = 4, ncols = 4, nlyrs = 3, xmin = 0, xmax = 4, ymin = 0, ymax = 4,
            crs = "EPSG:4326")
  values(r) <- cbind(rep(10, 16), rep(20, 16), rep(30, 16))
  time(r) <- dates

  boundary <- make_square_boundary(0, 0, 4, 4)  # covers the whole raster
  out_mat  <- tempfile(fileext = ".rds")
  out_meta <- tempfile(fileext = ".csv")

  extract_pixel_matrix(r, dates, boundary,
                        matrix_path = out_mat, metadata_path = out_meta,
                        matrix_label = "TEST")

  m    <- readRDS(out_mat)
  meta <- read.csv(out_meta)

  expect_equal(nrow(m), 16)
  expect_equal(colnames(m), c("2020-01-01", "2020-02-01"))
  expect_true(all(m[, "2020-01-01"] == 15))  # mean(10, 20)
  expect_true(all(m[, "2020-02-01"] == 30))
  expect_equal(nrow(meta), nrow(m))
})

test_that("extract_pixel_matrix crops/masks to the boundary, reducing valid pixel count", {
  # Use >= 2 layers: tapp() misreads a single-layer index
  dates <- as.Date(c("2020-01-15", "2020-01-20"))
  r <- rast(nrows = 4, ncols = 4, nlyrs = 2, xmin = 0, xmax = 4, ymin = 0, ymax = 4,
            crs = "EPSG:4326")
  values(r) <- cbind(rep(5, 16), rep(5, 16))
  time(r) <- dates

  boundary <- make_square_boundary(1, 1, 3, 3)  # covers only part of the raster

  m <- readRDS({
    out <- tempfile(fileext = ".rds")
    extract_pixel_matrix(r, dates, boundary, matrix_path = out, matrix_label = "TEST")
    out
  })

  expect_true(nrow(m) < 16)
  expect_true(nrow(m) > 0)
})

test_that("verify_metadata_path passes silently when the grid matches", {
  dates <- as.Date(c("2020-01-15", "2020-01-20"))
  r <- rast(nrows = 2, ncols = 2, nlyrs = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2,
            crs = "EPSG:4326")
  values(r) <- cbind(c(1, 2, 3, 4), c(1, 2, 3, 4))
  time(r) <- dates
  boundary <- make_square_boundary(0, 0, 2, 2)

  meta_path <- tempfile(fileext = ".csv")
  extract_pixel_matrix(r, dates, boundary, matrix_path = tempfile(fileext = ".rds"),
                        metadata_path = meta_path, matrix_label = "REF")

  expect_no_error(
    extract_pixel_matrix(r, dates, boundary, matrix_path = tempfile(fileext = ".rds"),
                          verify_metadata_path = meta_path, matrix_label = "TEST")
  )
})

test_that("verify_metadata_path errors when pixel count doesn't match", {
  dates <- as.Date(c("2020-01-15", "2020-01-20"))

  r_big <- rast(nrows = 4, ncols = 4, nlyrs = 2, xmin = 0, xmax = 4, ymin = 0, ymax = 4,
                crs = "EPSG:4326")
  values(r_big) <- cbind(rep(1, 16), rep(1, 16))
  time(r_big) <- dates
  meta_path <- tempfile(fileext = ".csv")
  extract_pixel_matrix(r_big, dates, make_square_boundary(0, 0, 4, 4),
                        matrix_path = tempfile(fileext = ".rds"),
                        metadata_path = meta_path, matrix_label = "REF")

  r_small <- rast(nrows = 2, ncols = 2, nlyrs = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2,
                  crs = "EPSG:4326")
  values(r_small) <- cbind(c(1, 2, 3, 4), c(1, 2, 3, 4))
  time(r_small) <- dates

  expect_error(
    extract_pixel_matrix(r_small, dates, make_square_boundary(0, 0, 2, 2),
                          matrix_path = tempfile(fileext = ".rds"),
                          verify_metadata_path = meta_path, matrix_label = "TEST"),
    "pixel count doesn't match"
  )
})

test_that("sample_only mode works with n_sample = 1 (regression test)", {
  # n_sample = 1 must not drop the matrix to a vector
  dates <- as.Date(c("2020-01-15", "2020-02-10"))
  r <- rast(nrows = 4, ncols = 4, nlyrs = 2, xmin = 0, xmax = 4, ymin = 0, ymax = 4,
            crs = "EPSG:4326")
  values(r) <- cbind(rep(1, 16), rep(2, 16))
  time(r) <- dates

  expect_error(
    extract_pixel_matrix(r, dates, make_square_boundary(0, 0, 4, 4),
                          matrix_path = tempfile(fileext = ".rds"),
                          matrix_label = "TEST", sample_only = TRUE, n_sample = 1),
    "Sample complete"
  )
})

# NDWI QC mask and formula (inline in 01_extract_pixel_ndwi.R)

test_that("QC bit-extraction (%% 4) matches real MODIS lookup table semantics", {
  # From the MOD09A1 QC lookup: 1073741824 = ideal, 1073741877 = less than ideal
  expect_equal(1073741824 %% 4, 0)
  expect_equal(1073741877 %% 4, 1)
})

test_that("NDWI formula computes correctly and QC-masks bad-quality cells", {
  b02 <- rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:4326")
  values(b02) <- c(0.3, 0.4, 0.2, 0.5)
  b05 <- rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:4326")
  values(b05) <- c(0.1, 0.1, 0.1, 0.1)
  qc <- rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:4326")
  values(qc) <- c(0, 1, 0, 3)  # cells 2 and 4 are bad quality (%% 4 != 0)

  good_qc    <- qc %% 4 == 0
  b02_masked <- mask(b02, good_qc, maskvalue = FALSE)
  b05_masked <- mask(b05, good_qc, maskvalue = FALSE)
  ndwi       <- (b02_masked - b05_masked) / (b02_masked + b05_masked)

  vals <- values(ndwi)[, 1]
  expect_true(is.na(vals[2]))
  expect_true(is.na(vals[4]))
  expect_equal(vals[1], (0.3 - 0.1) / (0.3 + 0.1))
  expect_equal(vals[3], (0.2 - 0.1) / (0.2 + 0.1))
})

# NDVI QC mask (inline in 01_extract_pixel_ndvi.R): drops only bits 11

test_that("VI_Quality bit-extraction (%% 4) matches real MOD13A1 values", {
  # VI_Quality values with bits 0-1 = 00, 01, 10, 11
  expect_equal(2116 %% 4, 0)
  expect_equal(2061 %% 4, 1)
  expect_equal(2062 %% 4, 2)
  expect_equal(2119 %% 4, 3)
})

test_that("NDVI QC mask keeps bits 00/01/10, drops bits 11", {
  ndvi <- rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:4326")
  values(ndvi) <- c(0.5, 0.6, 0.7, 0.8)
  qc <- rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2, crs = "EPSG:4326")
  values(qc) <- c(2116, 2061, 2062, 2119)  # good, check-QA, cloudy, not-produced

  good_qc     <- qc %% 4 <= 2
  ndvi_masked <- mask(ndvi, good_qc, maskvalue = FALSE)

  vals <- values(ndvi_masked)[, 1]
  expect_equal(vals[1], 0.5)   # good -- kept
  expect_equal(vals[2], 0.6)   # check-other-QA -- kept
  expect_equal(vals[3], 0.7)   # cloudy -- kept
  expect_true(is.na(vals[4]))  # not-produced -- masked
})
