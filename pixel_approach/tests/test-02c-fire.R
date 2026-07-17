library(testthat)
library(terra)

# 02c_pixel_fire.R is a top-level script, not a set of functions, so these
# tests replicate its key logic patterns on small synthetic rasters rather
# than sourcing the script against the real multi-GB MTBS files.

BURNED_SEVERITY <- c(3, 4)  # moderate, high -- same as the real script

test_that("severity reclassification keeps only moderate/high (3/4)", {
  # MTBS codes: 1/2 = low severity, 3 = moderate, 4 = high, 5 = increased
  # greenness, 6 = non-processing mask
  r <- rast(nrows = 1, ncols = 6, xmin = 0, xmax = 6, ymin = 0, ymax = 1, crs = "EPSG:32610")
  values(r) <- c(1, 2, 3, 4, 5, 6)

  r_bin <- app(r, function(v) as.numeric(v %in% BURNED_SEVERITY))
  vals  <- values(r_bin)[, 1]

  expect_equal(vals, c(0, 0, 1, 1, 0, 0))
})

test_that("explicit buffer polygon + extract(fun=max) aggregates correctly", {
  # 10x10 raster, cell values 1:100. A point at the center with a 250-unit
  # buffer should pick up the true max value among nearby cells, not just
  # the value at the point's own cell.
  r <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, crs = "EPSG:32610")
  values(r) <- 1:100
  pts <- vect(data.frame(x = 5, y = 5), geom = c("x", "y"), crs = "EPSG:32610")

  buf <- terra::buffer(pts, width = 250)
  val <- terra::extract(r, buf, fun = max)[, 2]

  expect_equal(val, 100)
})

test_that("build_year_raster-style reclassify+resample-to-template preserves burned cells", {
  # Mirrors 02c's build_year_raster(): reclassify to binary, then
  # crop/resample onto a common template grid.
  r <- rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4, crs = "EPSG:32610")
  values(r) <- c(1, 3, 2, 4,
                 1, 1, 1, 1,
                 1, 1, 1, 1,
                 1, 1, 1, 1)
  r_bin <- app(r, function(v) as.numeric(v %in% BURNED_SEVERITY))

  template <- rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4, crs = "EPSG:32610")
  r_bin <- crop(r_bin, template, extend = TRUE)
  r_bin <- resample(r_bin, template, method = "max")
  r_bin[is.na(r_bin)] <- 0

  vals <- values(r_bin)[, 1]
  expect_equal(sum(vals), 2)  # only the two cells with codes 3 and 4
})

test_that("build_burned_raster-style OR-across-years combines multiple binary rasters", {
  # Mirrors 02c's build_burned_raster(): max() across years == logical OR,
  # and a NULL (missing year) is skipped rather than propagating NA.
  r1 <- rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 1, crs = "EPSG:32610")
  values(r1) <- c(1, 0, 0)
  r2 <- rast(nrows = 1, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 1, crs = "EPSG:32610")
  values(r2) <- c(0, 1, 0)

  build_burned_raster <- function(rasters) {
    acc <- NULL
    for (r_bin in rasters) {
      if (is.null(r_bin)) next
      acc <- if (is.null(acc)) r_bin else max(acc, r_bin, na.rm = TRUE)
    }
    acc
  }

  combined <- build_burned_raster(list(r1, NULL, r2))
  vals <- values(combined)[, 1]
  expect_equal(vals, c(1, 1, 0))
})

test_that("extract_batched-style chunking gives the same result as one unbatched call", {
  r <- rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10, crs = "EPSG:32610")
  values(r) <- 1:100
  set.seed(1)
  pts <- vect(data.frame(x = runif(23, 0.5, 9.5), y = runif(23, 0.5, 9.5)),
              geom = c("x", "y"), crs = "EPSG:32610")

  extract_batched <- function(r, pts, batch_size) {
    n <- length(pts)
    out <- vector("numeric", n)
    for (start in seq(1, n, by = batch_size)) {
      end <- min(start + batch_size - 1, n)
      batch_buf <- terra::buffer(pts[start:end], width = 250)
      vals <- terra::extract(r, batch_buf, fun = max)[, 2]
      out[start:end] <- vals
    }
    out
  }

  unbatched <- extract_batched(r, pts, batch_size = 23)  # one batch = all points
  batched   <- extract_batched(r, pts, batch_size = 5)   # forces 5 batches

  expect_equal(batched, unbatched)
})

test_that("NA from no-overlap extraction is treated as FALSE (not-burned), matching 02c's edge-case handling", {
  pixel_meta <- data.frame(
    pixel_id = 1:3,
    burned_baseline_00_14 = c(TRUE, NA, FALSE),
    burned_2015 = c(NA, FALSE, TRUE)
  )
  all_cols <- c("burned_baseline_00_14", "burned_2015")
  for (col in all_cols) pixel_meta[[col]][is.na(pixel_meta[[col]])] <- FALSE

  expect_equal(pixel_meta$burned_baseline_00_14, c(TRUE, FALSE, FALSE))
  expect_equal(pixel_meta$burned_2015, c(FALSE, FALSE, TRUE))
  expect_false(anyNA(pixel_meta[all_cols]))
})
