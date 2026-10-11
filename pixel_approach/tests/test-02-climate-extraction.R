library(testthat)
library(terra)

# Tests for the logic in 02a_pixel_vpd.R and 02b_pixel_ppt.R, replicated on synthetic data.

target_dates <- seq(as.Date("2000-02-01"), as.Date("2016-12-01"), by = "month")
target_ym    <- format(target_dates, "%Y%m")

test_that("date-range filtering (%in% target_ym) selects the correct subset", {
  source_dates <- as.Date(c("1999-12-01", "2000-01-01", "2000-02-01",
                            "2005-06-01", "2016-12-01", "2017-01-01"))
  source_ym <- format(source_dates, "%Y%m")
  keep <- which(source_ym %in% target_ym)

  expect_equal(source_dates[keep],
               as.Date(c("2000-02-01", "2005-06-01", "2016-12-01")))
})

test_that("column dates derived from time() match layers even when source is out of order (regression test)", {
  # Labels must come from the raster's own dates, even when out of order
  dates_out_of_order <- as.Date(c("2020-03-01", "2020-02-01", "2020-01-01"))
  r <- rast(nrows = 1, ncols = 1, nlyrs = 3, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
            crs = "EPSG:4326")
  values(r) <- matrix(c(300, 200, 100), nrow = 1)
  time(r) <- dates_out_of_order

  keep_v <- which(format(time(r), "%Y%m") %in% format(dates_out_of_order, "%Y%m"))
  r <- r[[keep_v]]
  extracted_dates <- time(r)  # the fix: read dates off the subsetted raster itself

  pts <- vect(data.frame(lon = 0.5, lat = 0.5), geom = c("lon", "lat"), crs = "EPSG:4326")
  vals <- terra::extract(r, pts, ID = FALSE)
  colnames(vals) <- format(extracted_dates, "%Y-%m-%d")

  expect_equal(vals[["2020-01-01"]], 100)
  expect_equal(vals[["2020-02-01"]], 200)
  expect_equal(vals[["2020-03-01"]], 300)
})

test_that("PRISM filename regex extracts YYYYMM correctly", {
  extract_ym <- function(files) sub(".*_([0-9]{6})\\.bil$", "\\1", basename(files))

  expect_equal(
    extract_ym(c(
      "data_common/raw/prism/prismz/prism_ppt_us_25m_200002/prism_ppt_us_25m_200002.bil",
      "prism_ppt_us_25m_201612.bil"
    )),
    c("200002", "201612")
  )
})

test_that("sort-before-stack pattern (02b) produces chronologically ordered dates matching file order", {
  files_unsorted <- c("prism_ppt_us_25m_200003.bil",
                      "prism_ppt_us_25m_200001.bil",
                      "prism_ppt_us_25m_200002.bil")
  ym <- sub(".*_([0-9]{6})\\.bil$", "\\1", basename(files_unsorted))

  ord         <- order(ym)
  files_sorted <- files_unsorted[ord]
  ym_sorted   <- ym[ord]
  dates       <- as.Date(paste0(ym_sorted, "01"), "%Y%m%d")

  expect_equal(files_sorted, c("prism_ppt_us_25m_200001.bil",
                               "prism_ppt_us_25m_200002.bil",
                               "prism_ppt_us_25m_200003.bil"))
  expect_true(!is.unsorted(dates, strictly = TRUE))
  expect_equal(dates, as.Date(c("2000-01-01", "2000-02-01", "2000-03-01")))
})

test_that("bilinear extraction interpolates between known cell values", {
  # Halfway between cells of 10 and 20 should give 15
  r <- rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1, crs = "EPSG:4326")
  values(r) <- c(10, 20)

  pt <- vect(data.frame(lon = 1.0, lat = 0.5), geom = c("lon", "lat"), crs = "EPSG:4326")
  val <- terra::extract(r, pt, method = "bilinear", ID = FALSE)

  expect_equal(val[[1]], 15, tolerance = 1e-6)
})
