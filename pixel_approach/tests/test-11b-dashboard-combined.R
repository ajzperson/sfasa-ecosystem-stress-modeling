library(testthat)
library(tidyverse)
library(sf)
library(patchwork)

source("../src/clustering_helpers.R")

# Tests for 11b_dashboard_combined.R (map panel only).

# Globals present before this file; anything newer is removed at the end
.pre_existing_globals <- ls(.GlobalEnv)

# build_combined_dashboard() reads fixed panel paths, so it isn't tested
exprs    <- parse("../src/11b_dashboard_combined.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|st_read\\(|^clusters_path <-|^pixels <-|^K <-|^boundary_xy <-",
  src_text
)
stop_at <- which(grepl("^build_combined_dashboard\\(", src_text))[1]
stopifnot(!is.na(stop_at))

# Assigned into .GlobalEnv, where build_map_panel() looks them up
assign("K", 4, envir = .GlobalEnv)
assign("boundary_xy", data.frame(X = c(-121, -121, -118, -118, -121),
                                  Y = c(35, 40, 40, 35, 35)),
       envir = .GlobalEnv)

for (i in seq_len(stop_at - 1)) {
  if (!skip_lines[i]) eval(exprs[[i]], envir = .GlobalEnv)
}

test_that("lon_labels formats west/east/prime-meridian correctly", {
  expect_equal(lon_labels(-120), "120°W")
  expect_equal(lon_labels(120), "120°E")
  expect_equal(lon_labels(0), "0°")
})

test_that("lat_labels formats south/north/equator correctly", {
  expect_equal(lat_labels(-38), "38°S")
  expect_equal(lat_labels(38), "38°N")
})

test_that("build_map_panel returns a ggplot colored by all K=4 clusters, not a hardcoded 3", {
  d <- tibble(pixel_id = 1:8,
              lon = c(-120, -119.5, -119, -120, -119.5, -119, -120, -119),
              lat = c(36, 37, 38, 36, 37, 38, 36, 37),
              cluster = c(1, 1, 2, 2, 3, 3, 4, 4))

  p <- build_map_panel(d)
  expect_s3_class(p, "ggplot")

  built_colors <- unique(ggplot_build(p)$data[[1]]$colour)
  expect_equal(length(built_colors), 4)
})

rm(list = setdiff(ls(.GlobalEnv), .pre_existing_globals), envir = .GlobalEnv)
