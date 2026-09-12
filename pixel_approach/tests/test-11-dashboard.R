library(testthat)
library(tidyverse)
library(sf)
library(patchwork)

source("../src/clustering_helpers.R")

# Snapshot of .GlobalEnv before this file writes anything into it (see
# cleanup at the bottom -- same reasoning as test-06-clustering.R).
.pre_existing_globals <- ls(.GlobalEnv)

# Loads lon_labels/lat_labels/build_map_panel (plus their cluster_cols/
# clust_labels/boundary_xy dependencies) from 11_dashboard.R, stopping
# before the driver code that loads real ~250k-pixel data and calls
# build_dashboard() (which reads fire/beetle panels from hardcoded
# repo-root-relative paths -- pixel_approach/data/processed/panel_fire_
# <idx><suffix>.rds -- not parameterized via its `cfg` argument, so it
# can't be pointed at synthetic fixtures the way build_map_panel can; not
# covered here). sn_boundary/boundary_xy are skipped (real geojson) and
# replaced with a tiny synthetic square boundary instead.
exprs    <- parse("../src/11_dashboard.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|st_read\\(|^boundary_xy <-|^index_config <-",
  src_text
)
stop_at <- which(grepl("^for \\(idx in names\\(index_config\\)\\)", src_text))[1]
stopifnot(!is.na(stop_at))

# build_map_panel()'s closure is forced into .GlobalEnv below (it needs
# cluster_cols/clust_labels there too, same reasoning as
# test-06-clustering.R's run_clustering()), so boundary_xy must be
# assigned there explicitly -- plain `<-` would land in this file's own
# environment instead, invisible to a .GlobalEnv-closured function.
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
  expect_equal(lat_labels(0), "0°")
})

test_that("lon_labels/lat_labels drop trailing zeros and decimals like %g", {
  expect_equal(lon_labels(-120.50), "120.5°W")
  expect_equal(lat_labels(38.00), "38°N")
})

test_that("build_map_panel returns a ggplot colored by the pixel's own cluster column", {
  clusters_path <- tempfile(fileext = ".csv")
  write_csv(tibble(pixel_id = 1:6, lon = c(-120, -119.5, -119, -120, -119.5, -119),
                    lat = c(36, 37, 38, 36, 37, 38),
                    cluster = c(1, 1, 2, 2, 3, 3)),
            clusters_path)

  p <- build_map_panel(clusters_path, "TEST")
  expect_s3_class(p, "ggplot")

  built_colors <- unique(ggplot_build(p)$data[[1]]$colour)
  expect_equal(length(built_colors), 3)
})

rm(list = setdiff(ls(.GlobalEnv), .pre_existing_globals), envir = .GlobalEnv)
