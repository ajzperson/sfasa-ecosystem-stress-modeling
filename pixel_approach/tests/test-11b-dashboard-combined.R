library(testthat)
library(tidyverse)
library(sf)
library(patchwork)

source("../src/clustering_helpers.R")

# Snapshot of .GlobalEnv before this file writes anything into it (see
# cleanup at the bottom -- same reasoning as test-06-clustering.R).
.pre_existing_globals <- ls(.GlobalEnv)

# Loads lon_labels/lat_labels/build_map_panel (plus their K/cluster_cols/
# clust_labels/boundary_xy dependencies) from 11b_dashboard_combined.R,
# stopping before build_combined_dashboard() is called for real (which
# reads fire/beetle panels from hardcoded repo-root-relative paths, not
# parameterized -- same reason 11_dashboard.R's build_dashboard() isn't
# covered either; see test-11-dashboard.R).
exprs    <- parse("../src/11b_dashboard_combined.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  "^source\\(|st_read\\(|^clusters_path <-|^pixels <-|^K <-|^boundary_xy <-",
  src_text
)
stop_at <- which(grepl("^build_combined_dashboard\\(", src_text))[1]
stopifnot(!is.na(stop_at))

# build_map_panel()'s closure is forced into .GlobalEnv below, so K/
# boundary_xy must be assigned there explicitly -- plain `<-` would land
# in this file's own environment instead (see test-11-dashboard.R for the
# same issue with boundary_xy).
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
