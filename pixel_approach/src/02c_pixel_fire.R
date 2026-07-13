library(tidyverse)
library(terra)
library(sf)

# Extracts MTBS (Monitoring Trends in Burn Severity) burn exposure for every
# pixel in the shared grid (pixel_metadata.csv) -- independent of clustering,
# same role as 02a_pixel_vpd.R / 02b_pixel_ppt.R (a pixel-level covariate
# extraction step), not a validation step. 09_fire_validation.R consumes
# this output to cross-check it against cluster assignments.
#
# A pixel counts as "burned" only at moderate (3) or high (4) MTBS severity
# -- low-severity understory fire (1/2) and post-fire regrowth (5) don't
# count, since they likely don't meaningfully disturb canopy NDVI/NDWI.
#
# Validation years (2015-2019) are extracted individually, one column per
# year, rather than pre-collapsed into a single OR'd flag -- so downstream
# validation scripts (09_fire_validation.R) can pick any year range within
# 2015-2019 (e.g. 2015-2016, to match beetle's ADS window, whose
# LEGACY_TPA field is only populated for 2015-2016) without re-running
# this extraction. Baseline (2000-2014) stays a single OR'd flag since
# it's only ever used as an all-or-nothing exclusion criterion, never
# compared across sub-windows.
#
# Output: pixel_approach/data/processed/pixel_fire.csv
#   (pixel_id, burned_baseline_00_14, burned_2015, burned_2016, burned_2017,
#    burned_2018, burned_2019)

sf_use_s2(FALSE)

MTBS_DIR <- "data_common/raw/mtbs/tif"
BURNED_SEVERITY <- c(3, 4)  # moderate, high
CACHE_DIR <- "pixel_approach/data/processed"

# --------------------------------------------------------------------------
# 1. Build a common Sierra-Nevada-sized template grid (MTBS's own CRS/30m),
#    since each year's mosaic is cropped to that year's own fire extent and
#    the extents don't match across years.
# --------------------------------------------------------------------------
sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

sample_r <- rast(file.path(MTBS_DIR, "mtbs_CA_2015.tif"))
sn_albers <- st_transform(sn_boundary, crs(sample_r))
template <- rast(vect(sn_albers), resolution = 30, crs = crs(sample_r))

# --------------------------------------------------------------------------
# 2. Load each year, reclassify to binary (moderate/high severity = 1, else
#    0), and align to the common template. Baseline years (2000-2014) are
#    OR'd together into one raster; validation years (2015-2019) are kept
#    as separate per-year rasters (see header for why).
# --------------------------------------------------------------------------
VALIDATION_YEARS <- 2015:2019

build_year_raster <- function(yr) {
  f <- file.path(MTBS_DIR, sprintf("mtbs_CA_%d.tif", yr))
  if (!file.exists(f)) return(NULL)
  r <- rast(f)
  r_bin <- app(r, function(v) as.numeric(v %in% BURNED_SEVERITY))
  r_bin <- crop(r_bin, template, extend = TRUE)
  r_bin <- resample(r_bin, template, method = "max")
  r_bin[is.na(r_bin)] <- 0
  r_bin
}

build_burned_raster <- function(years) {
  acc <- NULL
  for (yr in years) {
    cat(sprintf("  Processing %d...\n", yr))
    r_bin <- build_year_raster(yr)
    if (is.null(r_bin)) next
    acc <- if (is.null(acc)) r_bin else max(acc, r_bin, na.rm = TRUE)
  }
  acc
}

baseline_cache <- file.path(CACHE_DIR, "mtbs_baseline_severity.tif")
if (file.exists(baseline_cache)) {
  cat("Loading cached baseline severity raster...\n")
  baseline_r <- rast(baseline_cache)
} else {
  cat("Building baseline (2000-2014) burned-severity raster...\n")
  baseline_r <- build_burned_raster(2000:2014)
  writeRaster(baseline_r, baseline_cache, overwrite = TRUE)
}

year_caches <- setNames(
  file.path(CACHE_DIR, sprintf("mtbs_severity_%d.tif", VALIDATION_YEARS)),
  as.character(VALIDATION_YEARS)
)
year_rasters <- list()
for (yr in as.character(VALIDATION_YEARS)) {
  if (file.exists(year_caches[[yr]])) {
    cat(sprintf("Loading cached %s severity raster...\n", yr))
    year_rasters[[yr]] <- rast(year_caches[[yr]])
  } else {
    cat(sprintf("Building %s burned-severity raster...\n", yr))
    r_bin <- build_year_raster(as.integer(yr))
    year_rasters[[yr]] <- r_bin
    if (!is.null(r_bin)) writeRaster(r_bin, year_caches[[yr]], overwrite = TRUE)
  }
}

# --------------------------------------------------------------------------
# 3. Extract burned status at every pixel in the shared grid. Each 500m
#    pixel covers roughly 17x17 of these 30m MTBS cells -- "any" of them at
#    moderate+ severity flags the whole pixel, via a 250m buffer around its
#    centroid (explicit buffer polygons, not extract()'s buffer= shorthand,
#    which was verified to return the buffer width itself instead of
#    aggregating raster values -- a real bug in how that argument behaves
#    combined with fun= in this terra version).
#
#    Processed in batches: extract() over all ~252k points in a single call
#    silently truncates partway through (verified -- the first ~70k rows
#    came back with real values, everything after was NA, exactly matching
#    where a cluster-sorted file's first block ends), so chunking avoids
#    whatever internal limit causes that.
# --------------------------------------------------------------------------
pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv",
                        show_col_types = FALSE)
cat(sprintf("✓ %d pixels in shared grid\n", nrow(pixel_meta)))

pts <- st_as_sf(pixel_meta, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs(baseline_r)) |>
  vect()

extract_batched <- function(r, pts, batch_size = 5000) {
  n <- length(pts)
  out <- vector("numeric", n)
  for (start in seq(1, n, by = batch_size)) {
    end <- min(start + batch_size - 1, n)
    cat(sprintf("  Extracting %d-%d / %d...\n", start, end, n))
    batch_buf <- terra::buffer(pts[start:end], width = 250)
    vals <- terra::extract(r, batch_buf, fun = max)[, 2]
    out[start:end] <- vals
  }
  out
}

cat("Extracting baseline (2000-2014) burn exposure...\n")
baseline_vals <- extract_batched(baseline_r, pts)
pixel_meta$burned_baseline_00_14 <- baseline_vals > 0

year_cols <- character(0)
for (yr in names(year_rasters)) {
  if (is.null(year_rasters[[yr]])) next
  cat(sprintf("Extracting %s burn exposure...\n", yr))
  vals <- extract_batched(year_rasters[[yr]], pts)
  col <- sprintf("burned_%s", yr)
  pixel_meta[[col]] <- vals > 0
  year_cols <- c(year_cols, col)
}

# A handful of pixels sit right at the template's edge, where the 250m
# buffer never overlaps a raster cell at all -- extract() returns NA for
# those, not FALSE. Treated as "not burned" (a reasonable default for a
# tiny edge case) rather than propagating NA through the rest of the
# pipeline, where mean()/sum() would silently turn into NA for an entire
# cluster group from just one NA row.
all_cols <- c("burned_baseline_00_14", year_cols)
n_edge_na <- sum(rowSums(is.na(pixel_meta[all_cols])) > 0)
cat(sprintf("  %d pixels at the template edge had no raster overlap -- treated as not burned\n", n_edge_na))
for (col in all_cols) pixel_meta[[col]][is.na(pixel_meta[[col]])] <- FALSE

cat(sprintf("\nPixels burned (moderate+) in baseline window (2000-2014): %d / %d (%.1f%%)\n",
            sum(pixel_meta$burned_baseline_00_14), nrow(pixel_meta),
            100 * mean(pixel_meta$burned_baseline_00_14)))
for (col in year_cols) {
  cat(sprintf("Pixels burned (moderate+) in %s: %d / %d (%.1f%%)\n",
              col, sum(pixel_meta[[col]]), nrow(pixel_meta), 100 * mean(pixel_meta[[col]])))
}

out_path <- file.path(CACHE_DIR, "pixel_fire.csv")
write_csv(pixel_meta |> select(pixel_id, all_of(all_cols)), out_path)
cat(sprintf("✓ Saved: %s\n", out_path))
