library(tidyverse)
library(terra)
library(sf)

# Flags each pixel burned (MTBS moderate/high severity) for 2000-2014
# combined and for each year 2015-2019.
#
# Output: pixel_approach/data/processed/pixel_fire.csv
#         pixel_approach/data/processed/mtbs_*severity*.tif (caches)

sf_use_s2(FALSE)

MTBS_ZIP <- "data_common/raw/mtbs/bsp_ugjds3epryqgif0ej5sb.zip"
MTBS_DIR <- "data_common/raw/mtbs/tif"
BURNED_SEVERITY <- c(3, 4)  # moderate, high
CACHE_DIR <- "pixel_approach/data/processed"

# Outer zip holds one zip per year
if (!dir.exists(MTBS_DIR) || length(list.files(MTBS_DIR, pattern = "\\.tif$")) == 0) {
  cat("Unzipping MTBS data (two levels: outer bundle -> per-year zips)...\n")
  tmp_dir <- tempfile("mtbs_extract_")
  dir.create(tmp_dir)
  unzip(MTBS_ZIP, exdir = tmp_dir)
  inner_zips <- list.files(tmp_dir, pattern = "\\.zip$", full.names = TRUE)
  dir.create(MTBS_DIR, showWarnings = FALSE, recursive = TRUE)
  for (z in inner_zips) unzip(z, exdir = MTBS_DIR)
  unlink(tmp_dir, recursive = TRUE)
  cat(sprintf("✓ Extracted %d years of MTBS data to %s\n", length(inner_zips), MTBS_DIR))
}

# Common 30m grid, since each year's mosaic has a different extent
sn_boundary <- st_read("data_common/processed/sierra_nevada_boundary.geojson", quiet = TRUE)

sample_r <- rast(file.path(MTBS_DIR, "mtbs_CA_2015.tif"))
sn_albers <- st_transform(sn_boundary, crs(sample_r))
template <- rast(vect(sn_albers), resolution = 30, crs = crs(sample_r))

# Burned (1/0) rasters: 2000-2014 OR'd into one, 2015-2019 one per year
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

# A pixel is burned if any 30m cell within 250m of its centroid is burned.
# Batched because one extract() over all points silently returns NAs.
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

# Edge pixels with no raster overlap (NA) count as not burned
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
