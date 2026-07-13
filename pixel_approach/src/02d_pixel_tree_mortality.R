library(tidyverse)
library(terra)
library(sf)

# Extracts LCMAP tree-cover-loss transitions for every pixel in the shared
# grid (pixel_metadata.csv) -- independent of clustering, same role as
# 02a_pixel_vpd.R / 02b_pixel_ppt.R / 02c_pixel_fire.R.
#
# Uses LCACHG (Annual Land Cover Change) only, not LCPRI -- LCACHG's own
# two-digit code already encodes both the prior and current year's land
# cover class, so a separate LCPRI lookup isn't needed for this purpose.
# LCACHG value = 1-8 means no change (single-digit class code). A two-digit
# value in 12-87 means change, encoded as (prior class)(current class).
# Tree Cover = class 4, verified against the official USGS LCMAP Data
# Format Control Book (LSDS-2074, Table 2-4). Tree cover loss is therefore
# any LCACHG value in {41, 42, 43, 45, 46, 47, 48} -- prior class 4 (tree),
# current class anything else.
#
# Years 2015-2019 are extracted individually, one column per year, rather
# than pre-collapsed into a single OR'd flag -- so downstream validation
# scripts (10_lcmap_validation.R) can pick any year range within 2015-2019
# (e.g. 2015-2016, to match beetle's ADS window, whose LEGACY_TPA field is
# only populated for 2015-2016) without re-running this extraction. Land
# cover change during the 2000-2014 baseline isn't an unrelated confound
# the way fire is, so there's no baseline-exclusion use case here (see
# DATAPIPELINE.md).
#
# Output: pixel_approach/data/processed/pixel_tree_mortality.csv
#   (pixel_id, lost_tree_cover_2015, lost_tree_cover_2016,
#    lost_tree_cover_2017, lost_tree_cover_2018, lost_tree_cover_2019)

sf_use_s2(FALSE)

LCMAP_DIR <- "data_common/raw/lcmap"
TREE_LOSS_CODES <- c(41, 42, 43, 45, 46, 47, 48)
CACHE_DIR <- "pixel_approach/data/processed"
YEARS <- 2015:2019

# --------------------------------------------------------------------------
# 1. Load each year's LCACHG and reclassify to binary (tree loss = 1, else
#    0), kept as separate per-year rasters (see header for why) -- no
#    alignment step needed, unlike MTBS, since LCMAP's CONUS files share
#    an identical fixed grid every year (verified: same extent/resolution/
#    CRS across years).
# --------------------------------------------------------------------------
cat("Building per-year tree-cover-loss rasters (2015-2019)...\n")
year_rasters <- list()
for (yr in as.character(YEARS)) {
  f <- file.path(LCMAP_DIR, sprintf("LCMAP_CU_%s_V13_LCACHG.tif", yr))
  cat(sprintf("  Processing %s...\n", yr))
  r <- rast(f)
  year_rasters[[yr]] <- app(r, function(v) as.numeric(v %in% TREE_LOSS_CODES))
}

# --------------------------------------------------------------------------
# 2. Extract tree-loss status at every pixel in the shared grid, same
#    250m-buffer aggregation approach as MTBS (each 500m study pixel covers
#    roughly 17x17 of these 30m LCMAP cells -- "any" of them showing tree
#    loss flags the whole pixel), processed in batches for memory safety.
# --------------------------------------------------------------------------
pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv",
                        show_col_types = FALSE)
cat(sprintf("✓ %d pixels in shared grid\n", nrow(pixel_meta)))

pts <- st_as_sf(pixel_meta, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs(year_rasters[[1]])) |>
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

year_cols <- character(0)
for (yr in names(year_rasters)) {
  cat(sprintf("Extracting %s tree-cover-loss status...\n", yr))
  vals <- extract_batched(year_rasters[[yr]], pts)
  col <- sprintf("lost_tree_cover_%s", yr)
  pixel_meta[[col]] <- vals > 0
  year_cols <- c(year_cols, col)
}

# A handful of pixels at the grid edge may have no raster overlap at all --
# treat as "no loss detected" rather than propagating NA (same reasoning as
# 02c_pixel_fire.R).
n_edge_na <- sum(rowSums(is.na(pixel_meta[year_cols])) > 0)
cat(sprintf("  %d pixels at the grid edge had no raster overlap -- treated as no loss\n", n_edge_na))
for (col in year_cols) pixel_meta[[col]][is.na(pixel_meta[[col]])] <- FALSE

for (col in year_cols) {
  cat(sprintf("Pixels with tree cover loss in %s: %d / %d (%.1f%%)\n",
              col, sum(pixel_meta[[col]]), nrow(pixel_meta), 100 * mean(pixel_meta[[col]])))
}

out_path <- file.path(CACHE_DIR, "pixel_tree_mortality.csv")
write_csv(pixel_meta |> select(pixel_id, all_of(year_cols)), out_path)
cat(sprintf("✓ Saved: %s\n", out_path))
