library(tidyverse)
library(terra)

# Extracts per-pixel SWE (snow water equivalent) at each NDVI pixel centroid
# using bilinear interpolation (downscales TerraClimate ~4km to 500m NDVI
# grid). Crop to Sierra Nevada extent before extraction for speed. Mirrors
# 02a_pixel_vpd.R's extraction pattern.
#
# Output:
#   pixel_approach/data/processed/swe_matrix.rds — 251971 × 203 months
# Row order matches pixel_metadata.csv (indexed by pixel_id).
# Column names: "YYYY-MM-DD", Feb 2000 – Dec 2016.

target_dates <- seq(as.Date("2000-02-01"), as.Date("2016-12-01"), by = "month")
target_ym    <- format(target_dates, "%Y%m")

cat("Loading pixel metadata...\n")
pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv",
                       show_col_types = FALSE)
pts <- vect(pixel_meta, geom = c("lon", "lat"), crs = "EPSG:4326")
cat(sprintf("✓ %d pixel centroids\n", nrow(pts)))

cat("Loading SWE (TerraClimate)...\n")
swe_rast <- rast("data_common/raw/swe/TerraClimate_swe_2000-2025.nc")
crs(swe_rast) <- "EPSG:4326"   # WGS84 — not embedded in the file

swe_ym   <- format(time(swe_rast), "%Y%m")
keep_s   <- which(swe_ym %in% target_ym)
swe_rast <- swe_rast[[keep_s]]
# Captured directly off the subsetted raster (not re-derived by filtering
# target_dates via %in%) so column labels are guaranteed to match the
# layers actually extracted below, regardless of source ordering —
# matches the pattern used in 02a_pixel_vpd.R.
swe_dates <- time(swe_rast)
cat(sprintf("✓ SWE layers: %d (%s to %s)\n", nlyr(swe_rast),
            format(swe_dates[1]),
            format(swe_dates[nlyr(swe_rast)])))

# Crop to Sierra Nevada extent before extraction
swe_rast <- crop(swe_rast, ext(pts) + 0.5)
cat("✓ Cropped to Sierra Nevada extent\n")

cat("Extracting SWE at pixel centroids (bilinear interpolation)...\n")
swe_vals <- terra::extract(swe_rast, pts, method = "bilinear", ID = FALSE)
colnames(swe_vals) <- format(swe_dates, "%Y-%m-%d")

swe_mat <- as.matrix(swe_vals)
cat(sprintf("✓ SWE range: %.3f to %.3f mm | NA count: %d\n",
            min(swe_mat, na.rm = TRUE), max(swe_mat, na.rm = TRUE),
            sum(is.na(swe_mat))))
cat(sprintf("✓ Dimensions: %d pixels × %d months\n", nrow(swe_mat), ncol(swe_mat)))

saveRDS(swe_mat, "pixel_approach/data/processed/swe_matrix.rds")
cat("✓ Saved: pixel_approach/data/processed/swe_matrix.rds\n")
