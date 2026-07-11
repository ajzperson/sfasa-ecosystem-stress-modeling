library(tidyverse)
library(terra)

# Extracts per-pixel VPD at each NDVI pixel centroid using bilinear
# interpolation (downscales TerraClimate ~4km to 500m NDVI grid).
# Crop to Sierra Nevada extent before extraction for speed.
#
# Output:
#   pixel_approach/data/processed/vpd_matrix.rds — 251971 × 203 months
# Row order matches pixel_metadata.csv (indexed by pixel_id).
# Column names: "YYYY-MM-DD", Feb 2000 – Dec 2016 (extended from Dec 2014 to
# cover the 2015-2016 validation window).

target_dates <- seq(as.Date("2000-02-01"), as.Date("2016-12-01"), by = "month")
target_ym    <- format(target_dates, "%Y%m")

cat("Loading pixel metadata...\n")
pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv",
                       show_col_types = FALSE)
pts <- vect(pixel_meta, geom = c("lon", "lat"), crs = "EPSG:4326")
cat(sprintf("✓ %d pixel centroids\n", nrow(pts)))

cat("Loading VPD (TerraClimate)...\n")
vpd_rast <- rast("data_common/raw/vpd/TerraClimate_vpd_2000-2025.nc")
crs(vpd_rast) <- "EPSG:4326"   # WGS84 — not embedded in the file

vpd_ym   <- format(time(vpd_rast), "%Y%m")
keep_v   <- which(vpd_ym %in% target_ym)
vpd_rast <- vpd_rast[[keep_v]]
cat(sprintf("✓ VPD layers: %d (%s to %s)\n", nlyr(vpd_rast),
            format(time(vpd_rast)[1]),
            format(time(vpd_rast)[nlyr(vpd_rast)])))

# Crop to Sierra Nevada extent before extraction
vpd_rast <- crop(vpd_rast, ext(pts) + 0.5)
cat("✓ Cropped to Sierra Nevada extent\n")

cat("Extracting VPD at pixel centroids (bilinear interpolation)...\n")
vpd_vals <- terra::extract(vpd_rast, pts, method = "bilinear", ID = FALSE)
colnames(vpd_vals) <- format(target_dates[target_ym %in% vpd_ym[keep_v]], "%Y-%m-%d")

vpd_mat <- as.matrix(vpd_vals)
cat(sprintf("✓ VPD range: %.3f to %.3f kPa | NA count: %d\n",
            min(vpd_mat, na.rm = TRUE), max(vpd_mat, na.rm = TRUE),
            sum(is.na(vpd_mat))))
cat(sprintf("✓ Dimensions: %d pixels × %d months\n", nrow(vpd_mat), ncol(vpd_mat)))

saveRDS(vpd_mat, "pixel_approach/data/processed/vpd_matrix.rds")
cat("✓ Saved: pixel_approach/data/processed/vpd_matrix.rds\n")
