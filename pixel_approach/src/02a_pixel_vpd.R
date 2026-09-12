library(tidyverse)
library(terra)

# Extracts per-pixel VPD (Feb 2000 - Dec 2016, monthly) at each NDVI pixel
# centroid, bilinear-interpolated from TerraClimate's ~4km grid down to the
# 500m NDVI grid. Row order matches pixel_metadata.csv (pixel_id).
#
# Output: pixel_approach/data/processed/vpd_matrix.rds

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
# Read off the subsetted raster, not re-derived from target_dates, so column
# labels always match the layers actually extracted below.
vpd_dates <- time(vpd_rast)
cat(sprintf("✓ VPD layers: %d (%s to %s)\n", nlyr(vpd_rast),
            format(vpd_dates[1]),
            format(vpd_dates[nlyr(vpd_rast)])))

vpd_rast <- crop(vpd_rast, ext(pts) + 0.5)
cat("✓ Cropped to Sierra Nevada extent\n")

cat("Extracting VPD at pixel centroids (bilinear interpolation)...\n")
vpd_vals <- terra::extract(vpd_rast, pts, method = "bilinear", ID = FALSE)
colnames(vpd_vals) <- format(vpd_dates, "%Y-%m-%d")

vpd_mat <- as.matrix(vpd_vals)
cat(sprintf("✓ VPD range: %.3f to %.3f kPa | NA count: %d\n",
            min(vpd_mat, na.rm = TRUE), max(vpd_mat, na.rm = TRUE),
            sum(is.na(vpd_mat))))
cat(sprintf("✓ Dimensions: %d pixels × %d months\n", nrow(vpd_mat), ncol(vpd_mat)))

saveRDS(vpd_mat, "pixel_approach/data/processed/vpd_matrix.rds")
cat("✓ Saved: pixel_approach/data/processed/vpd_matrix.rds\n")
