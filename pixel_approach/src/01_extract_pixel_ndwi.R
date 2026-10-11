library(terra)
library(sf)
library(tidyverse)
library(lubridate)

source("pixel_approach/src/extract_pixel_matrix.R")

# Builds the monthly per-pixel NDWI matrix (2000-2016) from MODIS MOD09A1
# bands 2 and 5 (Gao NDWI), on the same pixel grid as NDVI.
#
# Output: pixel_approach/data/processed/ndwi_matrix.rds

# SAMPLE_ONLY = TRUE previews a few pixels instead of the full run
if (!exists("SAMPLE_ONLY")) SAMPLE_ONLY <- FALSE
if (!exists("N_SAMPLE"))   N_SAMPLE    <- 10

nc_file <- "data_common/raw/ndwi/MOD09A1.061_500m_aid0001.nc"
cat("Loading MOD09A1 bands (b02, b05, QC)...\n")
r_b02 <- rast(nc_file, subds = "sur_refl_b02")
r_b05 <- rast(nc_file, subds = "sur_refl_b05")
r_qc  <- rast(nc_file, subds = "sur_refl_qc_500m")

dates_8day <- as.Date(time(r_b02))
stopifnot(all(dates_8day == as.Date(time(r_b05))),
          all(dates_8day == as.Date(time(r_qc))))
cat(sprintf("✓ Loaded %d 8-day layers (%s to %s)\n",
            nlyr(r_b02), min(dates_8day), max(dates_8day)))

keep      <- dates_8day <= as.Date("2016-12-31")
r_b02      <- r_b02[[which(keep)]]
r_b05      <- r_b05[[which(keep)]]
r_qc       <- r_qc[[which(keep)]]
dates_8day <- dates_8day[keep]
cat(sprintf("✓ Trimmed to 2000-2016: %d 8-day layers remaining\n", nlyr(r_b02)))

# Keep only MODLAND QA bits 0-1 == 00 ("ideal quality")
cat("Applying QC mask (MODLAND QA bits 0-1 == ideal quality)...\n")
good_qc <- r_qc %% 4 == 0
r_b02   <- mask(r_b02, good_qc, maskvalue = FALSE)
r_b05   <- mask(r_b05, good_qc, maskvalue = FALSE)

cat("Computing NDWI = (b02 - b05) / (b02 + b05)...\n")
r_ndwi <- (r_b02 - r_b05) / (r_b02 + r_b05)

sn_boundary <- st_read("data_common/processed/sierra_nevada_boundary.geojson", quiet = TRUE)

extract_pixel_matrix(
  r                     = r_ndwi,
  dates                 = dates_8day,
  sn_boundary           = sn_boundary,
  matrix_path           = "pixel_approach/data/processed/ndwi_matrix.rds",
  verify_metadata_path  = "pixel_approach/data/processed/pixel_metadata.csv",
  matrix_label          = "NDWI",
  sample_only           = SAMPLE_ONLY,
  n_sample              = N_SAMPLE
)
