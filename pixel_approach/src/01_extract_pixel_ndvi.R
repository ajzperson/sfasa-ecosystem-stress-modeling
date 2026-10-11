library(terra)
library(sf)
library(tidyverse)
library(lubridate)

source("pixel_approach/src/extract_pixel_matrix.R")

# Builds the monthly per-pixel NDVI matrix (2000-2016) from MODIS MOD13A1.
#
# Output: pixel_approach/data/processed/pixel_metadata.csv
#         pixel_approach/data/processed/ndvi_matrix.rds

# SAMPLE_ONLY = TRUE previews a few pixels instead of the full run
if (!exists("SAMPLE_ONLY")) SAMPLE_ONLY <- FALSE
if (!exists("N_SAMPLE"))   N_SAMPLE    <- 10

nc_file <- "data_common/raw/ndvi/MOD13A1.061_500m_aid0001.nc"
cat("Loading NDVI + VI_Quality rasters...\n")
r           <- rast(nc_file, subds = "_500m_16_days_NDVI")
r_qc        <- rast(nc_file, subds = "_500m_16_days_VI_Quality")
dates_16day <- as.Date(time(r))
stopifnot(all(dates_16day == as.Date(time(r_qc))))
cat(sprintf("✓ Loaded %d 16-day layers (%s to %s)\n",
            nlyr(r), min(dates_16day), max(dates_16day)))

keep       <- dates_16day <= as.Date("2016-12-31")
r           <- r[[which(keep)]]
r_qc        <- r_qc[[which(keep)]]
dates_16day <- dates_16day[keep]
cat(sprintf("✓ Trimmed to 2000-2016: %d 16-day layers remaining\n", nlyr(r)))

# Drop only VI_Quality bits 0-1 == 11 ("not produced")
cat("Applying QC mask (VI_Quality bits 0-1 <= 10: excluding not-produced only)...\n")
good_qc <- r_qc %% 4 <= 2
r       <- mask(r, good_qc, maskvalue = FALSE)

sn_boundary <- st_read("data_common/processed/sierra_nevada_boundary.geojson", quiet = TRUE)

extract_pixel_matrix(
  r             = r,
  dates         = dates_16day,
  sn_boundary   = sn_boundary,
  metadata_path = "pixel_approach/data/processed/pixel_metadata.csv",
  matrix_path   = "pixel_approach/data/processed/ndvi_matrix.rds",
  matrix_label  = "NDVI",
  sample_only   = SAMPLE_ONLY,
  n_sample      = N_SAMPLE
)
