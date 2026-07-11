library(terra)
library(sf)
library(tidyverse)
library(lubridate)

source("pixel_approach/src/extract_pixel_matrix.R")

# -----------------------------------------------------------------------
# SAMPLE_ONLY: set TRUE to print a few rows for verification, FALSE for
# the full 251,971-pixel run that writes pixel_metadata.csv +
# ndvi_matrix.rds
# -----------------------------------------------------------------------
if (!exists("SAMPLE_ONLY")) SAMPLE_ONLY <- FALSE
if (!exists("N_SAMPLE"))   N_SAMPLE    <- 10

# --- Load and prepare NDVI raster --------------------------------------
nc_file <- "data_common/raw/ndvi/MOD13A1.061_500m_aid0001.nc"
cat("Loading NDVI raster...\n")
r           <- rast(nc_file, subds = "_500m_16_days_NDVI")
dates_16day <- as.Date(time(r))
cat(sprintf("✓ Loaded %d 16-day layers (%s to %s)\n",
            nlyr(r), min(dates_16day), max(dates_16day)))

# Trim to 2000-2016
keep        <- dates_16day <= as.Date("2016-12-31")
r           <- r[[which(keep)]]
dates_16day <- dates_16day[keep]
cat(sprintf("✓ Trimmed to 2000-2016: %d 16-day layers remaining\n", nlyr(r)))

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

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
