library(terra)
library(sf)
library(tidyverse)
library(lubridate)

source("pixel_approach/src/extract_pixel_matrix.R")

# Extracts a per-pixel monthly NDWI time series (2000-2016, QC-masked) for
# every pixel within the Sierra Nevada boundary, from MODIS MOD09A1 8-day
# bands aggregated to monthly means.
#
# Output: pixel_approach/data/processed/ndwi_matrix.rds

# -----------------------------------------------------------------------
# SAMPLE_ONLY: set TRUE to print a few rows for verification, FALSE for
# the full run that writes ndwi_matrix.rds (verified against the existing
# pixel_metadata.csv, not rewritten)
# -----------------------------------------------------------------------
if (!exists("SAMPLE_ONLY")) SAMPLE_ONLY <- FALSE
if (!exists("N_SAMPLE"))   N_SAMPLE    <- 10

# --- Load MOD09A1 bands + QC --------------------------------------------
# NDWI (Gao) = (NIR - SWIR) / (NIR + SWIR), using MODIS band 2 (NIR,
# 841-876nm) and band 5 (SWIR, 1230-1250nm) — tracks canopy water
# content, a complementary drought-stress signal to NDVI's greenness.
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

# Trim to 2000-2016
keep       <- dates_8day <= as.Date("2016-12-31")
r_b02      <- r_b02[[which(keep)]]
r_b05      <- r_b05[[which(keep)]]
r_qc       <- r_qc[[which(keep)]]
dates_8day <- dates_8day[keep]
cat(sprintf("✓ Trimmed to 2000-2016: %d 8-day layers remaining\n", nlyr(r_b02)))

# --- QC mask: keep only "ideal quality" pixels (MODLAND QA bits 0-1 == 00) ---
cat("Applying QC mask (MODLAND QA bits 0-1 == ideal quality)...\n")
good_qc <- r_qc %% 4 == 0
r_b02   <- mask(r_b02, good_qc, maskvalue = FALSE)
r_b05   <- mask(r_b05, good_qc, maskvalue = FALSE)

# --- Compute NDWI --------------------------------------------------------
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
