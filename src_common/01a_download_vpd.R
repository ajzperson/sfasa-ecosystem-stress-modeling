library(climateR)
library(terra)
library(sf)
library(tidyverse)

cat("Downloading TerraClimate VPD for Sierra Nevada area...\n")
sn_aoi <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

# Download VPD from TerraClimate
cat("Downloading VPD (2000-2025) — this may take a few minutes...\n")
vpd_raw <- getTerraClim(
  AOI       = sn_aoi,
  varname   = "vpd",
  startDate = "2000-01-01",
  endDate   = "2025-12-31"
)

r <- vpd_raw$vpd
cat(sprintf("✓ Downloaded: %d monthly layers\n", nlyr(r)))
cat(sprintf("  Date range: %s to %s\n",
            min(as.Date(terra::time(r))), max(as.Date(terra::time(r)))))

# Save to disk
out_file <- "data_common/raw/vpd/TerraClimate_vpd_2000-2025.nc"
writeCDF(r, out_file, varname = "vpd", unit = "kPa", overwrite = TRUE)
cat(sprintf("✓ Saved: %s\n", out_file))
