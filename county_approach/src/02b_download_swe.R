library(climateR)
library(terra)
library(tidyverse)

cat("Downloading TerraClimate SWE for Sierra Nevada area...\n")
source("county_approach/src/cluster_polygons.R")
sn_aoi <- st_union(c(sequoia_poly, plumas_poly))

cat("Downloading SWE (2000-2025) — this may take a few minutes...\n")
swe_raw <- getTerraClim(
  AOI       = sn_aoi,
  varname   = "swe",
  startDate = "2000-01-01",
  endDate   = "2025-12-31"
)

r <- swe_raw$swe
cat(sprintf("✓ Downloaded: %d monthly layers\n", nlyr(r)))
cat(sprintf("  Date range: %s to %s\n",
            min(as.Date(terra::time(r))), max(as.Date(terra::time(r)))))

out_file <- "data_common/raw/swe/TerraClimate_swe_2000-2025.nc"
writeCDF(r, out_file, varname = "swe", unit = "mm", overwrite = TRUE)
cat(sprintf("✓ Saved: %s\n", out_file))
