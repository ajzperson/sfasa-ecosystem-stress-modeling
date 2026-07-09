library(terra)
library(tidyverse)

cat("Processing TerraClimate VPD for Sequoia and Plumas clusters...\n")
write_cluster_files <- FALSE
source("county_approach/src/cluster_polygons.R")

# Load raw VPD raster
raw_file <- "data_common/raw/vpd/TerraClimate_vpd_2000-2025.nc"
if (!file.exists(raw_file)) {
  stop("Raw VPD file not found. Run src_common/01a_download_vpd.R first.")
}
r <- rast(raw_file)
cat(sprintf("✓ Loaded: %d monthly layers\n", nlyr(r)))
cat(sprintf(
  "  VPD raw range: %.3f to %.3f\n",
  global(r, min, na.rm = TRUE)[[1]],
  global(r, max, na.rm = TRUE)[[1]]
))

# Extract mean VPD for each cluster polygon
clusters_v <- vect(st_transform(clusters_sf, crs(r)))
vals       <- terra::extract(r, clusters_v, fun = mean, na.rm = TRUE)

# Parse dates from NetCDF time dimension
dates <- as.Date(terra::time(r))

vpd_df <- bind_rows(
  data.frame(
    date    = dates,
    cluster = "Sequoia",
    vpd     = round(as.numeric(vals[1, -1]), 3)
  ),
  data.frame(
    date    = dates,
    cluster = "Plumas",
    vpd     = round(as.numeric(vals[2, -1]), 3)
  )
) |>
  arrange(cluster, date)

cat(sprintf(
  "\n  Records: %d (%d per cluster)\n",
  nrow(vpd_df), as.integer(nrow(vpd_df) / 2)
))
cat(sprintf(
  "  Date range: %s to %s\n", min(vpd_df$date), max(vpd_df$date)
))
cat(sprintf(
  "  VPD range — Sequoia: %.3f - %.3f kPa\n",
  min(vpd_df$vpd[vpd_df$cluster == "Sequoia"]),
  max(vpd_df$vpd[vpd_df$cluster == "Sequoia"])
))
cat(sprintf(
  "  VPD range — Plumas:  %.3f - %.3f kPa\n",
  min(vpd_df$vpd[vpd_df$cluster == "Plumas"]),
  max(vpd_df$vpd[vpd_df$cluster == "Plumas"])
))

write_csv(vpd_df, "county_approach/data/processed/Sierra_Nevada_VPD_2000-2025.csv")
cat("\n✓ Saved: county_approach/data/processed/Sierra_Nevada_VPD_2000-2025.csv\n")
cat("  Columns: date, cluster, vpd (kPa)\n")
cat("  Spatial method: county polygons ∩ Sierra Nevada boundary\n")
