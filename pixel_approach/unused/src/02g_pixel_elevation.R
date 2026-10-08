library(tidyverse)
library(elevatr)
library(sf)

# Extracts elevation at every pixel centroid in the shared grid
# (pixel_metadata.csv) -- independent of clustering, same role as
# 02a_pixel_vpd.R / 02b_pixel_ppt.R / etc. Elevation doesn't vary over time,
# so unlike the climate variables this is a single one-time pull, not a
# time-series matrix.
#
# Output: pixel_approach/data/processed/pixel_elevation.csv
#   (pixel_id, elevation_m)

pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv",
                       show_col_types = FALSE)
cat(sprintf("✓ %d pixel centroids\n", nrow(pixel_meta)))

pts <- st_as_sf(pixel_meta, coords = c("lon", "lat"), crs = 4326)

cat("Fetching elevation (USGS 3DEP via elevatr)...\n")
elev <- get_elev_point(pts, prj = 4326, src = "aws")

pixel_elev <- tibble(
  pixel_id    = pixel_meta$pixel_id,
  elevation_m = elev$elevation
)

cat(sprintf("✓ Elevation range: %.0f to %.0f m | NA count: %d\n",
            min(pixel_elev$elevation_m, na.rm = TRUE),
            max(pixel_elev$elevation_m, na.rm = TRUE),
            sum(is.na(pixel_elev$elevation_m))))

write_csv(pixel_elev, "pixel_approach/data/processed/pixel_elevation.csv")
cat("✓ Saved: pixel_approach/data/processed/pixel_elevation.csv\n")
