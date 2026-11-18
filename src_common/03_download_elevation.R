library(elevatr)
library(terra)
library(sf)

# Downloads elevation (AWS terrain tiles, USGS 3DEP-based, zoom 9, ~110m) for the
# Sierra Nevada AOI's bounding box, and saves as GeoTIFF.
#
# Output: data_common/raw/elevation/sierra_nevada_elevation_z9.tif

cat("Downloading elevation for Sierra Nevada area...\n")
sn_aoi <- st_read("data_common/processed/sierra_nevada_boundary.geojson", quiet = TRUE)

dem <- rast(get_elev_raster(sn_aoi, z = 9, src = "aws", clip = "bbox"))
cat(sprintf("✓ Downloaded: %d x %d cells, %.4f° resolution\n", nrow(dem), ncol(dem), res(dem)[1]))

dir.create("data_common/raw/elevation", showWarnings = FALSE, recursive = TRUE)
out_file <- "data_common/raw/elevation/sierra_nevada_elevation_z9.tif"
writeRaster(dem, out_file, overwrite = TRUE)
cat(sprintf("✓ Saved: %s\n", out_file))
