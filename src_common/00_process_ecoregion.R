library(sf)

# Filters the EPA Level III ecoregion shapefile to US_L3CODE == "5" (the
# region named "Sierra Nevada" in US_L3NAME), reprojects to EPSG:4326 (the
# GeoJSON spec, RFC 7946, mandates WGS84 -- output CRS reads as OGC:CRS84,
# same datum, GeoJSON's required axis-order label), and saves as GeoJSON.
# Area is computed via EPSG:3310 (NAD83 / California Albers, an equal-area
# projection for California): 53,102.7 km².
#
# Output: docs_common/sierra_nevada_boundary.geojson

sf_use_s2(TRUE)

zip_file  <- "data_common/raw/ecoregion/us_eco_l3.zip"
out_dir   <- "data_common/raw/ecoregion/us_eco_l3"
out_file  <- "docs_common/sierra_nevada_boundary.geojson"

if (!dir.exists(out_dir)) {
  cat("Unzipping EPA Level III ecoregion shapefile...\n")
  unzip(zip_file, exdir = out_dir)
}

shp <- file.path(out_dir, "us_eco_l3.shp")
eco <- st_read(shp, quiet = TRUE)

sierra <- eco |>
  dplyr::filter(US_L3CODE == "5") |>
  st_transform(4326) |>
  st_union()

area_km2 <- round(as.numeric(st_area(st_transform(sierra, 3310))) / 1e6, 1)
cat(sprintf("✓ Sierra Nevada ecoregion area: %s km²\n", format(area_km2, big.mark = ",")))

st_write(st_sf(name = "Sierra Nevada", geometry = sierra), out_file, delete_dsn = TRUE, quiet = TRUE)
cat(sprintf("✓ Saved: %s\n", out_file))
