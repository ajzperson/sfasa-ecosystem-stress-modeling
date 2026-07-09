library(terra)
library(sf)
library(tidyverse)

# Builds a binary "ever burned 2000-2016" mask on the same grid as the NDVI
# cube (06_process_ndvi.R), so it can be applied before the cluster spatial
# mean is computed. A wildfire NDVI crash is a different mechanism than the
# slow drought/beetle-driven biotic decline this pipeline is built to detect
# — leaving burned pixels in would let fire-driven variance (unexplained by
# VPD/ppt) contaminate both the cluster mean and the 2000-2011 climate-
# correction baseline regression in 11_1_csd_manual.R.

nc_file <- "data_common/raw/ndvi/MOD13A1.061_500m_aid0001.nc"
grid <- rast(nc_file, subds = "_500m_16_days_NDVI")[[1]]

fire_file <- "data_common/raw/fire_perimeters/California_Historic_Fire_Perimeters_-4934337877712688794.gpkg"
cat("Reading CAL FIRE FRAP perimeters...\n")
fires_raw <- st_read(fire_file, layer = "California_Fire_Perimeters__all_", quiet = TRUE)

sn_boundary <- st_union(st_transform(
  st_geometry(st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)),
  st_crs(fires_raw)
))

fires <- fires_raw |>
  filter(YEAR_ >= 2000, YEAR_ <= 2016) |>
  st_make_valid() |>
  st_filter(st_as_sf(sn_boundary), .predicate = st_intersects) |>
  st_transform(crs(grid))

cat(sprintf("✓ %d fire perimeters intersecting Sierra Nevada boundary, 2000-2016\n",
            nrow(fires)))

fire_mask <- rasterize(vect(fires), grid, field = 1, background = 0)
# rasterize(..., background = 0) fills the whole rectangular grid extent,
# overwriting the template's original NA cells (outside the valid AppEEARS
# footprint) with 0 — restore that NA structure so pixel counts/diagnostics
# reflect the actual study-area footprint, not the bounding box
fire_mask <- mask(fire_mask, grid)
names(fire_mask) <- "ever_burned_2000_2016"

n_burned <- as.numeric(global(fire_mask == 1, "sum", na.rm = TRUE))
n_total  <- as.numeric(global(!is.na(fire_mask), "sum", na.rm = TRUE))
cat(sprintf("Pixels ever burned 2000-2016 (full NDVI grid extent): %d / %d (%.1f%%)\n",
            n_burned, n_total, 100 * n_burned / n_total))

dir.create("county_approach/data/processed", showWarnings = FALSE, recursive = TRUE)
writeRaster(fire_mask, "county_approach/data/processed/Sierra_Nevada_FireMask_2000-2016.tif",
            overwrite = TRUE)
cat("✓ Saved: county_approach/data/processed/Sierra_Nevada_FireMask_2000-2016.tif\n")

# --- Acres burned by cluster and year (same shape as Sierra_Nevada_BeetleAcres
# and the Kaggle-derived fire validation CSV: cluster, year, acres) ---
write_cluster_files <- FALSE
source("county_approach/src/cluster_polygons.R")

# Assign to cluster using centroid (avoids double-counting boundary polygons)
centroids  <- st_centroid(st_transform(fires, st_crs(clusters_sf)))
in_sequoia <- st_intersects(centroids, sequoia_poly, sparse = FALSE)[, 1]
in_plumas  <- st_intersects(centroids, plumas_poly,  sparse = FALSE)[, 1]

fires_clustered <- fires |>
  mutate(cluster = case_when(
    in_sequoia ~ "Sequoia",
    in_plumas  ~ "Plumas",
    TRUE       ~ NA_character_
  )) |>
  filter(!is.na(cluster))

cat(sprintf("  Sequoia: %d fires | Plumas: %d fires\n",
            sum(fires_clustered$cluster == "Sequoia"),
            sum(fires_clustered$cluster == "Plumas")))

fire_acres_annual <- fires_clustered |>
  st_drop_geometry() |>
  group_by(cluster, year = YEAR_) |>
  summarise(n_fires = n(),
            acres_burned = round(sum(GIS_ACRES, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(cluster, year)

cat("\nFire acres by year:\n")
print(fire_acres_annual |>
        select(-n_fires) |>
        pivot_wider(names_from = cluster, values_from = acres_burned,
                    values_fill = 0))

out_file <- "county_approach/data/processed/Sierra_Nevada_FireAcres_2000-2016.csv"
write_csv(fire_acres_annual, out_file)
cat(sprintf("✓ Saved: %s\n", out_file))
