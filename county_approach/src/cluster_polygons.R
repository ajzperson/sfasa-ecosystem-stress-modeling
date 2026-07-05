library(sf)
library(tigris)

options(tigris_use_cache = TRUE)
sf_use_s2(TRUE)

hhz_counties     <- c("Amador", "Calaveras", "El Dorado", "Fresno", "Kern",
                      "Madera", "Mariposa", "Placer", "Tulare", "Tuolumne")
non_hhz_counties <- c("Nevada", "Plumas", "Sierra")

ca_counties <- counties(state = "CA", cb = TRUE, year = 2020,
                        progress_bar = FALSE)
sn_boundary <- st_union(st_transform(
  st_geometry(st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)),
  st_crs(ca_counties)
))

sequoia_poly <- st_intersection(
  st_union(ca_counties[ca_counties$NAME %in% hhz_counties, ]),
  sn_boundary
)
plumas_poly <- st_intersection(
  st_union(ca_counties[ca_counties$NAME %in% non_hhz_counties, ]),
  sn_boundary
)

clusters_sf <- rbind(
  st_sf(cluster = "Sequoia", geometry = sequoia_poly),
  st_sf(cluster = "Plumas",  geometry = plumas_poly)
)

cat("✓ Sequoia polygon: 10 HHZ counties ∩ Sierra Nevada boundary\n")
cat("✓ Plumas polygon:  3 non-HHZ counties ∩ Sierra Nevada boundary\n")

if (!exists("write_cluster_files")) write_cluster_files <- TRUE

if (write_cluster_files) {
  st_write(st_sf(cluster = "Sequoia", geometry = sequoia_poly),
           "county_approach/docs/cluster_sequoia.geojson",   delete_dsn = TRUE, quiet = TRUE)
  st_write(st_sf(cluster = "Plumas",  geometry = plumas_poly),
           "county_approach/docs/cluster_plumas.geojson",    delete_dsn = TRUE, quiet = TRUE)
  st_write(clusters_sf,
           "county_approach/docs/clusters_combined.geojson", delete_dsn = TRUE, quiet = TRUE)

  cat("✓ Saved: county_approach/docs/cluster_sequoia.geojson\n")
  cat("✓ Saved: county_approach/docs/cluster_plumas.geojson\n")
  cat("✓ Saved: county_approach/docs/clusters_combined.geojson\n")
}
