library(tidyverse)
library(sf)
source("county_approach/src/cluster_polygons.R")

zip_file <- "data_common/raw/beetle/CONUS_Region5_AllYears.gdb.zip"
out_dir  <- "data_common/raw/beetle/CONUS_Region5_AllYears.gdb"

if (!dir.exists(out_dir)) {
  cat("Unzipping beetle geodatabase...\n")
  unzip(zip_file, exdir = "data_common/raw/beetle")
  unlink("data_common/raw/beetle/__MACOSX", recursive = TRUE)
}

gdb   <- out_dir
layer <- "DAMAGE_AREAS_FLAT_AllYears_CONUS_Rgn5"

beetle_dca <- c(
  11000,  # bark beetles (generic — large historical catch-all code)
  11002,  # western pine beetle
  11004,  # Jeffrey pine beetle
  11006,  # mountain pine beetle
  11019,  # pinyon ips
  11027,  # California fivespined ips
  11029,  # pine engraver
  11030,  # ips engraver beetles
  11036,  # western cedar bark beetle
  11048,  # true fir bark beetles
  11050,  # fir engraver
  11900   # unknown bark beetle
)

sf_use_s2(FALSE)  # some IDS polygons have self-intersecting edges

cat("Reading beetle damage polygons from geodatabase...\n")
beetle_raw <- st_read(gdb, layer = layer, quiet = TRUE) |>
  filter(REGION_ID == 5,
         SURVEY_YEAR >= 2000, SURVEY_YEAR <= 2016,
         DCA_CODE %in% beetle_dca) |>
  # a polygon can appear once per DCA code with identical ACRES each time;
  # IDs are unique per year so deduplicating by ID+year is safe and explicit
  distinct(DAMAGE_AREA_ID, SURVEY_YEAR, .keep_all = TRUE)

cat(sprintf("✓ %d polygons after filtering (2000-2016, bark beetle)\n",
            nrow(beetle_raw)))

# Reproject to match cluster polygons
beetle_raw <- st_transform(beetle_raw, st_crs(sequoia_poly))

# Assign to cluster using centroid (avoids double-counting boundary polygons)
centroids  <- st_centroid(beetle_raw)
in_sequoia <- st_intersects(centroids, sequoia_poly, sparse = FALSE)[, 1]
in_plumas  <- st_intersects(centroids, plumas_poly,  sparse = FALSE)[, 1]

beetle_raw <- beetle_raw |>
  mutate(cluster = case_when(
    in_sequoia ~ "Sequoia",
    in_plumas  ~ "Plumas",
    TRUE       ~ NA_character_
  )) |>
  filter(!is.na(cluster))

cat(sprintf("  Sequoia: %d polygons | Plumas: %d polygons\n",
            sum(beetle_raw$cluster == "Sequoia"),
            sum(beetle_raw$cluster == "Plumas")))

# Aggregate to annual acres per cluster
beetle_annual <- beetle_raw |>
  st_drop_geometry() |>
  group_by(cluster, year = SURVEY_YEAR) |>
  summarise(beetle_acres = round(sum(ACRES, na.rm = TRUE)),
            .groups = "drop") |>
  arrange(cluster, year)

cat("\nBeetle acres by year:\n")
print(beetle_annual |>
        pivot_wider(names_from = cluster, values_from = beetle_acres,
                    values_fill = 0))

out_file <- "county_approach/data/processed/Sierra_Nevada_BeetleAcres_2000-2016.csv"
write_csv(beetle_annual, out_file)
cat(sprintf("✓ Saved: %s\n", out_file))
