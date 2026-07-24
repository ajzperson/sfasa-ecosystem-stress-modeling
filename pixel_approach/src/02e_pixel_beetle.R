library(tidyverse)
library(sf)

# Extracts USFS ADS bark-beetle mortality (TPA -- trees per acre killed) for
# every pixel in the shared grid (pixel_metadata.csv) -- independent of
# clustering, same role as 02c_pixel_fire.R / 02d_pixel_tree_mortality.R.
#
# ALL_YEARS stays 2015:2016 -- the raw ADS data itself only has usable
# LEGACY_TPA for 2015-2016 (near-zero from 2017, none from 2018 on --
# verified against the raw geodatabase).
#
# Output: pixel_approach/data/processed/pixel_beetle.csv
#   (pixel_id, tpa_2015, tpa_2016)

sf_use_s2(FALSE)   # IDS polygons have self-intersecting edges

ALL_YEARS <- 2015:2016
CACHE_DIR <- "pixel_approach/data/processed"

bark_beetle_dca <- c(
  11000, 11002, 11004, 11006, 11007, 11019, 11027,
  11029, 11030, 11035, 11036, 11039, 11041, 11048,
  11049, 11050, 11900
)

gdb_zip <- "data_common/raw/forest_damage/CONUS_Region5_AllYears.gdb.zip"
gdb_dir <- "data_common/raw/forest_damage/CONUS_Region5_AllYears.gdb"

if (!dir.exists(gdb_dir)) {
  cat("Unzipping beetle geodatabase...\n")
  unzip(gdb_zip, exdir = "data_common/raw/forest_damage")
  unlink("data_common/raw/forest_damage/__MACOSX", recursive = TRUE)
}

cat("Reading beetle polygons from GDB...\n")
beetle <- st_read(gdb_dir,
                  layer = "DAMAGE_AREAS_FLAT_AllYears_CONUS_Rgn5",
                  quiet = TRUE) |>
  filter(REGION_ID    == 5,
         SURVEY_YEAR  %in% ALL_YEARS,
         DCA_CODE     %in% bark_beetle_dca,
         DAMAGE_TYPE  == "Mortality") |>
  distinct(DAMAGE_AREA_ID, SURVEY_YEAR, .keep_all = TRUE) |>
  mutate(tpa = if_else(LEGACY_TPA <= 0, NA_real_, as.numeric(LEGACY_TPA))) |>
  st_transform(4326)

cat(sprintf("✓ %d beetle polygons (2015-2016, bark beetle)\n", nrow(beetle)))
cat(sprintf("  TPA range: %.1f to %.1f (NAs: %d / %d)\n",
            min(beetle$tpa, na.rm = TRUE), max(beetle$tpa, na.rm = TRUE),
            sum(is.na(beetle$tpa)), nrow(beetle)))

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE) |>
  st_transform(4326)

beetle_sn <- st_filter(beetle, sn_boundary)
cat(sprintf("✓ %d polygons within Sierra Nevada\n", nrow(beetle_sn)))

pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv",
                        show_col_types = FALSE)
cat(sprintf("✓ %d pixels in shared grid\n", nrow(pixel_meta)))

pts <- st_as_sf(pixel_meta, coords = c("lon", "lat"), crs = 4326)

cat("Joining pixel centroids to beetle polygons (this may take a minute)...\n")
joined_raw <- st_join(pts, beetle_sn |> select(tpa, SURVEY_YEAR), join = st_intersects) |>
  st_drop_geometry()

yearly <- joined_raw |>
  filter(!is.na(SURVEY_YEAR)) |>
  group_by(pixel_id, SURVEY_YEAR) |>
  summarise(tpa = max(tpa, na.rm = TRUE), .groups = "drop") |>
  mutate(tpa = if_else(is.infinite(tpa) | is.nan(tpa), 0, tpa)) |>
  pivot_wider(names_from = SURVEY_YEAR, values_from = tpa,
              names_prefix = "tpa_", values_fill = 0)

tpa_cols <- sprintf("tpa_%d", ALL_YEARS)
out <- pixel_meta |>
  select(pixel_id) |>
  left_join(yearly, by = "pixel_id")
for (col in tpa_cols) if (!col %in% names(out)) out[[col]] <- 0
out[tpa_cols][is.na(out[tpa_cols])] <- 0

for (col in tpa_cols) {
  cat(sprintf("Pixels with beetle mortality in %s: %d / %d (%.1f%%)\n",
              col, sum(out[[col]] > 0), nrow(out), 100 * mean(out[[col]] > 0)))
}

out_path <- file.path(CACHE_DIR, "pixel_beetle.csv")
write_csv(out, out_path)
cat(sprintf("✓ Saved: %s\n", out_path))
