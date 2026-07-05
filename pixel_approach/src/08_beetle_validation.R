library(tidyverse)
library(sf)
library(terra)

# Cross-validates K-means clusters against USFS ADS beetle-kill polygons.
# LEGACY_TPA (trees per acre killed) is ground-truth mortality.
# Strategy: spatial join of 500m pixel centroids to beetle polygons (2015-2016,
# bark beetle DCA codes), then zonal summary by cluster.
#
# Output:
#   pixel_approach/charts/beetle_tpa_by_cluster.pdf
#   pixel_approach/data/processed/pixel_beetle_tpa.csv

sf_use_s2(FALSE)   # IDS polygons have self-intersecting edges

# --------------------------------------------------------------------------
# 1. Load pixel clusters
# --------------------------------------------------------------------------
cat("Loading pixel clusters...\n")
pixels <- read_csv("pixel_approach/data/processed/pixel_clusters.csv",
                   show_col_types = FALSE) |>
  select(pixel_id, lon, lat, cluster)

pts <- st_as_sf(pixels, coords = c("lon", "lat"), crs = 4326)
cat(sprintf("✓ %d pixel centroids\n", nrow(pts)))

# --------------------------------------------------------------------------
# 2. Load and filter beetle polygons
# --------------------------------------------------------------------------
cat("Reading beetle polygons from GDB...\n")

bark_beetle_dca <- c(
  11000, 11002, 11004, 11006, 11019, 11027,
  11029, 11030, 11036, 11048, 11050, 11900
)

gdb_zip <- "data_common/raw/beetle/CONUS_Region5_AllYears.gdb.zip"
gdb_dir <- "data_common/raw/beetle/CONUS_Region5_AllYears.gdb"

if (!dir.exists(gdb_dir)) {
  cat("Unzipping beetle geodatabase...\n")
  unzip(gdb_zip, exdir = "data_common/raw/beetle")
  unlink("data_common/raw/beetle/__MACOSX", recursive = TRUE)
}

beetle <- st_read(gdb_dir,
                  layer = "DAMAGE_AREAS_FLAT_AllYears_CONUS_Rgn5",
                  quiet = TRUE) |>
  filter(REGION_ID    == 5,
         SURVEY_YEAR  %in% 2015:2016,
         DCA_CODE     %in% bark_beetle_dca) |>
  distinct(DAMAGE_AREA_ID, SURVEY_YEAR, .keep_all = TRUE) |>
  mutate(tpa = if_else(LEGACY_TPA <= 0, NA_real_, as.numeric(LEGACY_TPA))) |>
  st_transform(st_crs(pts))

cat(sprintf("✓ %d beetle polygons (2015-2016, bark beetle)\n", nrow(beetle)))
cat(sprintf("  TPA range: %.1f to %.1f (NAs: %d / %d)\n",
            min(beetle$tpa, na.rm=TRUE), max(beetle$tpa, na.rm=TRUE),
            sum(is.na(beetle$tpa)), nrow(beetle)))

# --------------------------------------------------------------------------
# 3. Clip beetle polygons to Sierra Nevada extent
# --------------------------------------------------------------------------
sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE) |>
  st_transform(st_crs(pts))

beetle_sn <- st_filter(beetle, sn_boundary)
cat(sprintf("✓ %d polygons within Sierra Nevada\n", nrow(beetle_sn)))

# --------------------------------------------------------------------------
# 4. Spatial join: pixels → beetle polygons
# --------------------------------------------------------------------------
cat("Joining pixel centroids to beetle polygons (this may take a minute)...\n")

# Each pixel may overlap multiple polygons — take max TPA per pixel
joined <- st_join(pts, beetle_sn |> select(tpa), join = st_intersects) |>
  st_drop_geometry() |>
  group_by(pixel_id, cluster) |>
  summarise(tpa = max(tpa, na.rm = TRUE), .groups = "drop") |>
  mutate(tpa = if_else(is.infinite(tpa) | is.nan(tpa), 0, tpa),
         any_mortality = tpa > 0)

cat(sprintf("✓ Join complete\n"))
cat(sprintf("  Pixels with beetle mortality: %d / %d (%.1f%%)\n",
            sum(joined$any_mortality), nrow(joined),
            100 * mean(joined$any_mortality)))

# --------------------------------------------------------------------------
# 5. Zonal statistics by cluster
# --------------------------------------------------------------------------
zonal <- joined |>
  group_by(cluster) |>
  summarise(
    n_pixels          = n(),
    n_with_mortality  = sum(any_mortality),
    pct_with_mortality = round(100 * mean(any_mortality), 1),
    mean_tpa_affected  = round(mean(tpa[any_mortality], na.rm = TRUE), 1),
    mean_tpa_all       = round(mean(tpa, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  arrange(as.integer(cluster))

cat("\n=== Zonal TPA by Cluster ===\n")
print(zonal)

# --------------------------------------------------------------------------
# 6. Save pixel-level beetle data
# --------------------------------------------------------------------------
write_csv(joined, "pixel_approach/data/processed/pixel_beetle_tpa.csv")
cat("✓ Saved: pixel_approach/data/processed/pixel_beetle_tpa.csv\n")

# --------------------------------------------------------------------------
# 7. Plot: mean TPA per cluster (pixels with any mortality only)
# --------------------------------------------------------------------------
cluster_cols <- setNames(
  colorRampPalette(c("#2166AC", "#92C5DE", "#F4A582", "#D73027"))(3),
  c("1", "2", "3")
)
cluster_labels <- c("1" = "1 — no CSD", "2" = "2 — partial", "3" = "3 — full CSD")

p_tpa <- ggplot(zonal, aes(x = factor(cluster), y = mean_tpa_all,
                            fill = factor(cluster))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f TPA\n(%s%% affected)",
                                mean_tpa_all, pct_with_mortality)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = cluster_cols, guide = "none") +
  scale_x_discrete(labels = cluster_labels) +
  labs(
    title    = "USFS ADS Beetle Mortality by CSD Cluster (2015–2016)",
    subtitle = "Mean LEGACY_TPA (trees per acre killed) per cluster — ground-truth validation",
    x = "Cluster", y = "Mean TPA (all pixels)"
  ) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))

ggsave("pixel_approach/charts/beetle_tpa_by_cluster.pdf",
       plot = p_tpa, width = 7, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/beetle_tpa_by_cluster.pdf\n")
