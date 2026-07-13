library(tidyverse)
library(sf)
library(terra)
library(patchwork)

# Cross-validates K-means clusters against USFS ADS beetle-kill polygons.
# LEGACY_TPA (trees per acre killed) is ground-truth mortality.
# Strategy: spatial join of 500m pixel centroids to beetle polygons (2015-2016,
# bark beetle DCA codes), keeping SURVEY_YEAR through the join so per-year
# TPA/mortality can be broken out separately (tpa_2015, tpa_2016, ...) --
# same pattern as pixel_fire.csv (burned_<year>) and pixel_tree_mortality.csv
# (lost_tree_cover_<year>) -- then OR'd into one combined VALIDATION_YEARS
# window for the zonal summary/chart. The raw ADS data itself only has
# usable LEGACY_TPA for 2015-2016 (near-zero from 2017, none from 2018 on --
# verified against the raw geodatabase), so ALL_YEARS stays 2015:2016
# regardless of what VALIDATION_YEARS is narrowed to.
#
# Output:
#   pixel_approach/charts/beetle_tpa_by_cluster.pdf (NDVI + NDWI side by side)
#   pixel_approach/data/processed/pixel_beetle_tpa_ndvi.csv
#   pixel_approach/data/processed/pixel_beetle_tpa_ndwi.csv

sf_use_s2(FALSE)   # IDS polygons have self-intersecting edges

ALL_YEARS        <- 2015:2016
VALIDATION_YEARS <- 2015:2016
year_label <- if (length(VALIDATION_YEARS) == 1) as.character(VALIDATION_YEARS) else
  sprintf("%d-%d", min(VALIDATION_YEARS), max(VALIDATION_YEARS))

OUT_CHART <- "pixel_approach/charts/beetle_tpa_by_cluster.pdf"

index_config <- list(
  ndvi = list(pixel_path = "pixel_approach/data/processed/ndvi_pixel_clusters.csv",
              label      = "NDVI",
              out_csv    = "pixel_approach/data/processed/pixel_beetle_tpa_ndvi.csv"),
  ndwi = list(pixel_path = "pixel_approach/data/processed/ndwi_pixel_clusters.csv",
              label      = "NDWI",
              out_csv    = "pixel_approach/data/processed/pixel_beetle_tpa_ndwi.csv")
)

# --------------------------------------------------------------------------
# 1. Load and filter beetle polygons (shared across both indices — the
#    beetle data itself doesn't depend on which vegetation index is used)
# --------------------------------------------------------------------------
cat("Reading beetle polygons from GDB...\n")

bark_beetle_dca <- c(
  11000, 11002, 11004, 11006, 11019, 11027,
  11029, 11030, 11036, 11048, 11050, 11900
)

gdb_zip <- "data_common/raw/forest_damage/CONUS_Region5_AllYears.gdb.zip"
gdb_dir <- "data_common/raw/forest_damage/CONUS_Region5_AllYears.gdb"

if (!dir.exists(gdb_dir)) {
  cat("Unzipping beetle geodatabase...\n")
  unzip(gdb_zip, exdir = "data_common/raw/forest_damage")
  unlink("data_common/raw/forest_damage/__MACOSX", recursive = TRUE)
}

beetle <- st_read(gdb_dir,
                  layer = "DAMAGE_AREAS_FLAT_AllYears_CONUS_Rgn5",
                  quiet = TRUE) |>
  filter(REGION_ID    == 5,
         SURVEY_YEAR  %in% ALL_YEARS,
         DCA_CODE     %in% bark_beetle_dca) |>
  distinct(DAMAGE_AREA_ID, SURVEY_YEAR, .keep_all = TRUE) |>
  mutate(tpa = if_else(LEGACY_TPA <= 0, NA_real_, as.numeric(LEGACY_TPA))) |>
  st_transform(4326)

cat(sprintf("✓ %d beetle polygons (2015-2016, bark beetle)\n", nrow(beetle)))
cat(sprintf("  TPA range: %.1f to %.1f (NAs: %d / %d)\n",
            min(beetle$tpa, na.rm=TRUE), max(beetle$tpa, na.rm=TRUE),
            sum(is.na(beetle$tpa)), nrow(beetle)))

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE) |>
  st_transform(4326)

beetle_sn <- st_filter(beetle, sn_boundary)
cat(sprintf("✓ %d polygons within Sierra Nevada\n", nrow(beetle_sn)))

# --------------------------------------------------------------------------
# Per-index: join pixels to beetle polygons, summarize, save, plot
# --------------------------------------------------------------------------
run_beetle_validation <- function(cfg, beetle_sn) {
  cat(sprintf("\n=== %s ===\n", cfg$label))
  cat("Loading pixel clusters...\n")
  pixels <- read_csv(cfg$pixel_path, show_col_types = FALSE) |>
    select(pixel_id, lon, lat, cluster)

  pts <- st_as_sf(pixels, coords = c("lon", "lat"), crs = 4326)
  cat(sprintf("✓ %d pixel centroids\n", nrow(pts)))

  cat("Joining pixel centroids to beetle polygons (this may take a minute)...\n")
  # Each pixel may overlap multiple polygons -- take max TPA per pixel,
  # per year (SURVEY_YEAR kept through the join, unlike before, so this
  # doesn't need re-running to change VALIDATION_YEARS later).
  joined_raw <- st_join(pts, beetle_sn |> select(tpa, SURVEY_YEAR), join = st_intersects) |>
    st_drop_geometry()

  yearly <- joined_raw |>
    filter(!is.na(SURVEY_YEAR)) |>
    group_by(pixel_id, cluster, SURVEY_YEAR) |>
    summarise(tpa = max(tpa, na.rm = TRUE), .groups = "drop") |>
    mutate(tpa = if_else(is.infinite(tpa) | is.nan(tpa), 0, tpa)) |>
    pivot_wider(names_from = SURVEY_YEAR, values_from = tpa,
                names_prefix = "tpa_", values_fill = 0)

  tpa_cols <- sprintf("tpa_%d", ALL_YEARS)
  joined <- pixels |>
    select(pixel_id, cluster) |>
    left_join(yearly, by = c("pixel_id", "cluster"))
  for (col in tpa_cols) if (!col %in% names(joined)) joined[[col]] <- 0
  joined[tpa_cols][is.na(joined[tpa_cols])] <- 0

  # Combine into one flag for the requested VALIDATION_YEARS window
  win_cols  <- sprintf("tpa_%d", VALIDATION_YEARS)
  joined$tpa <- do.call(pmax, c(unname(joined[win_cols]), na.rm = TRUE))
  joined$any_mortality <- joined$tpa > 0

  cat(sprintf("✓ Join complete\n"))
  cat(sprintf("  Pixels with beetle mortality (%s): %d / %d (%.1f%%)\n",
              year_label, sum(joined$any_mortality), nrow(joined),
              100 * mean(joined$any_mortality)))

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

  cat(sprintf("\n=== %s: Zonal TPA by Cluster ===\n", cfg$label))
  print(zonal)

  write_csv(joined, cfg$out_csv)
  cat(sprintf("✓ Saved: %s\n", cfg$out_csv))

  cluster_cols <- setNames(
    colorRampPalette(c("#2166AC", "#92C5DE", "#F4A582", "#D73027"))(3),
    c("1", "2", "3")
  )
  cluster_labels <- c("1" = "1 — no CSD", "2" = "2 — partial", "3" = "3 — full CSD")

  p_tpa <- ggplot(zonal, aes(x = factor(cluster), y = pct_with_mortality,
                              fill = factor(cluster))) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%\n(%.1f TPA)",
                                  pct_with_mortality, mean_tpa_all)),
              vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_x_discrete(labels = cluster_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = cfg$label,
      x = "Cluster", y = "% pixels with beetle mortality"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  list(zonal = zonal, plot = p_tpa)
}

results <- lapply(index_config, run_beetle_validation, beetle_sn = beetle_sn)

p_combined <- (results$ndvi$plot + results$ndwi$plot) +
  plot_annotation(
    title    = sprintf("USFS ADS Beetle Mortality by CSD Cluster (%s)", year_label),
    subtitle = "% of pixels with bark beetle mortality — ground-truth validation (mean TPA shown for reference)",
    theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
                      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))
  )

ggsave(OUT_CHART, plot = p_combined, width = 12, height = 5, dpi = 150)
cat(sprintf("✓ Saved: %s\n", OUT_CHART))

# Individual per-index panels, saved separately so 11_dashboard.R can
# assemble a one-page-per-index summary without duplicating this plotting code.
for (idx in names(index_config)) {
  panel_path <- sprintf("pixel_approach/data/processed/panel_beetle_%s.rds", idx)
  saveRDS(results[[idx]]$plot, panel_path)
  cat(sprintf("✓ Saved: %s\n", panel_path))
}
