library(tidyverse)
library(sf)

# Cross-validates K-means clusters against CAL FIRE historical fire perimeters.
# Two separate checks, same fire dataset, different time windows:
#   1. Baseline window (2000-2014): did any fire burn a pixel during the period
#      its pre-collapse warning-sign statistics were computed from? A "yes" here
#      is a confound flag, not validation -- a fire-caused NDVI/NDWI drop during
#      the baseline would corrupt the input signal, unrelated to CSD theory.
#   2. Validation window (2015-2019, matching the beetle-mortality window for
#      the same detection-lag/single-drought-event reasoning): did a fire burn
#      a pixel after the predicted-risk clusters were assigned? A "yes" here is
#      a second, independent ground-truth outcome, same role as beetle mortality.
#
# Output (per index):
#   pixel_approach/charts/fire_by_ndvi_cluster.pdf
#   pixel_approach/charts/fire_by_ndwi_cluster.pdf
#   pixel_approach/data/processed/pixel_fire_ndvi.csv
#   pixel_approach/data/processed/pixel_fire_ndwi.csv

sf_use_s2(FALSE)

index_config <- list(
  ndvi = list(pixel_path = "pixel_approach/data/processed/ndvi_pixel_clusters.csv",
              label      = "NDVI",
              out_csv    = "pixel_approach/data/processed/pixel_fire_ndvi.csv",
              out_chart  = "pixel_approach/charts/fire_by_ndvi_cluster.pdf"),
  ndwi = list(pixel_path = "pixel_approach/data/processed/ndwi_pixel_clusters.csv",
              label      = "NDWI",
              out_csv    = "pixel_approach/data/processed/pixel_fire_ndwi.csv",
              out_chart  = "pixel_approach/charts/fire_by_ndwi_cluster.pdf")
)

# --------------------------------------------------------------------------
# 1. Load fire perimeters, clip to Sierra Nevada (shared across both indices)
# --------------------------------------------------------------------------
cat("Loading fire perimeters...\n")
fires <- st_read("data_common/raw/fire_perimeters/California_Historic_Fire_Perimeters_-4934337877712688794.gpkg",
                  quiet = TRUE) |>
  filter(!is.na(YEAR_)) |>
  st_transform(4326)

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE) |>
  st_transform(4326)

fires_sn <- st_filter(fires, sn_boundary)
cat(sprintf("✓ %d fire perimeters within Sierra Nevada\n", nrow(fires_sn)))

baseline_fires    <- fires_sn |> filter(YEAR_ >= 2000, YEAR_ <= 2014)
validation_fires  <- fires_sn |> filter(YEAR_ >= 2015, YEAR_ <= 2019)
cat(sprintf("  Baseline window (2000-2014): %d fires\n", nrow(baseline_fires)))
cat(sprintf("  Validation window (2015-2019): %d fires\n", nrow(validation_fires)))

# --------------------------------------------------------------------------
# Per-index: flag pixels burned in each window, summarize, save, plot
# --------------------------------------------------------------------------
run_fire_validation <- function(cfg, baseline_fires, validation_fires) {
  cat(sprintf("\n=== %s ===\n", cfg$label))
  pixels <- read_csv(cfg$pixel_path, show_col_types = FALSE) |>
    select(pixel_id, lon, lat, cluster)
  pts <- st_as_sf(pixels, coords = c("lon", "lat"), crs = 4326)
  cat(sprintf("✓ %d pixel centroids\n", nrow(pts)))

  hit_baseline   <- st_intersects(pts, baseline_fires, sparse = TRUE)
  hit_validation <- st_intersects(pts, validation_fires, sparse = TRUE)

  pixels$burned_baseline_00_14   <- lengths(hit_baseline) > 0
  pixels$burned_validation_15_19 <- lengths(hit_validation) > 0

  cat(sprintf("  Pixels burned in baseline window (2000-2014, potential confound): %d / %d (%.1f%%)\n",
              sum(pixels$burned_baseline_00_14), nrow(pixels),
              100 * mean(pixels$burned_baseline_00_14)))
  cat(sprintf("  Pixels burned in validation window (2015-2019): %d / %d (%.1f%%)\n",
              sum(pixels$burned_validation_15_19), nrow(pixels),
              100 * mean(pixels$burned_validation_15_19)))

  zonal <- pixels |>
    group_by(cluster) |>
    summarise(
      n_pixels                  = n(),
      pct_burned_baseline_00_14 = round(100 * mean(burned_baseline_00_14), 1),
      pct_burned_validation     = round(100 * mean(burned_validation_15_19), 1),
      .groups = "drop"
    ) |>
    arrange(as.integer(cluster))

  cat(sprintf("\n=== %s: Fire Exposure by Cluster ===\n", cfg$label))
  print(zonal)

  write_csv(pixels, cfg$out_csv)
  cat(sprintf("✓ Saved: %s\n", cfg$out_csv))

  cluster_cols <- setNames(
    colorRampPalette(c("#2166AC", "#92C5DE", "#F4A582", "#D73027"))(3),
    c("1", "2", "3")
  )
  cluster_labels <- c("1" = "1 — no CSD", "2" = "2 — partial", "3" = "3 — full CSD")

  p_fire <- ggplot(zonal, aes(x = factor(cluster), y = pct_burned_validation,
                              fill = factor(cluster))) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", pct_burned_validation)),
              vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_x_discrete(labels = cluster_labels) +
    labs(
      title    = sprintf("CAL FIRE Perimeter Overlap by %s CSD Cluster (2015–2019)", cfg$label),
      subtitle = "% of pixels burned after cluster assignment — ground-truth validation",
      x = "Cluster", y = "% pixels burned"
    ) +
    theme_bw() +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))

  ggsave(cfg$out_chart, plot = p_fire, width = 7, height = 5, dpi = 150)
  cat(sprintf("✓ Saved: %s\n", cfg$out_chart))

  invisible(zonal)
}

for (idx in index_config) {
  run_fire_validation(idx, baseline_fires, validation_fires)
}
