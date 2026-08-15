library(tidyverse)
library(sf)
library(patchwork)

source("pixel_approach/src/clustering_helpers.R")

# One-page-per-index summary: cluster map on top, the two real-world
# validator bars (fire, beetle) below it. LCMAP was dropped as a validator --
# its detection method (Landsat break detection, tuned to abrupt disturbance)
# doesn't match what CSD is built to detect (slow, gradual pre-tipping-point
# decline); see DESIGN.md. Assembled from panels saved by
# 06_clustering.R/09_fire_validation.R/08_beetle_validation.R -- run those
# three first (in any order; this script only reads their saved outputs, it
# doesn't recompute anything itself).
#
# The cluster map panel is rebuilt here directly from ndvi_pixel_clusters.csv/
# ndwi_pixel_clusters.csv rather than loaded from a saved .rds -- a ggplot
# object embeds its full underlying data, and this map plots ~220k raw
# points, which would make a saved panel ~175MB; rebuilding from the CSV
# (which already persists) is cheap and avoids that bloat.
#
# Output:
#   pixel_approach/charts/dashboard_ndvi.pdf
#   pixel_approach/charts/dashboard_ndwi.pdf

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

style        <- build_cluster_style(3)
cluster_cols <- style$cols
clust_labels <- style$labels

build_map_panel <- function(clusters_path, label) {
  d <- read_csv(clusters_path, show_col_types = FALSE) |>
    mutate(cluster = factor(cluster, levels = 1:3))
  ggplot(d, aes(x = lon, y = lat, color = cluster)) +
    geom_point(size = 0.05, alpha = 0.5) +
    geom_sf(data = sn_boundary, fill = NA, color = "grey20",
            linewidth = 0.4, inherit.aes = FALSE) +
    scale_color_manual(values = cluster_cols, name = "Stress level",
                       labels = clust_labels) +
    coord_sf(xlim = range(d$lon), ylim = range(d$lat)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    labs(title = "Cluster Map", x = "Longitude", y = "Latitude",
         caption = "White gaps: pixels burned (MTBS moderate+ severity) at any point\n2000-2014, excluded before clustering") +
    guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
    theme_bw() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      plot.caption    = element_text(hjust = 0.5, color = "grey40", size = 8),
      # Legend inset into the panel's upper-right corner, which is genuinely
      # empty of data (lat > 38.0 & lon > -119.15 has zero pixels, verified
      # for both NDVI and NDWI) -- a side legend was tried and abandoned:
      # patchwork gives a coord_sf-locked panel a column sized to its own
      # aspect ratio and centers the (narrower) panel within it, leaving a
      # large dead gap before the legend that neither plot_layout() width
      # ratios, legend.box.spacing, nor patchwork::free() could close
      # without either leaving the gap or clipping the legend text.
      legend.position             = "inside",
      legend.position.inside      = c(0.97, 0.97),
      legend.justification.inside = c(1, 1),
      legend.background           = element_rect(fill = "white", color = "grey60"),
      legend.title                = element_text(size = 9),
      legend.text                 = element_text(size = 8),
      legend.key.size             = unit(0.35, "cm"),
      axis.text       = element_text(size = 7)
    )
}

index_config <- list(
  ndvi = list(clusters_path = "pixel_approach/data/processed/ndvi_pixel_clusters.csv",
              label = "NDVI", out = "pixel_approach/charts/dashboard_ndvi.pdf"),
  ndwi = list(clusters_path = "pixel_approach/data/processed/ndwi_pixel_clusters.csv",
              label = "NDWI", out = "pixel_approach/charts/dashboard_ndwi.pdf")
)

build_dashboard <- function(cfg, idx) {
  cat(sprintf("Building dashboard for %s...\n", cfg$label))

  map_panel <- build_map_panel(cfg$clusters_path, cfg$label)

  fire_panel   <- readRDS(sprintf("pixel_approach/data/processed/panel_fire_%s.rds", idx)) +
    labs(title = "Fire (MTBS)")
  beetle_panel <- readRDS(sprintf("pixel_approach/data/processed/panel_beetle_%s.rds", idx)) +
    labs(title = "Beetle (USFS ADS)")

  # Map is tall/narrow (Sierra Nevada's geographic extent), bars are wide/short
  # -- stacking them (map above bars) forces both into the same column width,
  # leaving the map row mostly blank on the sides. Placing the map beside a
  # stacked bar column instead pairs complementary shapes, filling a landscape
  # canvas with far less dead space.
  p <- ((map_panel) | (fire_panel / beetle_panel)) +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(
      title = sprintf("%s: Cluster Map + Real-World Validation Summary", cfg$label),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )

  ggsave(cfg$out, plot = p, width = 13, height = 7.5, dpi = 150)
  cat(sprintf("✓ Saved: %s\n", cfg$out))
}

for (idx in names(index_config)) {
  build_dashboard(index_config[[idx]], idx)
}
