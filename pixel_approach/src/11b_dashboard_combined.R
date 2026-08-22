library(tidyverse)
library(sf)
library(patchwork)

source("pixel_approach/src/clustering_helpers.R")

# One-page summary for the COMBINED (NDVI+NDWI 4-feature) clustering: cluster
# map beside the two real-world validator bars -- same role as
# 11_dashboard.R, but for the combined approach instead of per-index. Only
# one cluster assignment per pixel now, so one output file, not two.
# Assembled from panels saved by 06b_clustering_combined.R/
# 09b_fire_validation_combined.R/08b_beetle_validation_combined.R -- run
# those three first (in any order; this script only reads their saved
# outputs, it doesn't recompute anything itself).
#
# The cluster map panel is rebuilt here directly from
# combined_pixel_clusters.csv rather than loaded from a saved .rds -- a
# ggplot object embeds its full underlying data, and this map plots ~220k
# raw points, which would make a saved panel far larger than needed;
# rebuilding from the CSV (which already persists) is cheap and avoids that
# bloat. Same reasoning as 11_dashboard.R.
#
# Output:
#   pixel_approach/charts/dashboard_combined.pdf

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

clusters_path <- "pixel_approach/data/processed/combined_pixel_clusters.csv"
pixels <- read_csv(clusters_path, show_col_types = FALSE)
K <- max(pixels$cluster, na.rm = TRUE)

style        <- build_cluster_style(K)
cluster_cols <- style$cols
clust_labels <- style$labels

build_map_panel <- function(d) {
  d <- d |> mutate(cluster = factor(cluster, levels = seq_len(K)))
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
      # for the combined pixel population too, matching NDVI/NDWI) -- see
      # 11_dashboard.R for why a side legend was tried and abandoned.
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

cat("Building combined dashboard...\n")

map_panel <- build_map_panel(pixels)

fire_panel   <- readRDS("pixel_approach/data/processed/panel_fire_combined.rds") +
  labs(title = "Fire (MTBS)")
beetle_panel <- readRDS("pixel_approach/data/processed/panel_beetle_combined.rds") +
  labs(title = "Beetle (USFS ADS)")

# Map is tall/narrow (Sierra Nevada's geographic extent), bars are wide/short
# -- stacking them (map above bars) forces both into the same column width,
# leaving the map row mostly blank on the sides. Placing the map beside a
# stacked bar column instead pairs complementary shapes, filling a landscape
# canvas with far less dead space. Same reasoning as 11_dashboard.R.
p <- ((map_panel) | (fire_panel / beetle_panel)) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "NDVI+NDWI-combined: Cluster Map + Real-World Validation Summary",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

out_file <- "pixel_approach/charts/dashboard_combined.pdf"
ggsave(out_file, plot = p, width = 13, height = 7.5, dpi = 150)
cat(sprintf("✓ Saved: %s\n", out_file))
