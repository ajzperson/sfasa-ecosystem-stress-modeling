library(tidyverse)
library(sf)
library(patchwork)

source("pixel_approach/src/clustering_helpers.R")

# One page per index: cluster map beside the fire and beetle bar charts.
#
# Output: pixel_approach/charts/dashboard_{ndvi,ndwi}{,_ratio}.pdf

sn_boundary <- st_read("data_common/processed/sierra_nevada_boundary.geojson", quiet = TRUE)

style        <- build_cluster_style(3)
cluster_cols <- style$cols
clust_labels <- style$labels

# coord_fixed (not coord_sf) so the degree labels can be bold
boundary_xy <- as.data.frame(st_coordinates(sn_boundary))
lon_labels <- function(x) sprintf("%g°%s", abs(x), ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
lat_labels <- function(y) sprintf("%g°%s", abs(y), ifelse(y < 0, "S", ifelse(y > 0, "N", "")))

build_map_panel <- function(clusters_path, label) {
  d <- read_csv(clusters_path, show_col_types = FALSE) |>
    mutate(cluster = factor(cluster, levels = 1:3))
  ggplot(d, aes(x = lon, y = lat, color = cluster)) +
    geom_point(size = 0.05, alpha = 0.5) +
    geom_polygon(data = boundary_xy, aes(x = X, y = Y), fill = NA, color = "grey20",
                 linewidth = 0.4, inherit.aes = FALSE) +
    scale_color_manual(values = cluster_cols, name = "Stress level",
                       labels = clust_labels) +
    coord_fixed(ratio = 1 / cos(median(d$lat) * pi / 180),
                xlim = range(d$lon), ylim = range(d$lat)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3), labels = lon_labels) +
    scale_y_continuous(labels = lat_labels) +
    labs(title = "Cluster Map", x = "Longitude", y = "Latitude",
         caption = "White gaps: pixels burned (MTBS moderate+ severity) at any point\n2000-2014, excluded before clustering") +
    guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
    theme_bw() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
      plot.caption    = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"),
      # Legend inside the map's empty upper-right corner
      legend.position             = "inside",
      legend.position.inside      = c(0.97, 0.97),
      legend.justification.inside = c(1, 1),
      legend.background           = element_rect(fill = "white", color = "grey60"),
      legend.title                = element_text(size = 10, face = "bold", color = "#1B2A1F"),
      legend.text                 = element_text(size = 9, face = "bold", color = "#1B2A1F"),
      legend.key.size             = unit(0.35, "cm"),
      axis.title      = element_text(face = "bold", color = "#1B2A1F"),
      axis.text       = element_text(size = 9, face = "bold", color = "#1B2A1F")
    )
}

index_config <- list(
  ndvi = list(clusters_path = "pixel_approach/data/processed/ndvi_pixel_clusters.csv",
              label = "NDVI", out = "pixel_approach/charts/dashboard_ndvi.pdf",
              out_ratio = "pixel_approach/charts/dashboard_ndvi_ratio.pdf"),
  ndwi = list(clusters_path = "pixel_approach/data/processed/ndwi_pixel_clusters.csv",
              label = "NDWI", out = "pixel_approach/charts/dashboard_ndwi.pdf",
              out_ratio = "pixel_approach/charts/dashboard_ndwi_ratio.pdf")
)

build_dashboard <- function(cfg, idx, ratio = FALSE) {
  cat(sprintf("Building dashboard for %s%s...\n", cfg$label, if (ratio) " (ratio)" else ""))

  map_panel <- build_map_panel(cfg$clusters_path, cfg$label)

  suffix <- if (ratio) "_ratio" else ""

  # Extra space between the fire caption and the beetle title
  fire_panel   <- readRDS(sprintf("pixel_approach/data/processed/panel_fire_%s%s.rds", idx, suffix)) +
    labs(title = "Fire (MTBS)") +
    theme(plot.margin = margin(t = 5.5, r = 5.5, b = 34, l = 5.5))
  beetle_panel <- readRDS(sprintf("pixel_approach/data/processed/panel_beetle_%s%s.rds", idx, suffix)) +
    labs(title = "Beetle (USFS ADS)") +
    theme(plot.margin = margin(t = 28, r = 5.5, b = 5.5, l = 5.5))

  # Tall map on the left, the two bar charts stacked on the right
  title_suffix <- if (ratio) " (Ratio vs. Cluster 1)" else ""
  p <- ((map_panel) | (fire_panel / beetle_panel)) +
    plot_layout(widths = c(1, 1)) +
    plot_annotation(
      title = sprintf("%s: Cluster Map + Real-World Evaluation Summary%s", cfg$label, title_suffix),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "#1B2A1F"))
    )

  out <- if (ratio) cfg$out_ratio else cfg$out
  ggsave(out, plot = p, width = 13, height = 7.5, dpi = 150)
  cat(sprintf("✓ Saved: %s\n", out))
}

for (idx in names(index_config)) {
  build_dashboard(index_config[[idx]], idx, ratio = FALSE)
  build_dashboard(index_config[[idx]], idx, ratio = TRUE)
}
