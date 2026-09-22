library(tidyverse)
library(sf)
library(patchwork)

source("pixel_approach/src/clustering_helpers.R")

# One-page-per-index summary: cluster map beside the two real-world
# validator bars (fire, beetle), stacked. Assembled from panels saved by
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
#   pixel_approach/charts/dashboard_ndvi_ratio.pdf (same layout, fire/beetle
#     bars as ratio vs. Cluster 1 instead of raw %) -- generated in parallel
#     from panel_fire_<idx>_ratio.rds/panel_beetle_<idx>_ratio.rds, raw
#     dashboards untouched
#   pixel_approach/charts/dashboard_ndwi_ratio.pdf

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

style        <- build_cluster_style(3)
cluster_cols <- style$cols
clust_labels <- style$labels

# Plain polygon on coord_fixed(), not coord_sf() -- its graticule labels
# ignore theme's bold face, so degree labels are recreated manually.
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
      # Legend inset into the panel's upper-right corner, which is genuinely
      # empty of data (lat > 38.0 & lon > -119.15 has zero pixels, for both
      # NDVI and NDWI) -- a side legend was tried and abandoned:
      # patchwork gives a coord_sf-locked panel a column sized to its own
      # aspect ratio and centers the (narrower) panel within it, leaving a
      # large dead gap before the legend that neither plot_layout() width
      # ratios, legend.box.spacing, nor patchwork::free() could close
      # without either leaving the gap or clipping the legend text.
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

  # Extra bottom/top margin here (not baked into the saved panels, which are
  # also reused standalone in fire_by_cluster.pdf/beetle_tpa_by_cluster.pdf
  # where this gap isn't needed) -- gives the bold "Relative risk" caption at
  # the foot of fire_panel room to breathe before beetle_panel's title starts,
  # which patchwork's "/" stacking doesn't add on its own.
  fire_panel   <- readRDS(sprintf("pixel_approach/data/processed/panel_fire_%s%s.rds", idx, suffix)) +
    labs(title = "Fire (MTBS)") +
    theme(plot.margin = margin(t = 5.5, r = 5.5, b = 34, l = 5.5))
  beetle_panel <- readRDS(sprintf("pixel_approach/data/processed/panel_beetle_%s%s.rds", idx, suffix)) +
    labs(title = "Beetle (USFS ADS)") +
    theme(plot.margin = margin(t = 28, r = 5.5, b = 5.5, l = 5.5))

  # Map is tall/narrow (Sierra Nevada's geographic extent), bars are wide/short
  # -- stacking them (map above bars) forces both into the same column width,
  # leaving the map row mostly blank on the sides. Placing the map beside a
  # stacked bar column instead pairs complementary shapes, filling a landscape
  # canvas with far less dead space.
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
