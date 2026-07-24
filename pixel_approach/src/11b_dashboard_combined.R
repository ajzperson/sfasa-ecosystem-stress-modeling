library(tidyverse)
library(sf)
library(patchwork)

# One-page summary for the COMBINED (NDVI+NDWI 4-feature) clustering: cluster
# map on top, the two real-world validator bars (fire, beetle) below it --
# same role as 11_dashboard.R, but for the combined approach instead of
# per-index. Only one cluster assignment per pixel now, so one output file,
# not two. Assembled from panels saved by 06b_clustering_combined.R/
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

cluster_cols <- setNames(
  colorRampPalette(c("#2166AC", "#92C5DE", "#F4A582", "#D73027"))(K),
  as.character(seq_len(K))
)
clust_labels <- if (K == 1) {
  "1 - all"
} else {
  paste0(seq_len(K), c(" - no CSD", rep(" - partial", max(K - 2, 0)), " - full CSD"))
}

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
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_bw() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      plot.caption    = element_text(hjust = 0.5, color = "grey40", size = 8),
      legend.position = "right",
      axis.text       = element_text(size = 7)
    )
}

cat("Building combined dashboard...\n")

map_panel <- build_map_panel(pixels)

fire_panel   <- readRDS("pixel_approach/data/processed/panel_fire_combined.rds") +
  labs(title = "Fire (MTBS)")
beetle_panel <- readRDS("pixel_approach/data/processed/panel_beetle_combined.rds") +
  labs(title = "Beetle (USFS ADS)")

p <- (map_panel) / (fire_panel | beetle_panel) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "Combined (NDVI+NDWI): Cluster Map + Real-World Validation Summary",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

out_file <- "pixel_approach/charts/dashboard_combined.pdf"
ggsave(out_file, plot = p, width = 12, height = 12, dpi = 150)
cat(sprintf("✓ Saved: %s\n", out_file))
