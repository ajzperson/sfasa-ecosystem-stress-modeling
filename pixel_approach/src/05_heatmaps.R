library(tidyverse)
library(sf)

# Spatial heatmaps of per-pixel CSD and collapse metrics:
#   (1) tau_r1           — lag-1 autocorrelation trend
#   (2) tau_var           — variance trend
#   (3) tau_skew          — skewness trend
#   (4) Actual_Drop_15_16 — observed NDVI drop (baseline 2000-2011 mean minus
#                           min NDVI 2015-2016); positive = collapsed more

dir.create("pixel_approach/charts", showWarnings = FALSE, recursive = TRUE)

cat("Loading data...\n")
d <- read_csv("pixel_approach/data/processed/pixel_tau.csv",
              show_col_types = FALSE) |>
  filter(!is.na(tau_r1))

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

plot_heatmap <- function(data, col, title, low_col, high_col, out_file) {
  p <- ggplot(data, aes(x = lon, y = lat, color = .data[[col]])) +
    geom_point(size = 0.05, alpha = 0.6) +
    geom_sf(data = sn_boundary, fill = NA, color = "grey30",
            linewidth = 0.4, inherit.aes = FALSE) +
    scale_color_gradient2(
      low      = low_col,
      mid      = "white",
      high     = high_col,
      midpoint = 0,
      name     = col
    ) +
    coord_sf(
      xlim = range(data$lon),
      ylim = range(data$lat)
    ) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      axis.text       = element_text(size = 7)
    )
  ggsave(out_file, plot = p, width = 8, height = 7, dpi = 150)
  cat(sprintf("✓ Saved: %s\n", out_file))
}

plot_heatmap(d, "tau_r1",
             "Kendall's tau - Lag-1 Autocorrelation (r1) Trend",
             low_col  = "#2166AC",
             high_col = "#D73027",
             out_file = "pixel_approach/charts/heatmap_tau_r1.pdf")

plot_heatmap(d, "tau_var",
             "Kendall's tau - Variance Trend",
             low_col  = "#2166AC",
             high_col = "#D73027",
             out_file = "pixel_approach/charts/heatmap_tau_var.pdf")

plot_heatmap(d, "tau_skew",
             "Kendall's tau - Skewness Trend\n(negative = more stressed; red = worse)",
             low_col  = "#D73027",
             high_col = "#2166AC",
             out_file = "pixel_approach/charts/heatmap_tau_skew.pdf")

plot_heatmap(d, "Actual_Drop_15_16",
             "Actual NDVI Drop 2015-2016\n(baseline 2000-2011 mean - min NDVI 2015-2016)",
             low_col  = "#2166AC",
             high_col = "#D73027",
             out_file = "pixel_approach/charts/heatmap_actual_drop.pdf")
