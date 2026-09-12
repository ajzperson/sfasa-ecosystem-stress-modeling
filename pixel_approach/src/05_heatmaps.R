library(tidyverse)
library(sf)
library(patchwork)

# Spatial heatmaps of per-pixel CSD metrics, NDVI vs NDWI side by side in
# each PDF: tau_r1 (lag-1 autocorrelation trend), tau_var (variance trend),
# tau_skew (skewness trend). All three share a fixed [-1,1] color scale,
# since Kendall's tau means the same thing regardless of which series it's
# computed from.
#
# A 4th heatmap (Actual_Drop, see deprecated/src/04_actual_drop.R) is kept
# commented out below -- real-world validation against MTBS fire + USFS ADS
# beetle mortality (LCMAP was also tried and dropped -- see
# 13_cluster_agreement_validation.R) is now the trusted ground truth
# instead, and Actual_Drop disagrees with it on cluster severity.
#
# Output: pixel_approach/charts/heatmap_tau_r1.pdf
#         pixel_approach/charts/heatmap_tau_var.pdf
#         pixel_approach/charts/heatmap_tau_skew.pdf

dir.create("pixel_approach/charts", showWarnings = FALSE, recursive = TRUE)

cat("Loading data...\n")
d_ndvi <- read_csv("pixel_approach/data/processed/ndvi_pixel_tau.csv",
                    show_col_types = FALSE) |>
  filter(!is.na(tau_r1))
d_ndwi <- read_csv("pixel_approach/data/processed/ndwi_pixel_tau.csv",
                    show_col_types = FALSE) |>
  filter(!is.na(tau_r1))

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

# Plain polygon on coord_fixed(), not coord_sf() -- its graticule labels
# ignore theme's bold face, so degree labels are recreated manually.
boundary_xy <- as.data.frame(st_coordinates(sn_boundary))
lon_labels <- function(x) sprintf("%g°%s", abs(x), ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
lat_labels <- function(y) sprintf("%g°%s", abs(y), ifelse(y < 0, "S", ifelse(y > 0, "N", "")))

plot_heatmap <- function(data, col, panel_title, low_col, high_col,
                          limits = NULL, legend_name = col,
                          legend_breaks = waiver(), legend_labels = waiver(),
                          legend_reverse = FALSE, diverging = TRUE) {
  ggplot(data, aes(x = lon, y = lat, color = .data[[col]])) +
    geom_point(size = 0.05, alpha = 0.6) +
    geom_polygon(data = boundary_xy, aes(x = X, y = Y), fill = NA, color = "grey30",
                 linewidth = 0.4, inherit.aes = FALSE) +
    (if (diverging) {
      scale_color_gradient2(
        low      = low_col,
        mid      = "white",
        high     = high_col,
        midpoint = 0,
        limits   = limits,
        name     = legend_name,
        breaks   = legend_breaks,
        labels   = legend_labels,
        guide    = guide_colorbar(reverse = legend_reverse)
      )
    } else {
      # One-sided (e.g. |tau_skew|, 0..1): plain low->high gradient, no
      # midpoint -- there's no "opposite direction" to diverge toward.
      scale_color_gradient(
        low    = low_col,
        high   = high_col,
        limits = limits,
        name   = legend_name,
        breaks = legend_breaks,
        labels = legend_labels,
        guide  = guide_colorbar(reverse = legend_reverse)
      )
    }) +
    coord_fixed(
      ratio = 1 / cos(median(data$lat) * pi / 180),
      xlim  = range(data$lon),
      ylim  = range(data$lat)
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3), labels = lon_labels) +
    scale_y_continuous(labels = lat_labels) +
    labs(title = panel_title, x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = 11, color = "#1B2A1F"),
      legend.text     = element_text(face = "bold", size = 10, color = "#1B2A1F"),
      axis.title      = element_text(face = "bold", color = "#1B2A1F"),
      axis.text       = element_text(face = "bold", size = 10, color = "#1B2A1F")
    )
}

save_combined <- function(p_ndvi, p_ndwi, title, out_file, collect_guides = TRUE) {
  p <- (p_ndvi + p_ndwi) +
    plot_layout(guides = if (collect_guides) "collect" else "keep") +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"))
    )
  ggsave(out_file, plot = p, width = 14, height = 7, dpi = 150)
  cat(sprintf("✓ Saved: %s\n", out_file))
}

save_combined(
  plot_heatmap(d_ndvi, "tau_r1", "NDVI", "#2166AC", "#D73027", limits = c(-1, 1)),
  plot_heatmap(d_ndwi, "tau_r1", "NDWI", "#2166AC", "#D73027", limits = c(-1, 1)),
  title    = "Kendall's tau - Lag-1 Autocorrelation (r1)",
  out_file = "pixel_approach/charts/heatmap_tau_r1.pdf"
)

save_combined(
  plot_heatmap(d_ndvi, "tau_var", "NDVI", "#2166AC", "#D73027", limits = c(-1, 1)),
  plot_heatmap(d_ndwi, "tau_var", "NDWI", "#2166AC", "#D73027", limits = c(-1, 1)),
  title    = "Kendall's tau - Variance",
  out_file = "pixel_approach/charts/heatmap_tau_var.pdf"
)

# Simple two-color diverging scale, same structure as tau_r1/tau_var --
# this is a forest-only study area (no desert/dryland pixels), so the
# dense-forest mechanism (lower/more-negative tau_skew = worse) applies
# uniformly; no need to hedge the positive direction separately.
save_combined(
  plot_heatmap(d_ndvi, "tau_skew", "NDVI", "#D73027", "#2166AC", limits = c(-1, 1)),
  plot_heatmap(d_ndwi, "tau_skew", "NDWI", "#D73027", "#2166AC", limits = c(-1, 1)),
  title    = "Kendall's tau - Skewness\n(negative = more stressed; red = worse)",
  out_file = "pixel_approach/charts/heatmap_tau_skew.pdf"
)

# Deprioritized -- see header note.
# save_combined(
#   plot_heatmap(d_ndvi, "Actual_Drop_15_16", "NDVI", "#2166AC", "#D73027",
#                legend_name = "Actual Drop"),
#   plot_heatmap(d_ndwi, "Actual_NDWI_Drop_15_16", "NDWI", "#2166AC", "#D73027",
#                legend_name = "Actual Drop"),
#   title          = "Actual Drop 2015-2016\n(baseline 2000-2011 mean - min 2015-2016)",
#   out_file       = "pixel_approach/charts/heatmap_actual_drop.pdf",
#   collect_guides = FALSE
# )
