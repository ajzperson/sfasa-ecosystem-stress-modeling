library(tidyverse)
library(patchwork)

source("pixel_approach/src/clustering_helpers.R")

# Cross-validates K-means clusters against USFS ADS beetle-kill mortality
# (extracted in 02e_pixel_beetle.R, independent of clustering) -- same role
# as 09_fire_validation.R.
#
# Validation window (VALIDATION_YEARS below, default 2015-2016) -- the raw
# ADS data itself only has usable LEGACY_TPA for 2015-2016 (see 02e for why
# ALL_YEARS there stays fixed at 2015:2016 regardless of this window).
#
# Output:
#   pixel_approach/charts/beetle_tpa_by_cluster.pdf (NDVI + NDWI side by side)
#   pixel_approach/charts/beetle_tpa_by_cluster_ratio.pdf (same, bars as ratio
#     vs. Cluster 1 instead of raw %) -- generated in parallel, raw
#     chart/panels untouched
#   pixel_approach/data/processed/pixel_beetle_tpa_ndvi.csv
#   pixel_approach/data/processed/pixel_beetle_tpa_ndwi.csv
#   pixel_approach/data/processed/panel_beetle_ndvi.rds (per-index panel,
#     reused by 11_dashboard.R)
#   pixel_approach/data/processed/panel_beetle_ndvi_ratio.rds
#   pixel_approach/data/processed/panel_beetle_ndwi.rds
#   pixel_approach/data/processed/panel_beetle_ndwi_ratio.rds

VALIDATION_YEARS <- 2015:2016
year_label <- if (length(VALIDATION_YEARS) == 1) as.character(VALIDATION_YEARS) else
  sprintf("%d-%d", min(VALIDATION_YEARS), max(VALIDATION_YEARS))

OUT_CHART       <- "pixel_approach/charts/beetle_tpa_by_cluster.pdf"
OUT_CHART_RATIO <- "pixel_approach/charts/beetle_tpa_by_cluster_ratio.pdf"

index_config <- list(
  ndvi = list(pixel_path = "pixel_approach/data/processed/ndvi_pixel_clusters.csv",
              label      = "NDVI",
              out_csv    = "pixel_approach/data/processed/pixel_beetle_tpa_ndvi.csv"),
  ndwi = list(pixel_path = "pixel_approach/data/processed/ndwi_pixel_clusters.csv",
              label      = "NDWI",
              out_csv    = "pixel_approach/data/processed/pixel_beetle_tpa_ndwi.csv")
)

beetle_flags <- read_csv("pixel_approach/data/processed/pixel_beetle.csv", show_col_types = FALSE)

# OR together (via pmax on continuous TPA) whichever per-year tpa_<year>
# columns fall in VALIDATION_YEARS into one combined validation-window value.
win_cols <- sprintf("tpa_%d", VALIDATION_YEARS)
stopifnot(all(win_cols %in% names(beetle_flags)))
beetle_flags$tpa <- do.call(pmax, c(unname(beetle_flags[win_cols]), na.rm = TRUE))
beetle_flags$any_mortality <- beetle_flags$tpa > 0

run_beetle_validation <- function(cfg, beetle_flags) {
  cat(sprintf("\n=== %s ===\n", cfg$label))
  pixels <- read_csv(cfg$pixel_path, show_col_types = FALSE) |>
    select(pixel_id, cluster) |>
    left_join(beetle_flags, by = "pixel_id")

  cat(sprintf("  Pixels with beetle mortality (%s): %d / %d (%.1f%%)\n",
              year_label, sum(pixels$any_mortality), nrow(pixels),
              100 * mean(pixels$any_mortality)))

  zonal <- pixels |>
    group_by(cluster) |>
    summarise(
      n_pixels                 = n(),
      n_with_mortality         = sum(any_mortality),
      pct_with_mortality_exact = 100 * mean(any_mortality),
      pct_with_mortality       = round(pct_with_mortality_exact, 1),
      mean_tpa_affected        = round(mean(tpa[any_mortality], na.rm = TRUE), 1),
      mean_tpa_all             = round(mean(tpa, na.rm = TRUE), 1),
      .groups = "drop"
    ) |>
    arrange(as.integer(cluster))

  cat(sprintf("\n=== %s: Zonal TPA by Cluster ===\n", cfg$label))
  print(zonal)

  # Relative risk and ratio_to_c1 are computed from the unrounded rate
  # (pct_with_mortality_exact), not the 1-decimal display value -- dividing
  # already-rounded percentages compounds their rounding error into the
  # ratio (e.g. a 2.9%/1.8% pair that's really 2.85/1.85 would read 1.61x
  # instead of the true 1.54x).
  rr_c3_vs_c1 <- zonal$pct_with_mortality_exact[zonal$cluster == 3] /
    zonal$pct_with_mortality_exact[zonal$cluster == 1]
  cat(sprintf("  Relative risk (Cluster 3 vs 1): %.2fx\n", rr_c3_vs_c1))

  # Ratio vs. Cluster 1 (least-stressed), for the ratio-chart variant below --
  # Cluster 1 is 1.00x by construction.
  zonal <- zonal |>
    mutate(ratio_to_c1 = pct_with_mortality_exact / pct_with_mortality_exact[cluster == 1])

  write_csv(pixels, cfg$out_csv)
  cat(sprintf("✓ Saved: %s\n", cfg$out_csv))

  style          <- build_cluster_style(3)
  cluster_cols   <- style$cols
  cluster_labels <- setNames(style$labels, names(style$cols))

  p_tpa <- ggplot(zonal, aes(x = factor(cluster), y = pct_with_mortality,
                              fill = factor(cluster))) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%\n(%.1f TPA)",
                                  pct_with_mortality, mean_tpa_all)),
              vjust = -0.3, size = 3.5, fontface = "bold", color = "#1B2A1F") +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_x_discrete(labels = cluster_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
    labs(
      title = cfg$label,
      x = "Cluster", y = "% pixels with beetle mortality",
      caption = sprintf("Relative risk (Cluster 3 vs 1): %.1f×", rr_c3_vs_c1)
    ) +
    theme_bw() +
    theme(plot.title   = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
          axis.title.y = element_text(face = "bold", color = "#1B2A1F"),
          axis.title.x = element_text(face = "bold", color = "#1B2A1F", margin = margin(t = 10)),
          axis.text    = element_text(face = "bold", size = 10, color = "#1B2A1F"),
          plot.caption = element_text(hjust = 0.5, size = 14, face = "bold", color = "#1B2A1F",
                                       margin = margin(t = 14)))

  p_tpa_ratio <- ggplot(zonal, aes(x = factor(cluster), y = ratio_to_c1,
                              fill = factor(cluster))) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.2f×", ratio_to_c1)),
              vjust = -0.3, size = 3.5, fontface = "bold", color = "#1B2A1F") +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_x_discrete(labels = cluster_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = cfg$label,
      x = "Cluster", y = "% pixels w/ beetle mortality, ratio vs. Cluster 1"
    ) +
    theme_bw() +
    theme(plot.title   = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
          axis.title.y = element_text(face = "bold", color = "#1B2A1F"),
          axis.title.x = element_text(face = "bold", color = "#1B2A1F", margin = margin(t = 10)),
          axis.text    = element_text(face = "bold", size = 10, color = "#1B2A1F"))

  list(zonal = zonal, plot = p_tpa, plot_ratio = p_tpa_ratio)
}

results <- lapply(index_config, run_beetle_validation, beetle_flags = beetle_flags)

p_combined <- (results$ndvi$plot + results$ndwi$plot) +
  plot_annotation(
    title    = sprintf("USFS ADS Beetle Mortality by CSD Cluster (%s)", year_label),
    subtitle = "% of pixels with bark beetle mortality (mean TPA shown for reference) — ground-truth evaluation, after cluster assignment",
    theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
                      plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"))
  )

ggsave(OUT_CHART, plot = p_combined, width = 12, height = 5, dpi = 150)
cat(sprintf("✓ Saved: %s\n", OUT_CHART))

p_combined_ratio <- (results$ndvi$plot_ratio + results$ndwi$plot_ratio) +
  plot_annotation(
    title    = sprintf("USFS ADS Beetle Mortality by CSD Cluster (%s)", year_label),
    subtitle = "% of pixels with bark beetle mortality, as a ratio vs. Cluster 1 (least-stressed) — ground-truth evaluation, after cluster assignment",
    theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
                      plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"))
  )

ggsave(OUT_CHART_RATIO, plot = p_combined_ratio, width = 12, height = 5, dpi = 150)
cat(sprintf("✓ Saved: %s\n", OUT_CHART_RATIO))

# Individual per-index panels, saved separately so 11_dashboard.R can
# assemble a one-page-per-index summary without duplicating this plotting code.
for (idx in names(index_config)) {
  panel_path <- sprintf("pixel_approach/data/processed/panel_beetle_%s.rds", idx)
  saveRDS(results[[idx]]$plot, panel_path)
  cat(sprintf("✓ Saved: %s\n", panel_path))

  panel_path_ratio <- sprintf("pixel_approach/data/processed/panel_beetle_%s_ratio.rds", idx)
  saveRDS(results[[idx]]$plot_ratio, panel_path_ratio)
  cat(sprintf("✓ Saved: %s\n", panel_path_ratio))
}