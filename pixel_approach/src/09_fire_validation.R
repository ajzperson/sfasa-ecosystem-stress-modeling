library(tidyverse)
library(patchwork)

source("pixel_approach/src/clustering_helpers.R")

# Cross-validates K-means clusters against MTBS burn exposure (extracted in
# 02c_pixel_fire.R, independent of clustering). A pixel counts as "burned"
# only at moderate+ MTBS severity -- see 02c for the extraction details.
#
# Two checks, same fire data, different time windows:
#   1. Baseline window (2000-2014): a sanity check, not a real cross-tab --
#      06_clustering.R already excludes any pixel burned in this window
#      before clustering, so every cluster's rate here should read ~0%.
#   2. Validation window (VALIDATION_YEARS below, default 2015-2016): the
#      real ground-truth outcome, same role as beetle mortality
#      (08_beetle_validation.R) -- default matches beetle's window so both
#      validators are compared on equal footing. pixel_fire.csv has
#      one column per year (burned_2015..burned_2019), so this window can
#      be widened/narrowed here without re-running extraction.
#
# Output:
#   pixel_approach/charts/fire_by_cluster.pdf (NDVI + NDWI side by side)
#   pixel_approach/charts/fire_by_cluster_ratio.pdf (same, bars as ratio vs.
#     Cluster 1 instead of raw %) -- generated in parallel, raw chart/panels
#     untouched
#   pixel_approach/data/processed/pixel_fire_ndvi.csv
#   pixel_approach/data/processed/pixel_fire_ndwi.csv

VALIDATION_YEARS <- 2015:2016
year_label <- if (length(VALIDATION_YEARS) == 1) as.character(VALIDATION_YEARS) else
  sprintf("%d-%d", min(VALIDATION_YEARS), max(VALIDATION_YEARS))

OUT_CHART       <- "pixel_approach/charts/fire_by_cluster.pdf"
OUT_CHART_RATIO <- "pixel_approach/charts/fire_by_cluster_ratio.pdf"

index_config <- list(
  ndvi = list(pixel_path = "pixel_approach/data/processed/ndvi_pixel_clusters.csv",
              label      = "NDVI",
              out_csv    = "pixel_approach/data/processed/pixel_fire_ndvi.csv"),
  ndwi = list(pixel_path = "pixel_approach/data/processed/ndwi_pixel_clusters.csv",
              label      = "NDWI",
              out_csv    = "pixel_approach/data/processed/pixel_fire_ndwi.csv")
)

fire_flags <- read_csv("pixel_approach/data/processed/pixel_fire.csv", show_col_types = FALSE)

# OR together whichever per-year burned_<year> columns fall in
# VALIDATION_YEARS into one combined validation-window flag.
year_cols <- sprintf("burned_%d", VALIDATION_YEARS)
stopifnot(all(year_cols %in% names(fire_flags)))
fire_flags$burned_validation_15_16 <- rowSums(fire_flags[year_cols]) > 0

run_fire_validation <- function(cfg, fire_flags) {
  cat(sprintf("\n=== %s ===\n", cfg$label))
  pixels <- read_csv(cfg$pixel_path, show_col_types = FALSE) |>
    select(pixel_id, lon, lat, cluster) |>
    left_join(fire_flags, by = "pixel_id")

  cat(sprintf("  Pixels burned (moderate+) in baseline window (2000-2014, potential confound): %d / %d (%.1f%%)\n",
              sum(pixels$burned_baseline_00_14), nrow(pixels),
              100 * mean(pixels$burned_baseline_00_14)))
  cat(sprintf("  Pixels burned (moderate+) in validation window (%s): %d / %d (%.1f%%)\n",
              year_label, sum(pixels$burned_validation_15_16), nrow(pixels),
              100 * mean(pixels$burned_validation_15_16)))

  zonal <- pixels |>
    group_by(cluster) |>
    summarise(
      n_pixels                    = n(),
      pct_burned_baseline_00_14   = round(100 * mean(burned_baseline_00_14), 1),
      pct_burned_validation_exact = 100 * mean(burned_validation_15_16),
      pct_burned_validation       = round(pct_burned_validation_exact, 1),
      .groups = "drop"
    ) |>
    arrange(as.integer(cluster))

  cat(sprintf("\n=== %s: Fire Exposure by Cluster (MTBS moderate+ severity) ===\n", cfg$label))
  print(zonal)

  # Relative risk and ratio_to_c1 are computed from the unrounded rate
  # (pct_burned_validation_exact), not the 1-decimal display value --
  # dividing already-rounded percentages compounds their rounding error
  # into the ratio, and fire rates are small enough (~2-3%) that this is
  # not negligible (e.g. 1.8%/3.3% rounded reads 1.83x vs. the true 1.79x).
  rr_c3_vs_c1 <- zonal$pct_burned_validation_exact[zonal$cluster == 3] /
    zonal$pct_burned_validation_exact[zonal$cluster == 1]
  cat(sprintf("  Relative risk (Cluster 3 vs 1): %.2fx\n", rr_c3_vs_c1))

  # Ratio vs. Cluster 1 (least-stressed), for the ratio-chart variant below --
  # Cluster 1 is 1.00x by construction.
  zonal <- zonal |>
    mutate(ratio_to_c1 = pct_burned_validation_exact / pct_burned_validation_exact[cluster == 1])

  write_csv(pixels, cfg$out_csv)
  cat(sprintf("✓ Saved: %s\n", cfg$out_csv))

  style          <- build_cluster_style(3)
  cluster_cols   <- style$cols
  cluster_labels <- setNames(style$labels, names(style$cols))

  p_fire <- ggplot(zonal, aes(x = factor(cluster), y = pct_burned_validation,
                              fill = factor(cluster))) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", pct_burned_validation)),
              vjust = -0.3, size = 3.5, fontface = "bold", color = "#1B2A1F") +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_x_discrete(labels = cluster_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = cfg$label,
      x = "Cluster", y = "% pixels burned",
      caption = sprintf("Relative risk (Cluster 3 vs 1): %.1f×", rr_c3_vs_c1)
    ) +
    theme_bw() +
    theme(plot.title   = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
          axis.title.y = element_text(face = "bold", color = "#1B2A1F"),
          axis.title.x = element_text(face = "bold", color = "#1B2A1F", margin = margin(t = 10)),
          axis.text    = element_text(face = "bold", size = 10, color = "#1B2A1F"),
          plot.caption = element_text(hjust = 0.5, size = 14, face = "bold", color = "#1B2A1F",
                                       margin = margin(t = 14)))

  p_fire_ratio <- ggplot(zonal, aes(x = factor(cluster), y = ratio_to_c1,
                                    fill = factor(cluster))) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.2f×", ratio_to_c1)),
              vjust = -0.3, size = 3.5, fontface = "bold", color = "#1B2A1F") +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_x_discrete(labels = cluster_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = cfg$label,
      x = "Cluster", y = "% pixels burned, ratio vs. Cluster 1"
    ) +
    theme_bw() +
    theme(plot.title   = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
          axis.title.y = element_text(face = "bold", color = "#1B2A1F"),
          axis.title.x = element_text(face = "bold", color = "#1B2A1F", margin = margin(t = 10)),
          axis.text    = element_text(face = "bold", size = 10, color = "#1B2A1F"))

  list(zonal = zonal, plot = p_fire, plot_ratio = p_fire_ratio)
}

results <- lapply(index_config, run_fire_validation, fire_flags = fire_flags)

p_combined <- (results$ndvi$plot + results$ndwi$plot) +
  plot_annotation(
    title    = sprintf("MTBS Moderate+ Burn Severity by CSD Cluster (%s)", year_label),
    subtitle = "% of pixels burned (severity >= moderate) — ground-truth evaluation, after cluster assignment",
    theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
                      plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"))
  )

ggsave(OUT_CHART, plot = p_combined, width = 12, height = 5, dpi = 150)
cat(sprintf("✓ Saved: %s\n", OUT_CHART))

p_combined_ratio <- (results$ndvi$plot_ratio + results$ndwi$plot_ratio) +
  plot_annotation(
    title    = sprintf("MTBS Moderate+ Burn Severity by CSD Cluster (%s)", year_label),
    subtitle = "% of pixels burned (severity >= moderate), as a ratio vs. Cluster 1 (least-stressed) — ground-truth evaluation, after cluster assignment",
    theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
                      plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"))
  )

ggsave(OUT_CHART_RATIO, plot = p_combined_ratio, width = 12, height = 5, dpi = 150)
cat(sprintf("✓ Saved: %s\n", OUT_CHART_RATIO))

# Individual per-index panels, saved separately so 11_dashboard.R can
# assemble a one-page-per-index summary without duplicating this plotting code.
for (idx in names(index_config)) {
  panel_path <- sprintf("pixel_approach/data/processed/panel_fire_%s.rds", idx)
  saveRDS(results[[idx]]$plot, panel_path)
  cat(sprintf("✓ Saved: %s\n", panel_path))

  panel_path_ratio <- sprintf("pixel_approach/data/processed/panel_fire_%s_ratio.rds", idx)
  saveRDS(results[[idx]]$plot_ratio, panel_path_ratio)
  cat(sprintf("✓ Saved: %s\n", panel_path_ratio))
}