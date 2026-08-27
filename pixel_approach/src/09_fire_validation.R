library(tidyverse)
library(patchwork)

source("pixel_approach/src/clustering_helpers.R")

# Cross-validates K-means clusters against MTBS burn exposure (extracted in
# 02c_pixel_fire.R, independent of clustering). A pixel counts as "burned"
# only at moderate+ MTBS severity -- see 02c for the extraction details.
#
# Two checks, same fire data, different time windows:
#   1. Baseline window (2000-2014): informational only (clustering no longer
#      excludes on this -- see 06_clustering.R revert). Flags pixels whose
#      pre-collapse tau_r1/tau_var statistics may be fire-contaminated.
#   2. Validation window (VALIDATION_YEARS below, default 2015-2016, not
#      the full 2015-2019 MTBS availability): a second, independent
#      ground-truth outcome, same role as beetle mortality (08_beetle_
#      validation.R) -- default matches beetle's window (whose LEGACY_TPA
#      data effectively ends in 2016) so all three validators are compared
#      on equal footing. pixel_fire.csv now has one column per year
#      (burned_2015..burned_2019, from 02c_pixel_fire.R), so this window
#      can be widened/narrowed here without re-running extraction.
#
# Output:
#   pixel_approach/charts/fire_by_cluster.pdf (NDVI + NDWI side by side)
#   pixel_approach/data/processed/pixel_fire_ndvi.csv
#   pixel_approach/data/processed/pixel_fire_ndwi.csv

VALIDATION_YEARS <- 2015:2016
year_label <- if (length(VALIDATION_YEARS) == 1) as.character(VALIDATION_YEARS) else
  sprintf("%d-%d", min(VALIDATION_YEARS), max(VALIDATION_YEARS))

OUT_CHART <- "pixel_approach/charts/fire_by_cluster.pdf"

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
      n_pixels                  = n(),
      pct_burned_baseline_00_14 = round(100 * mean(burned_baseline_00_14), 1),
      pct_burned_validation     = round(100 * mean(burned_validation_15_16), 1),
      .groups = "drop"
    ) |>
    arrange(as.integer(cluster))

  cat(sprintf("\n=== %s: Fire Exposure by Cluster (MTBS moderate+ severity) ===\n", cfg$label))
  print(zonal)

  rr_c3_vs_c1 <- zonal$pct_burned_validation[zonal$cluster == 3] /
    zonal$pct_burned_validation[zonal$cluster == 1]
  cat(sprintf("  Relative risk (Cluster 3 vs 1): %.2fx\n", rr_c3_vs_c1))

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

  list(zonal = zonal, plot = p_fire)
}

results <- lapply(index_config, run_fire_validation, fire_flags = fire_flags)

p_combined <- (results$ndvi$plot + results$ndwi$plot) +
  plot_annotation(
    title    = sprintf("MTBS Moderate+ Burn Severity by CSD Cluster (%s)", year_label),
    subtitle = "% of pixels burned (severity >= moderate) — ground-truth validation, after cluster assignment",
    theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
                      plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"))
  )

ggsave(OUT_CHART, plot = p_combined, width = 12, height = 5, dpi = 150)
cat(sprintf("✓ Saved: %s\n", OUT_CHART))

# Individual per-index panels, saved separately so 11_dashboard.R can
# assemble a one-page-per-index summary without duplicating this plotting code.
for (idx in names(index_config)) {
  panel_path <- sprintf("pixel_approach/data/processed/panel_fire_%s.rds", idx)
  saveRDS(results[[idx]]$plot, panel_path)
  cat(sprintf("✓ Saved: %s\n", panel_path))
}