library(tidyverse)
library(patchwork)

# Cross-validates K-means clusters against LCMAP tree-cover-loss transitions
# (extracted in 02d_pixel_tree_mortality.R, independent of clustering) --
# an independent mortality signal alongside the bark beetle ADS data
# (08_beetle_validation.R), using a genuinely different detection
# methodology (Landsat CCDC spectral time-series break detection, not a
# vegetation-index trend or an aerial survey).
#
# Validation window (VALIDATION_YEARS below, default 2015-2016, not the
# full 2015-2019 LCMAP availability) matches beetle's ADS window by
# default; no baseline-confound check applies here the way it does for
# fire. pixel_tree_mortality.csv now has one column per year
# (lost_tree_cover_2015..lost_tree_cover_2019, from 02d_pixel_tree_
# mortality.R), so this window can be widened/narrowed here without
# re-running extraction.
#
# Output:
#   pixel_approach/charts/lcmap_by_cluster.pdf (NDVI + NDWI side by side)
#   pixel_approach/data/processed/pixel_lcmap_ndvi.csv
#   pixel_approach/data/processed/pixel_lcmap_ndwi.csv

VALIDATION_YEARS <- 2015:2016
year_label <- if (length(VALIDATION_YEARS) == 1) as.character(VALIDATION_YEARS) else
  sprintf("%d-%d", min(VALIDATION_YEARS), max(VALIDATION_YEARS))

OUT_CHART <- "pixel_approach/charts/lcmap_by_cluster.pdf"

index_config <- list(
  ndvi = list(pixel_path = "pixel_approach/data/processed/ndvi_pixel_clusters.csv",
              label      = "NDVI",
              out_csv    = "pixel_approach/data/processed/pixel_lcmap_ndvi.csv"),
  ndwi = list(pixel_path = "pixel_approach/data/processed/ndwi_pixel_clusters.csv",
              label      = "NDWI",
              out_csv    = "pixel_approach/data/processed/pixel_lcmap_ndwi.csv")
)

mortality_flags <- read_csv("pixel_approach/data/processed/pixel_tree_mortality.csv",
                             show_col_types = FALSE)

# OR together whichever per-year lost_tree_cover_<year> columns fall in
# VALIDATION_YEARS into one combined validation-window flag.
year_cols <- sprintf("lost_tree_cover_%d", VALIDATION_YEARS)
stopifnot(all(year_cols %in% names(mortality_flags)))
mortality_flags$lost_tree_cover_15_16 <- rowSums(mortality_flags[year_cols]) > 0

run_lcmap_validation <- function(cfg, mortality_flags) {
  cat(sprintf("\n=== %s ===\n", cfg$label))
  pixels <- read_csv(cfg$pixel_path, show_col_types = FALSE) |>
    select(pixel_id, lon, lat, cluster) |>
    left_join(mortality_flags, by = "pixel_id")

  cat(sprintf("  Pixels with tree cover loss (%s): %d / %d (%.1f%%)\n",
              year_label, sum(pixels$lost_tree_cover_15_16), nrow(pixels),
              100 * mean(pixels$lost_tree_cover_15_16)))

  zonal <- pixels |>
    group_by(cluster) |>
    summarise(
      n_pixels                = n(),
      pct_lost_tree_cover     = round(100 * mean(lost_tree_cover_15_16), 1),
      .groups = "drop"
    ) |>
    arrange(as.integer(cluster))

  cat(sprintf("\n=== %s: LCMAP Tree Cover Loss by Cluster ===\n", cfg$label))
  print(zonal)

  write_csv(pixels, cfg$out_csv)
  cat(sprintf("✓ Saved: %s\n", cfg$out_csv))

  cluster_cols <- setNames(
    colorRampPalette(c("#2166AC", "#92C5DE", "#F4A582", "#D73027"))(3),
    c("1", "2", "3")
  )
  cluster_labels <- c("1" = "1 — no CSD", "2" = "2 — partial", "3" = "3 — full CSD")

  p_lcmap <- ggplot(zonal, aes(x = factor(cluster), y = pct_lost_tree_cover,
                                fill = factor(cluster))) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", pct_lost_tree_cover)),
              vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_x_discrete(labels = cluster_labels) +
    labs(
      title = cfg$label,
      x = "Cluster", y = "% pixels with tree cover loss"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  list(zonal = zonal, plot = p_lcmap)
}

results <- lapply(index_config, run_lcmap_validation, mortality_flags = mortality_flags)

p_combined <- (results$ndvi$plot + results$ndwi$plot) +
  plot_annotation(
    title    = sprintf("LCMAP Tree Cover Loss by CSD Cluster (%s)", year_label),
    subtitle = "% of pixels with a tree-cover-loss transition (LCACHG) — ground-truth validation, after cluster assignment",
    theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
                      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))
  )

ggsave(OUT_CHART, plot = p_combined, width = 12, height = 5, dpi = 150)
cat(sprintf("✓ Saved: %s\n", OUT_CHART))

# Individual per-index panels, saved separately so 11_dashboard.R can
# assemble a one-page-per-index summary without duplicating this plotting code.
for (idx in names(index_config)) {
  panel_path <- sprintf("pixel_approach/data/processed/panel_lcmap_%s.rds", idx)
  saveRDS(results[[idx]]$plot, panel_path)
  cat(sprintf("✓ Saved: %s\n", panel_path))
}
