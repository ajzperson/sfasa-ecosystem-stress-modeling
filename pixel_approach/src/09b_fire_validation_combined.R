library(tidyverse)

# Cross-validates the COMBINED (NDVI+NDWI 4-feature) clustering (06b) against
# MTBS burn exposure (02c) -- same role as 09_fire_validation.R, but for the
# combined-feature clustering approach instead of per-index clustering.
# Only one cluster assignment per pixel now, so one chart, not two panels.
#
# Output:
#   pixel_approach/charts/combined_fire_by_cluster.pdf
#   pixel_approach/data/processed/combined_pixel_fire.csv

VALIDATION_YEARS <- 2015:2016
year_label <- if (length(VALIDATION_YEARS) == 1) as.character(VALIDATION_YEARS) else
  sprintf("%d-%d", min(VALIDATION_YEARS), max(VALIDATION_YEARS))

fire_flags <- read_csv("pixel_approach/data/processed/pixel_fire.csv", show_col_types = FALSE)
year_cols <- sprintf("burned_%d", VALIDATION_YEARS)
stopifnot(all(year_cols %in% names(fire_flags)))
fire_flags$burned_validation_15_16 <- rowSums(fire_flags[year_cols]) > 0

pixels <- read_csv("pixel_approach/data/processed/combined_pixel_clusters.csv", show_col_types = FALSE) |>
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

cat("\n=== Combined: Fire Exposure by Cluster (MTBS moderate+ severity) ===\n")
print(zonal)

write_csv(pixels, "pixel_approach/data/processed/combined_pixel_fire.csv")
cat("✓ Saved: pixel_approach/data/processed/combined_pixel_fire.csv\n")

K <- max(pixels$cluster, na.rm = TRUE)
cluster_cols <- setNames(
  colorRampPalette(c("#2166AC", "#92C5DE", "#F4A582", "#D73027"))(K),
  as.character(seq_len(K))
)
cluster_labels <- setNames(
  paste0(seq_len(K), c(" — no CSD", rep(" — partial", max(K - 2, 0)), " — full CSD")),
  as.character(seq_len(K))
)

p_fire <- ggplot(zonal, aes(x = factor(cluster), y = pct_burned_validation,
                            fill = factor(cluster))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct_burned_validation)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = cluster_cols, guide = "none") +
  scale_x_discrete(labels = cluster_labels) +
  labs(
    title    = sprintf("MTBS Moderate+ Burn Severity by Combined Cluster (%s)", year_label),
    subtitle = "% of pixels burned (severity >= moderate) — ground-truth validation, after cluster assignment",
    x = "Cluster", y = "% pixels burned"
  ) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))

ggsave("pixel_approach/charts/combined_fire_by_cluster.pdf", plot = p_fire, width = 7, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_fire_by_cluster.pdf\n")
