library(tidyverse)

source("pixel_approach/src/clustering_helpers.R")

# Share of pixels burned in 2015-2016 in each combined cluster.
#
# Output: pixel_approach/charts/combined_fire_by_cluster{,_ratio}.pdf
#         pixel_approach/data/processed/combined_pixel_fire.csv
#         pixel_approach/data/processed/panel_fire_combined{,_ratio}.rds

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
    n_pixels                    = n(),
    pct_burned_baseline_00_14   = round(100 * mean(burned_baseline_00_14), 1),
    pct_burned_validation_exact = 100 * mean(burned_validation_15_16),
    pct_burned_validation       = round(pct_burned_validation_exact, 1),
    .groups = "drop"
  ) |>
  arrange(as.integer(cluster))

cat("\n=== Combined: Fire Exposure by Cluster (MTBS moderate+ severity) ===\n")
print(zonal)

# Ratios use unrounded rates
K <- max(zonal$cluster, na.rm = TRUE)
rr_worst_vs_best <- zonal$pct_burned_validation_exact[zonal$cluster == K] /
  zonal$pct_burned_validation_exact[zonal$cluster == 1]
cat(sprintf("  Relative risk (Cluster %d vs 1): %.2fx\n", K, rr_worst_vs_best))

zonal <- zonal |>
  mutate(ratio_to_c1 = pct_burned_validation_exact / pct_burned_validation_exact[cluster == 1])

write_csv(pixels, "pixel_approach/data/processed/combined_pixel_fire.csv")
cat("✓ Saved: pixel_approach/data/processed/combined_pixel_fire.csv\n")

K <- max(pixels$cluster, na.rm = TRUE)
style          <- build_cluster_style(K)
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
  labs(x = "Cluster", y = "% pixels burned",
       caption = sprintf("Relative risk (Cluster %d vs 1): %.1f×", K, rr_worst_vs_best)) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
        plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"),
        axis.title.y  = element_text(face = "bold", color = "#1B2A1F"),
        axis.title.x  = element_text(face = "bold", color = "#1B2A1F", margin = margin(t = 10)),
        axis.text     = element_text(face = "bold", size = 10, color = "#1B2A1F"),
        plot.caption  = element_text(hjust = 0.5, size = 14, face = "bold", color = "#1B2A1F",
                                      margin = margin(t = 14)))

# Title only on the standalone chart; the saved panel stays untitled for 11b
p_fire_standalone <- p_fire +
  labs(title    = sprintf("MTBS Moderate+ Burn Severity by Combined Cluster (%s)", year_label),
       subtitle = "% of pixels burned (severity >= moderate) — ground-truth evaluation, after cluster assignment")

ggsave("pixel_approach/charts/combined_fire_by_cluster.pdf", plot = p_fire_standalone, width = 7, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_fire_by_cluster.pdf\n")

p_fire_ratio <- ggplot(zonal, aes(x = factor(cluster), y = ratio_to_c1,
                                  fill = factor(cluster))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f×", ratio_to_c1)),
            vjust = -0.3, size = 3.5, fontface = "bold", color = "#1B2A1F") +
  scale_fill_manual(values = cluster_cols, guide = "none") +
  scale_x_discrete(labels = cluster_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Cluster", y = "% pixels burned, ratio vs. Cluster 1") +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
        plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"),
        axis.title.y  = element_text(face = "bold", color = "#1B2A1F"),
        axis.title.x  = element_text(face = "bold", color = "#1B2A1F", margin = margin(t = 10)),
        axis.text     = element_text(face = "bold", size = 10, color = "#1B2A1F"))

p_fire_ratio_standalone <- p_fire_ratio +
  labs(title    = sprintf("MTBS Moderate+ Burn Severity by Combined Cluster (%s)", year_label),
       subtitle = "% of pixels burned (severity >= moderate), as a ratio vs. Cluster 1 (least-stressed) — ground-truth evaluation, after cluster assignment")

ggsave("pixel_approach/charts/combined_fire_by_cluster_ratio.pdf", plot = p_fire_ratio_standalone, width = 7, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_fire_by_cluster_ratio.pdf\n")

# Panels for 11b_dashboard_combined.R
saveRDS(p_fire, "pixel_approach/data/processed/panel_fire_combined.rds")
cat("✓ Saved: pixel_approach/data/processed/panel_fire_combined.rds\n")

saveRDS(p_fire_ratio, "pixel_approach/data/processed/panel_fire_combined_ratio.rds")
cat("✓ Saved: pixel_approach/data/processed/panel_fire_combined_ratio.rds\n")
