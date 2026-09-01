library(tidyverse)

source("pixel_approach/src/clustering_helpers.R")

# Cross-validates the COMBINED (NDVI+NDWI 4-feature) clustering (06b) against
# USFS ADS beetle-kill mortality (02e) -- same role as
# 08_beetle_validation.R, but for the combined-feature clustering approach
# instead of per-index clustering. Only one cluster assignment per pixel
# now, so one chart, not two panels.
#
# Output:
#   pixel_approach/charts/combined_beetle_tpa_by_cluster.pdf
#   pixel_approach/charts/combined_beetle_tpa_by_cluster_ratio.pdf (same,
#     bars as ratio vs. Cluster 1 instead of raw %) -- generated in parallel,
#     raw chart/panel untouched
#   pixel_approach/data/processed/combined_pixel_beetle_tpa.csv
#   pixel_approach/data/processed/panel_beetle_combined.rds
#   pixel_approach/data/processed/panel_beetle_combined_ratio.rds

VALIDATION_YEARS <- 2015:2016
year_label <- if (length(VALIDATION_YEARS) == 1) as.character(VALIDATION_YEARS) else
  sprintf("%d-%d", min(VALIDATION_YEARS), max(VALIDATION_YEARS))

beetle_flags <- read_csv("pixel_approach/data/processed/pixel_beetle.csv", show_col_types = FALSE)
win_cols <- sprintf("tpa_%d", VALIDATION_YEARS)
stopifnot(all(win_cols %in% names(beetle_flags)))
beetle_flags$tpa <- do.call(pmax, c(unname(beetle_flags[win_cols]), na.rm = TRUE))
beetle_flags$any_mortality <- beetle_flags$tpa > 0

pixels <- read_csv("pixel_approach/data/processed/combined_pixel_clusters.csv", show_col_types = FALSE) |>
  select(pixel_id, lon, lat, cluster) |>
  left_join(beetle_flags, by = "pixel_id")

cat(sprintf("  Pixels with beetle mortality (%s): %d / %d (%.1f%%)\n",
            year_label, sum(pixels$any_mortality), nrow(pixels),
            100 * mean(pixels$any_mortality)))

zonal <- pixels |>
  group_by(cluster) |>
  summarise(
    n_pixels           = n(),
    n_with_mortality   = sum(any_mortality),
    pct_with_mortality = round(100 * mean(any_mortality), 1),
    mean_tpa_affected  = round(mean(tpa[any_mortality], na.rm = TRUE), 1),
    mean_tpa_all       = round(mean(tpa, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  arrange(as.integer(cluster))

cat("\n=== Combined: Zonal TPA by Cluster ===\n")
print(zonal)

K <- max(zonal$cluster, na.rm = TRUE)
rr_worst_vs_best <- zonal$pct_with_mortality[zonal$cluster == K] /
  zonal$pct_with_mortality[zonal$cluster == 1]
cat(sprintf("  Relative risk (Cluster %d vs 1): %.2fx\n", K, rr_worst_vs_best))

# Ratio vs. Cluster 1 (least-stressed), for the ratio-chart variant below --
# Cluster 1 is 1.00x by construction.
zonal <- zonal |>
  mutate(ratio_to_c1 = pct_with_mortality / pct_with_mortality[cluster == 1])

write_csv(pixels, "pixel_approach/data/processed/combined_pixel_beetle_tpa.csv")
cat("✓ Saved: pixel_approach/data/processed/combined_pixel_beetle_tpa.csv\n")

K <- max(pixels$cluster, na.rm = TRUE)
style          <- build_cluster_style(K)
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
  labs(x = "Cluster", y = "% pixels with beetle mortality",
       caption = sprintf("Relative risk (Cluster %d vs 1): %.1f×", K, rr_worst_vs_best)) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
        plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"),
        axis.title.y  = element_text(face = "bold", color = "#1B2A1F"),
        axis.title.x  = element_text(face = "bold", color = "#1B2A1F", margin = margin(t = 10)),
        axis.text     = element_text(face = "bold", size = 10, color = "#1B2A1F"),
        plot.caption  = element_text(hjust = 0.5, size = 14, face = "bold", color = "#1B2A1F",
                                      margin = margin(t = 14)))

# Title/subtitle only on the standalone chart -- kept off the base p_tpa so
# the panel saved below stays clean for 11b_dashboard_combined.R (matches
# 08_beetle_validation.R, where the per-index panel only ever gets a short
# cfg$label title and the elaborate title/subtitle is added separately, at
# ggsave time, for that script's own two-panel combined chart).
p_tpa_standalone <- p_tpa +
  labs(title    = sprintf("USFS ADS Beetle Mortality by Combined Cluster (%s)", year_label),
       # Wrapped -- at this chart's narrower standalone width (7in, vs the
       # two-panel version's 12in that this same subtitle was sized for),
       # the unwrapped string overflows past the right edge of the canvas.
       subtitle = str_wrap("% of pixels with bark beetle mortality (mean TPA shown for reference) — ground-truth validation, after cluster assignment", width = 70))

ggsave("pixel_approach/charts/combined_beetle_tpa_by_cluster.pdf", plot = p_tpa_standalone, width = 7, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_beetle_tpa_by_cluster.pdf\n")

p_tpa_ratio <- ggplot(zonal, aes(x = factor(cluster), y = ratio_to_c1,
                            fill = factor(cluster))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f×", ratio_to_c1)),
            vjust = -0.3, size = 3.5, fontface = "bold", color = "#1B2A1F") +
  scale_fill_manual(values = cluster_cols, guide = "none") +
  scale_x_discrete(labels = cluster_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Cluster", y = "% pixels w/ beetle mortality, ratio vs. Cluster 1") +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
        plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"),
        axis.title.y  = element_text(face = "bold", color = "#1B2A1F"),
        axis.title.x  = element_text(face = "bold", color = "#1B2A1F", margin = margin(t = 10)),
        axis.text     = element_text(face = "bold", size = 10, color = "#1B2A1F"))

p_tpa_ratio_standalone <- p_tpa_ratio +
  labs(title    = sprintf("USFS ADS Beetle Mortality by Combined Cluster (%s)", year_label),
       subtitle = str_wrap("% of pixels with bark beetle mortality, as a ratio vs. Cluster 1 (least-stressed) — ground-truth validation, after cluster assignment", width = 70))

ggsave("pixel_approach/charts/combined_beetle_tpa_by_cluster_ratio.pdf", plot = p_tpa_ratio_standalone, width = 7, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_beetle_tpa_by_cluster_ratio.pdf\n")

# Saved separately so 11b_dashboard_combined.R can assemble a one-page
# summary without duplicating this plotting code.
saveRDS(p_tpa, "pixel_approach/data/processed/panel_beetle_combined.rds")
cat("✓ Saved: pixel_approach/data/processed/panel_beetle_combined.rds\n")

saveRDS(p_tpa_ratio, "pixel_approach/data/processed/panel_beetle_combined_ratio.rds")
cat("✓ Saved: pixel_approach/data/processed/panel_beetle_combined_ratio.rds\n")
