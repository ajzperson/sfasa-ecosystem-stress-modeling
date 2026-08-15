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
#   pixel_approach/data/processed/combined_pixel_beetle_tpa.csv
#   pixel_approach/data/processed/panel_beetle_combined.rds

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
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = cluster_cols, guide = "none") +
  scale_x_discrete(labels = cluster_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = "Cluster", y = "% pixels with beetle mortality") +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))

# Title/subtitle only on the standalone chart -- kept off the base p_tpa so
# the panel saved below stays clean for 11b_dashboard_combined.R (matches
# 08_beetle_validation.R, where the per-index panel only ever gets a short
# cfg$label title and the elaborate title/subtitle is added separately, at
# ggsave time, for that script's own two-panel combined chart).
p_tpa_standalone <- p_tpa +
  labs(title    = sprintf("USFS ADS Beetle Mortality by Combined Cluster (%s)", year_label),
       subtitle = "% of pixels with bark beetle mortality (mean TPA shown for reference) — ground-truth validation, after cluster assignment")

ggsave("pixel_approach/charts/combined_beetle_tpa_by_cluster.pdf", plot = p_tpa_standalone, width = 7, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_beetle_tpa_by_cluster.pdf\n")

# Saved separately so 11b_dashboard_combined.R can assemble a one-page
# summary without duplicating this plotting code.
saveRDS(p_tpa, "pixel_approach/data/processed/panel_beetle_combined.rds")
cat("✓ Saved: pixel_approach/data/processed/panel_beetle_combined.rds\n")
