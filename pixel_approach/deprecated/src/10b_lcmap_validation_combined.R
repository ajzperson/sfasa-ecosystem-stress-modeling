library(tidyverse)

# Cross-validates the COMBINED (NDVI+NDWI 4-feature) clustering (06b) against
# LCMAP tree-cover-loss transitions (02d) -- same role as
# 10_lcmap_validation.R, but for the combined-feature clustering approach
# instead of per-index clustering. Only one cluster assignment per pixel
# now, so one chart, not two panels.
#
# Output:
#   pixel_approach/charts/combined_lcmap_by_cluster.pdf
#   pixel_approach/data/processed/combined_pixel_lcmap.csv

VALIDATION_YEARS <- 2015:2016
year_label <- if (length(VALIDATION_YEARS) == 1) as.character(VALIDATION_YEARS) else
  sprintf("%d-%d", min(VALIDATION_YEARS), max(VALIDATION_YEARS))

mortality_flags <- read_csv("pixel_approach/data/processed/pixel_tree_mortality.csv",
                             show_col_types = FALSE)
year_cols <- sprintf("lost_tree_cover_%d", VALIDATION_YEARS)
stopifnot(all(year_cols %in% names(mortality_flags)))
mortality_flags$lost_tree_cover_15_16 <- rowSums(mortality_flags[year_cols]) > 0

pixels <- read_csv("pixel_approach/data/processed/combined_pixel_clusters.csv", show_col_types = FALSE) |>
  select(pixel_id, lon, lat, cluster, base_level_ndvi, base_level_ndwi) |>
  left_join(mortality_flags, by = "pixel_id")

cat(sprintf("  Pixels with tree cover loss (%s): %d / %d (%.1f%%)\n",
            year_label, sum(pixels$lost_tree_cover_15_16), nrow(pixels),
            100 * mean(pixels$lost_tree_cover_15_16)))

# Raw pixel-level tree-cover-loss detection is strongly driven by how much
# forest a pixel had to begin with (a pixel with little tree cover has few
# ways to trigger the LCACHG transition flag) -- verified: correlates 0.24
# with base_level_ndvi. That's the same density confound already removed
# from tau_r1/tau_skew via residualization, but it was never applied to
# this ground-truth outcome. Residualizing lost_tree_cover against both
# indices' baseline density isolates the tree-loss signal not explained by
# starting density, so it can be compared across clusters on equal footing.
dens_fit <- lm(as.numeric(lost_tree_cover_15_16) ~ base_level_ndvi + base_level_ndwi,
               data = pixels)
pixels$resid_lost_tree_cover <- resid(dens_fit)

zonal <- pixels |>
  group_by(cluster) |>
  summarise(
    n_pixels                 = n(),
    pct_lost_tree_cover      = round(100 * mean(lost_tree_cover_15_16), 1),
    mean_resid_lost_tree_cover = round(mean(resid_lost_tree_cover), 4),
    .groups = "drop"
  ) |>
  arrange(as.integer(cluster))

cat("\n=== Combined: LCMAP Tree Cover Loss by Cluster (raw vs density-adjusted) ===\n")
print(zonal)

write_csv(pixels, "pixel_approach/data/processed/combined_pixel_lcmap.csv")
cat("✓ Saved: pixel_approach/data/processed/combined_pixel_lcmap.csv\n")

K <- max(pixels$cluster, na.rm = TRUE)
cluster_cols <- setNames(
  colorRampPalette(c("#2166AC", "#92C5DE", "#F4A582", "#D73027"))(K),
  as.character(seq_len(K))
)
cluster_labels <- setNames(
  paste0(seq_len(K), c(" — no CSD", rep(" — partial", max(K - 2, 0)), " — full CSD")),
  as.character(seq_len(K))
)

p_lcmap <- ggplot(zonal, aes(x = factor(cluster), y = pct_lost_tree_cover,
                              fill = factor(cluster))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(adj %+.4f)",
                                pct_lost_tree_cover, mean_resid_lost_tree_cover)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = cluster_cols, guide = "none") +
  scale_x_discrete(labels = cluster_labels) +
  labs(
    title    = sprintf("LCMAP Tree Cover Loss by Combined Cluster (%s)", year_label),
    subtitle = "% of pixels with a tree-cover-loss transition (LCACHG); adj = mean residual after removing the\nbaseline-density confound — ground-truth validation, after cluster assignment",
    x = "Cluster", y = "% pixels with tree cover loss"
  ) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))

ggsave("pixel_approach/charts/combined_lcmap_by_cluster.pdf", plot = p_lcmap, width = 7, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_lcmap_by_cluster.pdf\n")
