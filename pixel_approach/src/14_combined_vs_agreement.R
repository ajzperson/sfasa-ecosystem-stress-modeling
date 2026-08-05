library(tidyverse)

# Compares the combined NDVI+NDWI clustering's "most stressed" group against
# the pixels where NDVI-alone and NDWI-alone's INDEPENDENT clusterings both
# separately agreed on "full CSD" (the "Both" group from
# 13_cluster_agreement_validation.R) -- tests whether the more complex
# combined 4-feature approach is rediscovering the same pixels a much
# simpler two-index agreement check already finds, or flagging a genuinely
# different set. Where they disagree, checks real 2015-2016 beetle/fire
# outcomes to see which method is more predictive in the disagreement zone.
#
# Output:
#   pixel_approach/charts/combined_vs_agreement.pdf
#   pixel_approach/data/processed/combined_vs_agreement.csv

combined <- read_csv("pixel_approach/data/processed/combined_pixel_clusters.csv", show_col_types = FALSE) |>
  select(pixel_id, combined_cluster = cluster)
K_combined <- max(as.integer(combined$combined_cluster))
cat(sprintf("Combined clustering K = %d\n", K_combined))

ndvi <- read_csv("pixel_approach/data/processed/ndvi_pixel_clusters.csv", show_col_types = FALSE) |>
  select(pixel_id, ndvi_cluster = cluster)
ndwi <- read_csv("pixel_approach/data/processed/ndwi_pixel_clusters.csv", show_col_types = FALSE) |>
  select(pixel_id, ndwi_cluster = cluster)

d <- combined |>
  inner_join(ndvi, by = "pixel_id") |>
  inner_join(ndwi, by = "pixel_id") |>
  mutate(
    in_combined_top   = as.integer(combined_cluster) == K_combined,
    in_agreement_both = as.integer(ndvi_cluster) == 3 & as.integer(ndwi_cluster) == 3,
    group = case_when(
      in_combined_top  & in_agreement_both  ~ "Both methods agree",
      in_combined_top  & !in_agreement_both ~ "Combined only",
      !in_combined_top & in_agreement_both  ~ "Agreement only",
      TRUE                                  ~ "Neither"
    )
  )
cat(sprintf("✓ %d pixels present in all three clusterings\n", nrow(d)))

cat("\n=== Overlap group sizes ===\n")
print(table(d$group))

n_a     <- sum(d$in_combined_top)
n_b     <- sum(d$in_agreement_both)
n_both  <- sum(d$in_combined_top & d$in_agreement_both)
n_union <- sum(d$in_combined_top | d$in_agreement_both)
jaccard <- n_both / n_union
cat(sprintf("\nCombined top cluster (cluster %d): %d pixels\n", K_combined, n_a))
cat(sprintf("Independent agreement (both alone flag full CSD): %d pixels\n", n_b))
cat(sprintf("Overlap: %d pixels\n", n_both))
cat(sprintf("Jaccard index (overlap / union): %.3f\n", jaccard))

# --- Real-world validation in the disagreement zone -------------------------
beetle <- read_csv("pixel_approach/data/processed/pixel_beetle.csv", show_col_types = FALSE)
beetle$any_mortality <- (!is.na(beetle$tpa_2015) & beetle$tpa_2015 > 0) |
                         (!is.na(beetle$tpa_2016) & beetle$tpa_2016 > 0)

fire <- read_csv("pixel_approach/data/processed/pixel_fire.csv", show_col_types = FALSE)
fire$burned_validation <- fire$burned_2015 | fire$burned_2016

d <- d |>
  left_join(beetle |> select(pixel_id, any_mortality), by = "pixel_id") |>
  left_join(fire |> select(pixel_id, burned_validation), by = "pixel_id")

group_order <- c("Neither", "Combined only", "Agreement only", "Both methods agree")
summary_tbl <- d |>
  mutate(group = factor(group, levels = group_order)) |>
  group_by(group) |>
  summarise(
    n           = n(),
    pct_beetle  = round(100 * mean(any_mortality, na.rm = TRUE), 1),
    pct_fire    = round(100 * mean(burned_validation, na.rm = TRUE), 1),
    .groups = "drop"
  ) |>
  arrange(group)

cat("\n=== Real-world 2015-2016 damage rates by group ===\n")
print(summary_tbl)

write_csv(d, "pixel_approach/data/processed/combined_vs_agreement.csv")
cat("\n✓ Saved: pixel_approach/data/processed/combined_vs_agreement.csv\n")

# --- Chart -------------------------------------------------------------
plot_df <- summary_tbl |>
  pivot_longer(cols = c(pct_beetle, pct_fire), names_to = "validator", values_to = "pct") |>
  mutate(validator = recode(validator, pct_beetle = "Beetle", pct_fire = "Fire"),
         validator = factor(validator, levels = c("Fire", "Beetle")))

group_n <- setNames(summary_tbl$n, as.character(summary_tbl$group))
group_x_labels <- sprintf("%s\n(n=%s)", group_order, format(group_n[group_order], big.mark = ","))

p <- ggplot(plot_df, aes(x = group, y = pct, fill = validator)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.75), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c(Beetle = "#4DAF4A", Fire = "#E41A1C"), name = NULL) +
  scale_x_discrete(labels = setNames(group_x_labels, group_order)) +
  labs(
    title    = "Combined Clustering vs. Independent NDVI/NDWI Agreement",
    subtitle = "Real 2015-2016 damage rates by whether each method flags a pixel as most-stressed",
    x = NULL, y = "% of pixels"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
    legend.position = "top"
  )

ggsave("pixel_approach/charts/combined_vs_agreement.pdf", plot = p, width = 9, height = 6, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_vs_agreement.pdf\n")
