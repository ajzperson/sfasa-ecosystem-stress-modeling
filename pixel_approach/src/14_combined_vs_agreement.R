library(tidyverse)
library(patchwork)
library(ggtext)

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
      in_combined_top  & in_agreement_both  ~ "Both",
      in_combined_top  & !in_agreement_both ~ "\"NDVI+NDWI_4\" only",
      !in_combined_top & in_agreement_both  ~ "\"NDVI_3 & NDWI_3\" only",
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

group_order <- c("Neither", "\"NDVI+NDWI_4\" only", "\"NDVI_3 & NDWI_3\" only", "Both")
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
  scale_fill_manual(values = c(Beetle = "#A0522D", Fire = "#ED6623"), name = NULL) +
  scale_x_discrete(labels = setNames(group_x_labels, group_order)) +
  labs(
    title    = "NDVI+NDWI_4 vs. NDVI_3 & NDWI_3",
    subtitle = "Real 2015-2016 damage rates for two different definitions of \"most-stressed\"",
    x = NULL, y = "% of pixels"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
    legend.position = "top"
  )

# --- Two circles above the bar chart: each shows what ONE definition of
# "most-stressed" looks like on its own, before the bar chart crosses them
# together. Left = the alone-indices definition (NDVI_3 & NDWI_3 overlap);
# right = the combined-clustering definition (NDVI+NDWI_4). Same total
# population (n = nrow(d)) in both, so the two circles are directly
# comparable in size.
build_pie <- function(n_flagged, n_total, flagged_label, panel_title) {
  pie_df <- tibble(
    slice = factor(c(flagged_label, "Rest"), levels = c(flagged_label, "Rest")),
    n     = c(n_flagged, n_total - n_flagged)
  ) |>
    mutate(pct = 100 * n / sum(n))

  # Only the short %/n go inside the wedge -- the descriptive slice name
  # (e.g. "NDVI_3 & NDWI_3") is too long to fit inside the smaller wedge
  # without overflowing past its edge, so it moves to the subtitle instead,
  # with the wedge's own color as the visual link back to which slice it names.
  ggplot(pie_df, aes(x = "", y = n, fill = slice)) +
    geom_col(width = 1, color = "white", linewidth = 0.6) +
    coord_polar(theta = "y", start = 0) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", pct, format(n, big.mark = ","))),
               position = position_stack(vjust = 0.5), size = 3.6, fontface = "bold",
               color = c("black", "grey30")) +
    # Red -- matches the "most-stressed = red" color used for every cluster
    # map/legend elsewhere in the deck (build_cluster_style()). Fire's bar
    # color was changed to orange instead, so this stays consistent with
    # the rest of the deck rather than being the one-off exception.
    scale_fill_manual(values = setNames(c("#D73027", "grey85"), c(flagged_label, "Rest")), guide = "none") +
    labs(title = panel_title,
         subtitle = sprintf("Red = **%s**", flagged_label)) +
    theme_void() +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_markdown(hjust = 0.5, color = "#D73027", size = 9.5))
}

n_total <- nrow(d)
pie_alone    <- build_pie(n_b, n_total, "NDVI_3 & NDWI_3", "Separate\nClustering\n(NDVI-alone, NDWI-alone)")
pie_combined <- build_pie(n_a, n_total, "NDVI+NDWI_4",     "Combined\nClustering\n(NDVI + NDWI together)")

# Small bordered notes flanking the two circles, spelling out what the
# cluster-number notation means -- "3" and "4" refer to different K's
# (separate clustering is K=3, combined is K=4), so both are each
# clustering's own top/most-stressed tier, not directly comparable numbers.
# geom_textbox (ggtext) instead of plain annotate("text", ...) -- it auto-wraps
# to the given width AND understands markdown, so the term itself (NDVI_3,
# NDWI_3, NDVI+NDWI_4) can be bolded with **term** while the rest of the
# definition stays regular weight, all within one text block.
build_note <- function(text_md) {
  ggplot() +
    geom_textbox(aes(x = 0.5, y = 0.5, label = text_md),
                 width = unit(1.5, "in"), box.color = NA, fill = NA,
                 hjust = 0.5, vjust = 0.5, halign = 0.5, size = 3.3) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),
      plot.margin  = margin(10, 8, 10, 8)
    )
}

wrap_def <- function(term, rhs, cluster_line) {
  sprintf("**%s** = %s<br>%s", term, rhs, cluster_line)
}

note_left <- build_note(paste(
  wrap_def("NDVI_3", "Separate Clustering (NDVI-alone)", "Cluster 3 (Most-stressed)"),
  wrap_def("NDWI_3", "Separate Clustering (NDWI-alone)", "Cluster 3 (Most-stressed)"),
  sep = "<br><br>"
))
note_right <- build_note(
  wrap_def("NDVI+NDWI_4", "Combined Clustering (NDVI+NDWI)", "Cluster 4 (Most-stressed)")
)

p_full <- (note_left | pie_alone | pie_combined | note_right) / p +
  plot_layout(heights = c(1, 1.6), widths = c(0.8, 1.4, 1.4, 0.8))

ggsave("pixel_approach/charts/combined_vs_agreement.pdf", plot = p_full, width = 9, height = 9.5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_vs_agreement.pdf\n")
