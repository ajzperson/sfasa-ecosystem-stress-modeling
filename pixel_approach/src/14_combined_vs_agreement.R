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
      !in_combined_top & in_agreement_both  ~ "\"NDVI_3 ∩ NDWI_3\" only",
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

group_order <- c("Neither", "\"NDVI+NDWI_4\" only", "\"NDVI_3 ∩ NDWI_3\" only", "Both")
summary_tbl <- d |>
  mutate(group = factor(group, levels = group_order)) |>
  group_by(group) |>
  summarise(
    n                 = n(),
    pct_beetle_exact  = 100 * mean(any_mortality, na.rm = TRUE),
    pct_beetle        = round(pct_beetle_exact, 1),
    pct_fire_exact    = 100 * mean(burned_validation, na.rm = TRUE),
    pct_fire          = round(pct_fire_exact, 1),
    .groups = "drop"
  ) |>
  arrange(group)

cat("\n=== Real-world 2015-2016 damage rates by group ===\n")
print(summary_tbl)

# Computed from the unrounded rate (pct_*_exact), not the 1-decimal display
# value -- dividing already-rounded percentages compounds their rounding
# error into the ratio.
rr_fire   <- summary_tbl$pct_fire_exact[summary_tbl$group == "Both"] /
  summary_tbl$pct_fire_exact[summary_tbl$group == "Neither"]
rr_beetle <- summary_tbl$pct_beetle_exact[summary_tbl$group == "Both"] /
  summary_tbl$pct_beetle_exact[summary_tbl$group == "Neither"]
cat(sprintf("  Relative risk (Both vs Neither): Fire %.2fx, Beetle %.2fx\n", rr_fire, rr_beetle))

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
            position = position_dodge(width = 0.75), vjust = -0.4, size = 3, fontface = "bold",
            color = "#1B2A1F") +
  scale_fill_manual(values = c(Beetle = "#A0522D", Fire = "#ED6623"), name = NULL) +
  scale_x_discrete(labels = setNames(group_x_labels, group_order)) +
  labs(
    title   = "Real 2015-2016 damage rates",
    caption = sprintf("**Both** = (NDVI_3 ∩ NDWI_3) ∩ NDVI+NDWI_4<br>**Neither** = in neither (NDVI_3 ∩ NDWI_3) nor NDVI+NDWI_4<br><span style='font-size:18pt'>**Relative risk (Both vs Neither): Fire %.1f×, Beetle %.1f×**</span>", rr_fire, rr_beetle),
    x = NULL, y = "% of pixels"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
    plot.caption    = element_markdown(hjust = 0.5, color = "#1B2A1F", size = 13, face = "bold",
                                        lineheight = 1.8, margin = margin(t = 14)),
    axis.title      = element_text(face = "bold", color = "#1B2A1F"),
    axis.text       = element_text(face = "bold", size = 10, color = "#1B2A1F"),
    legend.text     = element_text(face = "bold", size = 10, color = "#1B2A1F"),
    legend.position = "top"
  )

# --- Two circles above the bar chart: each shows what ONE definition of
# "most-stressed" looks like on its own, before the bar chart crosses them
# together. Left = the alone-indices definition (NDVI_3 & NDWI_3 overlap);
# right = the combined-clustering definition (NDVI+NDWI_4). Same total
# population (n = nrow(d)) in both, so the two circles are directly
# comparable in size.
build_pie <- function(n_flagged, n_total, flagged_label, wedge_color, color_label) {
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
               color = "#1B2A1F") +
    # Red family -- #D73027 (the "most-stressed" red used for NDVI+NDWI_4
    # everywhere else in the deck, slides 12-14) is reserved for that one
    # wedge only; NDVI_3 ∩ NDWI_3 gets a lighter tint of the same red so the
    # two pies read as one consistent color family without reusing the
    # reserved hue.
    scale_fill_manual(values = setNames(c(wedge_color, "grey85"), c(flagged_label, "Rest")), guide = "none") +
    labs(subtitle = sprintf("%s = **%s**", color_label, flagged_label)) +
    theme_void() +
    theme(plot.subtitle = element_markdown(hjust = 0.5, color = wedge_color, size = 10, face = "bold"))
}

n_total <- nrow(d)
pie_alone    <- build_pie(n_b, n_total, "NDVI_3 ∩ NDWI_3", "#E2655E", "Light red")
pie_combined <- build_pie(n_a, n_total, "NDVI+NDWI_4",     "#D73027", "Red")

# Small bordered notes flanking the two circles, spelling out what the
# cluster-number notation means -- "3" and "4" refer to different K's
# (separate clustering is K=3, combined is K=4), so both are each
# clustering's own top/most-stressed tier, not directly comparable numbers.
# geom_textbox (ggtext) instead of plain annotate("text", ...) -- it auto-wraps
# to the given width and understands markdown, so the term itself (NDVI_3,
# NDWI_3, NDVI+NDWI_4) can be marked with **term**; fontface="bold" below
# makes the whole block (term and definition alike) bold, matching every
# other text element in the deck.
build_note <- function(text_md) {
  ggplot() +
    geom_textbox(aes(x = 0.5, y = 0.5, label = text_md),
                 width = unit(1.5, "in"), box.color = NA, fill = NA,
                 hjust = 0.5, vjust = 0.5, halign = 0.5, size = 3.3,
                 fontface = "bold", color = "#1B2A1F") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.6),
      plot.margin  = margin(10, 8, 10, 8)
    )
}

wrap_def <- function(term, rhs) {
  sprintf("**%s = %s**", term, rhs)
}

note_left <- build_note(paste(
  wrap_def("NDVI_3", "NDVI Most Stressed Cluster"),
  wrap_def("NDWI_3", "NDWI Most Stressed Cluster"),
  sep = "<br><br>"
))
note_right <- build_note(
  wrap_def("NDVI+NDWI_4", "NDVI+NDWI-combined Most Stressed Cluster")
)

p_full <- (note_left | pie_alone | pie_combined | note_right) / p +
  plot_layout(heights = c(1, 1.6), widths = c(0.8, 1.4, 1.4, 0.8))

ggsave("pixel_approach/charts/combined_vs_agreement.pdf", plot = p_full, width = 9, height = 9.5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_vs_agreement.pdf\n")

# Also saved as PNG: the "∩" in the group labels above is outside the PDF
# device's built-in font encoding (falls back to "." when embedded as a
# vector font), but rasterizes correctly through the PNG device -- this PNG
# is the version actually used for the slide embed.
ggsave("pixel_approach/charts/combined_vs_agreement.png", plot = p_full, width = 9, height = 9.5, dpi = 300)
cat("✓ Saved: pixel_approach/charts/combined_vs_agreement.png\n")
