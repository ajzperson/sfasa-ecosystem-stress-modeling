library(tidyverse)

# A 2-circle Venn diagram (NDVI Most-stressed, NDWI Most-stressed) for
# slide 15 ("Where NDVI and NDWI independently agree, the signal is
# strongest") -- same visual language as 15_overlap_venn.R (slide 16), but
# simpler: 13_cluster_agreement_validation.R only ever compares NDVI-alone
# vs NDWI-alone directly (no Combined circle, and "Neither" is deliberately
# dropped -- "not the question this script asks", per that script's own
# comment), so this Venn only needs the two circles plus their overlap, no
# G4/outer boundary.
#
# Colors are red-family tints, matching 15_overlap_venn.R's NDVI_3/NDWI_3
# tiers exactly (same clusters, so same color in both places) -- #D73027,
# the deck's "most-stressed" red from slides 12-14, is reserved for
# NDVI+NDWI_4 only and never appears here.
#
# Output:
#   pixel_approach/charts/agreement_venn.pdf
#   pixel_approach/charts/agreement_venn.png

d <- read_csv("pixel_approach/data/processed/combined_vs_agreement.csv", show_col_types = FALSE)

rate <- function(subset, col) {
  n <- nrow(subset)
  if (n == 0) return(c(n = 0, pct = 0))
  c(n = n, pct = round(100 * mean(subset[[col]] == TRUE, na.rm = TRUE), 1))
}

ndvi_only <- filter(d, ndvi_cluster == 3, ndwi_cluster != 3)
ndwi_only <- filter(d, ndvi_cluster != 3, ndwi_cluster == 3)
both      <- filter(d, ndvi_cluster == 3, ndwi_cluster == 3)

ndvi_n <- rate(ndvi_only, "any_mortality"); ndvi_fire <- rate(ndvi_only, "burned_validation")
ndwi_n <- rate(ndwi_only, "any_mortality"); ndwi_fire <- rate(ndwi_only, "burned_validation")
both_n <- rate(both, "any_mortality");      both_fire <- rate(both, "burned_validation")

cat(sprintf("NDVI-3 only n=%.0f beetle=%.1f%% fire=%.1f%%\n", ndvi_n["n"], ndvi_n["pct"], ndvi_fire["pct"]))
cat(sprintf("NDWI-3 only n=%.0f beetle=%.1f%% fire=%.1f%%\n", ndwi_n["n"], ndwi_n["pct"], ndwi_fire["pct"]))
cat(sprintf("Both n=%.0f beetle=%.1f%% fire=%.1f%%\n", both_n["n"], both_n["pct"], both_fire["pct"]))
cat(sprintf("NDVI-3 total = %.0f | NDWI-3 total = %.0f\n",
            ndvi_n["n"] + both_n["n"], ndwi_n["n"] + both_n["n"]))

# --- Circle geometry ---------------------------------------------------
circle_pts <- function(cx, cy, r, n = 200) {
  theta <- seq(0, 2 * pi, length.out = n)
  tibble(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

r <- 1.35
offset <- 0.7
circ_ndvi <- circle_pts(-offset, 0, r)
circ_ndwi <- circle_pts( offset, 0, r)

# No G1/G2/G3 shading system here (that's for the 3-circle case in
# 15_overlap_venn.R, which needs to rank four regions). With only two
# circles there's nothing to rank -- NDVI-3-only and NDWI-3-only should
# just be their plain base colors, with the natural alpha-overlap giving
# "Both" its own shade. A wash layer here would darken the base colors
# away from 15_overlap_venn.R's NDVI_3/NDWI_3 tiers, breaking the
# same-cluster-same-color match between the two venns.
p <- ggplot() +
  geom_polygon(data = circ_ndvi, aes(x, y), fill = "#E2655E", alpha = 0.95, color = "black", linewidth = 0.5) +
  geom_polygon(data = circ_ndwi, aes(x, y), fill = "#EC9B97", alpha = 0.95, color = "black", linewidth = 0.5) +
  annotate("text", x = -1.3, y = 0, label = "NDVI-3 only", size = 4, fontface = "bold", color = "white") +
  annotate("text", x = 0,    y = 0, label = "Both",        size = 5, fontface = "bold", color = "white") +
  annotate("text", x = 1.3,  y = 0, label = "NDWI-3 only", size = 4, fontface = "bold", color = "black") +
  coord_fixed(clip = "off") +
  xlim(-2.3, 2.3) + ylim(-1.6, 1.6) +
  theme_void()

ggsave("pixel_approach/charts/agreement_venn.pdf", plot = p, width = 4.6, height = 3.2, dpi = 150)
cat("✓ Saved: pixel_approach/charts/agreement_venn.pdf\n")

ggsave("pixel_approach/charts/agreement_venn.png", plot = p, width = 4.6, height = 3.2, dpi = 300)
cat("✓ Saved: pixel_approach/charts/agreement_venn.png\n")
