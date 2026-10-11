library(tidyverse)

# Two-circle Venn of NDVI vs NDWI most-stressed pixels; prints the group rates.
#
# Output: pixel_approach/charts/agreement_venn.pdf

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

circle_pts <- function(cx, cy, r, n = 200) {
  theta <- seq(0, 2 * pi, length.out = n)
  tibble(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

r <- 1.35
offset <- 0.7
circ_ndvi <- circle_pts(-offset, 0, r)
circ_ndwi <- circle_pts( offset, 0, r)

# Same NDVI/NDWI colors as 15_overlap_venn.R
p <- ggplot() +
  geom_polygon(data = circ_ndvi, aes(x, y), fill = "#E2655E", alpha = 0.95, color = "black", linewidth = 0.5) +
  geom_polygon(data = circ_ndwi, aes(x, y), fill = "#EC9B97", alpha = 0.95, color = "black", linewidth = 0.5) +
  annotate("text", x = -1.3, y = 0, label = "NDVI-3 only", size = 4, fontface = "bold", color = "white") +
  annotate("text", x = 0,    y = 0, label = "Both",        size = 5, fontface = "bold", color = "white") +
  annotate("text", x = 1.3,  y = 0, label = "NDWI-3 only", size = 4, fontface = "bold", color = "#1B2A1F") +
  coord_fixed(clip = "off") +
  xlim(-2.3, 2.3) + ylim(-1.6, 1.6) +
  theme_void()

ggsave("pixel_approach/charts/agreement_venn.pdf", plot = p, width = 4.6, height = 3.2, dpi = 150)
cat("✓ Saved: pixel_approach/charts/agreement_venn.pdf\n")
