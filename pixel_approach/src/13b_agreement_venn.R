library(tidyverse)
library(sf)

# A 2-circle Venn diagram (NDVI Most-stressed, NDWI Most-stressed) for
# slide 15 ("Where NDVI and NDWI independently agree, the signal is
# strongest") -- same visual language as 15_overlap_venn.R (slide 17), but
# simpler: 13_cluster_agreement_validation.R only ever compares NDVI-alone
# vs NDWI-alone directly (no Combined circle, and "Neither" is deliberately
# dropped -- "not the question this script asks", per that script's own
# comment), so this Venn only needs G1/G2/G3, no G4/outer boundary.
#
# Output:
#   pixel_approach/charts/agreement_venn.pdf

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

poly_A <- st_polygon(list(cbind(c(circ_ndvi$x, circ_ndvi$x[1]), c(circ_ndvi$y, circ_ndvi$y[1]))))
poly_B <- st_polygon(list(cbind(c(circ_ndwi$x, circ_ndwi$x[1]), c(circ_ndwi$y, circ_ndwi$y[1]))))

g1_geom <- st_intersection(poly_A, poly_B)  # Both -- strongest
g2_geom <- st_difference(poly_A, poly_B)    # NDVI-3 only -- 2nd strongest
g3_geom <- st_difference(poly_B, poly_A)    # NDWI-3 only -- 3rd strongest

to_df <- function(geom, group) {
  m <- st_coordinates(geom)
  tibble(x = m[, "X"], y = m[, "Y"], piece = if ("L2" %in% colnames(m)) m[, "L2"] else 1, group = group)
}
regions <- bind_rows(to_df(g2_geom, "G2"), to_df(g3_geom, "G3"), to_df(g1_geom, "G1"))

shade_alpha <- c(G1 = 0.42, G2 = 0.22, G3 = 0.12)

p <- ggplot() +
  geom_polygon(data = circ_ndvi, aes(x, y), fill = "#2166AC", alpha = 0.5, color = NA) +
  geom_polygon(data = circ_ndwi, aes(x, y), fill = "#92C5DE", alpha = 0.5, color = NA) +
  geom_polygon(data = regions, aes(x, y, group = interaction(group, piece), alpha = group),
               fill = "black", color = "black", linewidth = 0.35) +
  scale_alpha_manual(values = shade_alpha, guide = "none") +
  annotate("text", x = 0,        y = 0, label = "G1", size = 5, fontface = "bold", color = "white") +
  annotate("text", x = -1.3,     y = 0, label = "G2", size = 5, fontface = "bold", color = "black") +
  annotate("text", x = 1.3,      y = 0, label = "G3", size = 5, fontface = "bold", color = "black") +
  coord_fixed(clip = "off") +
  xlim(-2.3, 2.3) + ylim(-1.6, 1.6) +
  theme_void()

ggsave("pixel_approach/charts/agreement_venn.pdf", plot = p, width = 4.6, height = 3.2, dpi = 150)
cat("✓ Saved: pixel_approach/charts/agreement_venn.pdf\n")
