library(tidyverse)
library(sf)

# Three-circle Venn of NDVI, NDWI and combined most-stressed pixels, with
# regions G1-G4; prints the group rates.
#
# Output: pixel_approach/charts/overlap_venn.pdf

d <- read_csv("pixel_approach/data/processed/combined_vs_agreement.csv", show_col_types = FALSE)

rate <- function(subset, col) {
  n <- nrow(subset)
  if (n == 0) return(c(n = 0, pct = 0))
  c(n = n, pct = round(100 * mean(subset[[col]] == TRUE, na.rm = TRUE), 1))
}

ndvi3 <- filter(d, ndvi_cluster == 3)
ndwi3 <- filter(d, ndwi_cluster == 3)
O1    <- filter(d, in_agreement_both)          # NDVI_3 & NDWI_3
O2    <- filter(d, in_combined_top)            # NDVI+NDWI_4
both  <- filter(d, in_agreement_both & in_combined_top)

n_ndvi3 <- nrow(ndvi3)
n_ndwi3 <- nrow(ndwi3)
o1_n    <- rate(O1, "any_mortality");   o1_fire  <- rate(O1, "burned_validation")
o2_n    <- rate(O2, "any_mortality");   o2_fire  <- rate(O2, "burned_validation")
both_n  <- rate(both, "any_mortality"); both_fire <- rate(both, "burned_validation")

cat(sprintf("NDVI-3 n=%d | NDWI-3 n=%d\n", n_ndvi3, n_ndwi3))
cat(sprintf("O1 n=%.0f beetle=%.1f%% fire=%.1f%%\n", o1_n["n"], o1_n["pct"], o1_fire["pct"]))
cat(sprintf("O2 n=%.0f beetle=%.1f%% fire=%.1f%%\n", o2_n["n"], o2_n["pct"], o2_fire["pct"]))
cat(sprintf("Both n=%.0f beetle=%.1f%% fire=%.1f%%\n", both_n["n"], both_n["pct"], both_fire["pct"]))

# The slide legend's numbers come from this printout

circle_pts <- function(cx, cy, r, n = 200) {
  theta <- seq(0, 2 * pi, length.out = n)
  tibble(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

r <- 1.35
dist <- 0.78
deg2rad <- function(deg) deg * pi / 180
centers <- list(
  ndvi     = c(dist * cos(deg2rad(90)),  dist * sin(deg2rad(90))),
  ndwi     = c(dist * cos(deg2rad(210)), dist * sin(deg2rad(210))),
  combined = c(dist * cos(deg2rad(330)), dist * sin(deg2rad(330)))
)

circ_ndvi     <- circle_pts(centers$ndvi[1],     centers$ndvi[2],     r)
circ_ndwi     <- circle_pts(centers$ndwi[1],     centers$ndwi[2],     r)
circ_combined <- circle_pts(centers$combined[1], centers$combined[2], r)

# G4 = all three, G3 = NDVI ∩ NDWI only, G2 = rest of Combined, G1 = everything else
poly_A <- st_polygon(list(cbind(c(circ_ndvi$x, circ_ndvi$x[1]),     c(circ_ndvi$y, circ_ndvi$y[1]))))
poly_B <- st_polygon(list(cbind(c(circ_ndwi$x, circ_ndwi$x[1]),     c(circ_ndwi$y, circ_ndwi$y[1]))))
poly_C <- st_polygon(list(cbind(c(circ_combined$x, circ_combined$x[1]), c(circ_combined$y, circ_combined$y[1]))))

AB <- st_intersection(poly_A, poly_B)
g4_geom <- st_intersection(AB, poly_C)                # A ∩ B ∩ C
g3_geom <- st_difference(AB, poly_C)                   # A ∩ B, excl. C
g2_geom <- st_difference(poly_C, g4_geom)              # C, excl. the G4 center

to_df <- function(geom, group) {
  # Keep disjoint pieces (L2) as separate polygons
  m <- st_coordinates(geom)
  tibble(x = m[, "X"], y = m[, "Y"], piece = if ("L2" %in% colnames(m)) m[, "L2"] else 1, group = group)
}
regions <- bind_rows(to_df(g2_geom, "G2"), to_df(g3_geom, "G3"), to_df(g4_geom, "G4"))

shade_alpha <- c(G4 = 0.42, G3 = 0.22, G2 = 0.12)
outer_center <- c(mean(c(centers$ndvi[1], centers$ndwi[1], centers$combined[1])),
                   mean(c(centers$ndvi[2], centers$ndwi[2], centers$combined[2])))
outer_r <- max(sapply(centers, function(c) sqrt(sum((c - outer_center)^2)))) + r + 0.45
circ_outer <- circle_pts(outer_center[1], outer_center[2], outer_r)

p <- ggplot() +
  geom_polygon(data = circ_outer, aes(x, y), fill = NA, color = "grey40", linewidth = 0.5, linetype = "dashed") +
  geom_polygon(data = circ_ndvi,     aes(x, y), fill = "#E2655E", alpha = 0.95, color = NA) +
  geom_polygon(data = circ_ndwi,     aes(x, y), fill = "#EC9B97", alpha = 0.95, color = NA) +
  geom_polygon(data = circ_combined, aes(x, y), fill = "#D73027", alpha = 0.95, color = NA) +
  geom_polygon(data = regions, aes(x, y, group = interaction(group, piece), alpha = group),
               fill = "black", color = "black", linewidth = 0.35) +
  scale_alpha_manual(values = shade_alpha, guide = "none") +
  annotate("text", x = 0,     y = -0.05, label = "G4", size = 5, fontface = "bold", color = "white") +
  annotate("text", x = -0.65, y = 0.5,   label = "G3", size = 5, fontface = "bold", color = "#1B2A1F") +
  annotate("text", x = 1.15,  y = -0.75, label = "G2", size = 5, fontface = "bold", color = "#1B2A1F") +
  annotate("text", x = 0.65,  y = 0.5,   label = "G2", size = 5, fontface = "bold", color = "#1B2A1F") +
  annotate("text", x = 0,     y = -1.0,  label = "G2", size = 5, fontface = "bold", color = "#1B2A1F") +
  annotate("text", x = 0,     y = 1.55,  label = "G1", size = 5, fontface = "bold", color = "#1B2A1F") +
  annotate("text", x = -1.3,  y = -0.9,  label = "G1", size = 5, fontface = "bold", color = "#1B2A1F") +
  annotate("text", x = 0,     y = -2.15, label = "G1", size = 5, fontface = "bold", color = "#1B2A1F") +
  coord_fixed(clip = "off") +
  xlim(-outer_r - 0.2, outer_r + 0.2) + ylim(-outer_r - 0.2, outer_r + 0.2) +
  theme_void()

ggsave("pixel_approach/charts/overlap_venn.pdf", plot = p, width = 4.2, height = 4.2, dpi = 150)
cat("✓ Saved: pixel_approach/charts/overlap_venn.pdf\n")
