library(tidyverse)
library(sf)
library(patchwork)

# 2015-2016 damage rates and maps for pixels flagged most-stressed by NDVI
# only, NDWI only, or both.
#
# Output: pixel_approach/charts/cluster_agreement_validation.pdf
#         pixel_approach/charts/cluster_agreement_map.pdf

VALIDATION_YEARS <- 2015:2016

sn_boundary <- st_read("data_common/processed/sierra_nevada_boundary.geojson", quiet = TRUE)

ndvi <- read_csv("pixel_approach/data/processed/ndvi_pixel_clusters.csv", show_col_types = FALSE) |>
  select(pixel_id, lon, lat, cluster_ndvi = cluster)
ndwi <- read_csv("pixel_approach/data/processed/ndwi_pixel_clusters.csv", show_col_types = FALSE) |>
  select(pixel_id, cluster_ndwi = cluster)

fire <- read_csv("pixel_approach/data/processed/pixel_fire.csv", show_col_types = FALSE)
year_cols <- sprintf("burned_%d", VALIDATION_YEARS)
fire$burned_validation <- rowSums(fire[year_cols]) > 0
fire <- fire |> select(pixel_id, burned_validation)

# NDVI's copy of the beetle flags; NDWI's is identical on shared pixels
beetle <- read_csv("pixel_approach/data/processed/pixel_beetle_tpa_ndvi.csv", show_col_types = FALSE) |>
  select(pixel_id, any_mortality)

joined <- ndvi |>
  inner_join(ndwi, by = "pixel_id") |>
  left_join(fire, by = "pixel_id") |>
  left_join(beetle, by = "pixel_id") |>
  mutate(group = case_when(
    cluster_ndvi == 3 & cluster_ndwi == 3 ~ "Both",
    cluster_ndvi == 3 & cluster_ndwi != 3 ~ "NDVI-3 only",
    cluster_ndvi != 3 & cluster_ndwi == 3 ~ "NDWI-3 only",
    TRUE                                  ~ "Neither"
  )) |>
  filter(group != "Neither")

group_order <- c("NDWI-3 only", "NDVI-3 only", "Both")

summary_tbl <- joined |>
  group_by(group) |>
  summarise(
    n             = n(),
    Fire          = 100 * mean(burned_validation, na.rm = TRUE),
    Beetle        = 100 * mean(any_mortality, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(group = factor(group, levels = group_order)) |>
  arrange(group)

cat("=== Cluster-3 agreement vs. real-world validation rates ===\n")
print(summary_tbl)

rr_fire   <- summary_tbl$Fire[summary_tbl$group == "Both"] /
  summary_tbl$Fire[summary_tbl$group == "NDWI-3 only"]
rr_beetle <- summary_tbl$Beetle[summary_tbl$group == "Both"] /
  summary_tbl$Beetle[summary_tbl$group == "NDWI-3 only"]
cat(sprintf("  Relative risk (Both vs NDWI-3 only): Fire %.2fx, Beetle %.2fx\n", rr_fire, rr_beetle))

plot_df <- summary_tbl |>
  pivot_longer(cols = c(Fire, Beetle), names_to = "validator", values_to = "pct") |>
  mutate(validator = factor(validator, levels = c("Fire", "Beetle")))

validator_cols <- c(Fire = "#ED6623", Beetle = "#A0522D")

group_n <- setNames(summary_tbl$n, as.character(summary_tbl$group))
group_x_labels <- sprintf("%s\n(n=%s)", group_order, format(group_n[group_order], big.mark = ","))

p <- ggplot(plot_df, aes(x = group, y = pct, fill = validator)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.75), vjust = -0.4, size = 3, fontface = "bold",
            color = "#1B2A1F") +
  scale_fill_manual(values = validator_cols, name = NULL) +
  scale_x_discrete(labels = setNames(group_x_labels, group_order)) +
  labs(
    title    = "Real-World Damage Rates by NDVI/NDWI Cluster-3 (Most-stressed) Agreement",
    subtitle = sprintf("%% of pixels with each real-world outcome (%s), by whether NDVI and/or NDWI flagged the pixel \"most-stressed\"",
                        if (length(VALIDATION_YEARS) == 1) as.character(VALIDATION_YEARS)
                        else sprintf("%d-%d", min(VALIDATION_YEARS), max(VALIDATION_YEARS))),
    x = NULL, y = "% of pixels",
    caption = sprintf("Relative risk (Both vs NDWI-3 only): Fire %.2f×, Beetle %.2f×", rr_fire, rr_beetle)
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
    plot.subtitle   = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"),
    plot.caption    = element_text(hjust = 0.5, size = 13, face = "bold", color = "#1B2A1F",
                                    margin = margin(t = 14)),
    axis.title      = element_text(face = "bold", color = "#1B2A1F"),
    axis.text       = element_text(face = "bold", size = 10, color = "#1B2A1F"),
    legend.text     = element_text(face = "bold", size = 10, color = "#1B2A1F"),
    legend.position = "top"
  )

out_file <- "pixel_approach/charts/cluster_agreement_validation.pdf"
ggsave(out_file, plot = p, width = 8, height = 6, dpi = 150)
cat(sprintf("\n✓ Saved: %s\n", out_file))
group_cols <- c(
  "NDWI-3 only" = "#377EB8",
  "NDVI-3 only" = "#FF7F00",
  "Both"        = "#984EA3"
)

LAT_STEP <- 0.004167
# Tiles drawn 3x size so they stay visible at full-map scale
TILE_SCALE <- 3
map_df <- joined |>
  mutate(
    group = factor(group, levels = group_order),
    tile_w = LAT_STEP / cos(lat * pi / 180) * TILE_SCALE
  )

group_n <- setNames(summary_tbl$n, as.character(summary_tbl$group))

# coord_fixed (not coord_sf) so the degree labels can be bold
boundary_xy <- as.data.frame(st_coordinates(sn_boundary))
map_lon_labels <- function(x) sprintf("%g°%s", abs(x), ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
map_lat_labels <- function(y) sprintf("%g°%s", abs(y), ifelse(y < 0, "S", ifelse(y > 0, "N", "")))

build_group_panel <- function(g) {
  ggplot(filter(map_df, group == g), aes(x = lon, y = lat, width = tile_w, height = LAT_STEP * TILE_SCALE)) +
    geom_tile(fill = group_cols[[g]]) +
    geom_polygon(data = boundary_xy, aes(x = X, y = Y), fill = NA, color = "grey20",
                 linewidth = 0.4, inherit.aes = FALSE) +
    coord_fixed(ratio = 1 / cos(median(map_df$lat) * pi / 180),
                xlim = range(map_df$lon), ylim = range(map_df$lat)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3), labels = map_lon_labels) +
    scale_y_continuous(labels = map_lat_labels) +
    labs(title = sprintf("%s (n=%s)", g, format(group_n[[g]], big.mark = ",")),
         x = NULL, y = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11, color = "#1B2A1F"),
      axis.text  = element_text(size = 6, color = "#1B2A1F", face = "bold")
    )
}

panels <- lapply(group_order, build_group_panel)

p_map <- (panels[[1]] | panels[[2]] | panels[[3]]) +
  plot_annotation(
    title    = "Where NDVI and NDWI Agree / Disagree on Cluster 3 (\"most-stressed\")",
    subtitle = "Each pixel is one real 500m grid cell -- see cluster_agreement_validation.pdf for damage rates",
    theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
                      plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"))
  )

map_out_file <- "pixel_approach/charts/cluster_agreement_map.pdf"
ggsave(map_out_file, plot = p_map, width = 13, height = 7, dpi = 300)
cat(sprintf("✓ Saved: %s\n", map_out_file))