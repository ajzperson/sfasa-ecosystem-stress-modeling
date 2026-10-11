library(tidyverse)
library(sf)
library(patchwork)
library(cowplot)

# Map of the study area and pixel count, with a zoomed inset of 500m pixels.
#
# Output: pixel_approach/charts/study_area_pixels.{pdf,png}

dir.create("pixel_approach/charts", showWarnings = FALSE, recursive = TRUE)

sn_boundary <- st_read("data_common/processed/sierra_nevada_boundary.geojson", quiet = TRUE)
pixel_meta  <- read_csv("pixel_approach/data/processed/pixel_metadata.csv", show_col_types = FALSE)
n_pixels    <- nrow(pixel_meta)
area_km2    <- round(as.numeric(st_area(st_transform(sn_boundary, 3310))) / 1e6, 1)

cat(sprintf("✓ %s pixels across %s km²\n", format(n_pixels, big.mark = ","), format(area_km2, big.mark = ",")))

# Zoom on the median centroid so the inset isn't empty
LAT_STEP <- 0.004167
center_lon <- median(pixel_meta$lon)
center_lat <- median(pixel_meta$lat)
zoom_half_deg <- 0.03
ZOOM_COLOR <- "#D73027"

zoom_pixels <- pixel_meta |>
  filter(lon > center_lon - zoom_half_deg, lon < center_lon + zoom_half_deg,
         lat > center_lat - zoom_half_deg, lat < center_lat + zoom_half_deg) |>
  mutate(tile_w = LAT_STEP / cos(lat * pi / 180))
cat(sprintf("✓ Inset window: %d pixels\n", nrow(zoom_pixels)))

# coord_fixed (not coord_sf) so the degree labels can be bold
boundary_xy <- as.data.frame(st_coordinates(sn_boundary))

lon_labels <- function(x) sprintf("%g°%s", abs(x), ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
lat_labels <- function(y) sprintf("%g°%s", abs(y), ifelse(y < 0, "S", ifelse(y > 0, "N", "")))

p_main <- ggplot() +
  geom_polygon(data = boundary_xy, aes(x = X, y = Y),
               fill = "#DCEEDC", color = "grey20", linewidth = 0.5) +
  annotate("rect",
           xmin = center_lon - zoom_half_deg, xmax = center_lon + zoom_half_deg,
           ymin = center_lat - zoom_half_deg, ymax = center_lat + zoom_half_deg,
           fill = NA, color = ZOOM_COLOR, linewidth = 0.6) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3), labels = lon_labels) +
  scale_y_continuous(labels = lat_labels) +
  coord_fixed(ratio = 1 / cos(center_lat * pi / 180)) +
  labs(
    title    = "Sierra Nevada Ecoregion Study Area",
    subtitle = sprintf("%s km² · %s pixels (\"500m\" MODIS grid cells, actually ~463m)",
                       format(area_km2, big.mark = ","), format(n_pixels, big.mark = ",")),
    x = "Longitude", y = "Latitude"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
    plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", face = "bold"),
    axis.title    = element_text(face = "bold", color = "#1B2A1F"),
    axis.text     = element_text(face = "bold", color = "#1B2A1F")
  )

p_inset <- ggplot(zoom_pixels, aes(x = lon, y = lat, width = tile_w, height = LAT_STEP)) +
  # Thicker gridlines survive slide image compression
  geom_tile(fill = "#4DAF4A", color = "white", linewidth = 0.4) +
  labs(
    title = "Zoomed: individual pixels (\"500m\", actually ~463m)",
    subtitle = "Each tile = one pixel, identified by pixel_id\nand a fixed lon/lat centroid",
    x = "Longitude", y = "Latitude"
  ) +
  # 1/cos(lat) so square ground tiles render square
  coord_fixed(ratio = 1 / cos(center_lat * pi / 180)) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 12, color = "#1B2A1F"),
    plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", face = "bold", size = 11),
    axis.title    = element_text(face = "bold", color = "#1B2A1F"),
    axis.text     = element_text(face = "bold", size = 10, color = "#1B2A1F"),
    panel.border  = element_rect(color = ZOOM_COLOR, linewidth = 1, fill = NA)
  )

p <- (p_main | p_inset) +
  plot_layout(widths = c(1.3, 1))

# Arrow from the locator box to the inset (npc coordinates)
p_final <- ggdraw() +
  draw_plot(p) +
  draw_line(
    x = c(0.305, 0.60), y = c(0.545, 0.545),
    color = ZOOM_COLOR, linewidth = 1,
    arrow = arrow(length = unit(0.12, "inches"), type = "closed")
  )

out_file <- "pixel_approach/charts/study_area_pixels.pdf"
ggsave(out_file, plot = p_final, width = 12, height = 7, dpi = 150)
cat(sprintf("✓ Saved: %s\n", out_file))

# Higher-dpi PNG for the slides
out_file_png <- "pixel_approach/charts/study_area_pixels.png"
ggsave(out_file_png, plot = p_final, width = 12, height = 7, dpi = 500)
cat(sprintf("✓ Saved: %s\n", out_file_png))