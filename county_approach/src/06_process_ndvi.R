library(terra)
library(tidyverse)
library(lubridate)

cat("Processing NDVI data from AppEEARS NetCDF...\n")
write_cluster_files <- FALSE
source("county_approach/src/cluster_polygons.R")

nc_file <- "data_common/raw/ndvi/MOD13A1.061_500m_aid0001.nc"
if (!file.exists(nc_file)) {
  stop("NDVI NetCDF not found. Place file at ", nc_file)
}

# Load NDVI raster
cat("Loading NDVI layer from NetCDF...\n")
r <- rast(nc_file, subds = "_500m_16_days_NDVI")
dates_16day <- as.Date(time(r))
cat(sprintf(
  "✓ Loaded %d time steps (%s to %s)\n",
  length(dates_16day), min(dates_16day), max(dates_16day)
))

# Check scale factor — AppEEARS should apply it, but verify
sample_vals <- as.numeric(values(r[[1]], na.rm = TRUE))
raw_max <- max(abs(sample_vals), na.rm = TRUE)
if (raw_max > 2) {
  cat(sprintf(
    "  Scale factor not applied (max=%.0f) — dividing by 10000\n", raw_max
  ))
  r <- r / 10000
} else {
  cat(sprintf("  Scale factor already applied (range ~%.2f to %.2f)\n",
              min(sample_vals, na.rm = TRUE), max(sample_vals, na.rm = TRUE)))
}

# Extract mean NDVI for each cluster polygon — process in batches for progress
clusters_v <- vect(st_transform(clusters_sf, crs(r)))
n_layers   <- nlyr(r)
batch_size <- 50
batches    <- seq(1, n_layers, by = batch_size)
results    <- list()

cat(sprintf(
  "Extracting mean NDVI per cluster polygon (%d time steps)...\n", n_layers
))

for (i in seq_along(batches)) {
  idx_start <- batches[i]
  idx_end   <- min(idx_start + batch_size - 1, n_layers)
  v <- terra::extract(r[[idx_start:idx_end]], clusters_v,
                      fun = mean, na.rm = TRUE)
  results[[i]] <- v[, -1, drop = FALSE]
  cat(sprintf("  [%d / %d] layers processed\n", idx_end, n_layers))
}

vals_mat <- do.call(cbind, results)

# Build tidy dataframe
ndvi_16day <- bind_rows(
  data.frame(
    date    = dates_16day,
    cluster = "Sequoia",
    ndvi    = round(as.numeric(vals_mat[1, ]), 4)
  ),
  data.frame(
    date    = dates_16day,
    cluster = "Plumas",
    ndvi    = round(as.numeric(vals_mat[2, ]), 4)
  )
) |>
  filter(!is.na(ndvi), ndvi > -0.3) |>
  arrange(cluster, date)

cat(sprintf(
  "✓ Extracted %d 16-day records (%d per cluster)\n",
  nrow(ndvi_16day), as.integer(nrow(ndvi_16day) / 2)
))

# Aggregate to monthly
ndvi_monthly <- ndvi_16day |>
  mutate(date = floor_date(date, "month")) |>
  group_by(date, cluster) |>
  summarise(ndvi = round(mean(ndvi, na.rm = TRUE), 4), .groups = "drop") |>
  arrange(cluster, date)

write_csv(ndvi_16day,   "county_approach/data/processed/Sierra_Nevada_NDVI_16day_2000-2026.csv")
write_csv(ndvi_monthly, "county_approach/data/processed/Sierra_Nevada_NDVI_monthly_2000-2026.csv")

cat("\n✓ 16-day NDVI saved:")
cat("  county_approach/data/processed/Sierra_Nevada_NDVI_16day_2000-2026.csv\n")
cat("✓ Monthly NDVI saved:")
cat("  county_approach/data/processed/Sierra_Nevada_NDVI_monthly_2000-2026.csv\n")
cat(sprintf(
  "  Total monthly records: %d (%d per cluster)\n",
  nrow(ndvi_monthly), as.integer(nrow(ndvi_monthly) / 2)
))
cat(sprintf(
  "  Monthly NDVI range: %.4f to %.4f\n",
  min(ndvi_monthly$ndvi), max(ndvi_monthly$ndvi)
))
cat("  Spatial method: county polygons ∩ Sierra Nevada boundary\n")
