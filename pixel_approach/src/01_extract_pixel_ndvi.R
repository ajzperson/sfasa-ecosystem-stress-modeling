library(terra)
library(sf)
library(tidyverse)
library(lubridate)

# -----------------------------------------------------------------------
# SAMPLE_ONLY: set TRUE to print a few rows for verification, FALSE for
# the full 251,971-pixel run that writes pixel_metadata.csv +
# ndvi_matrix.rds
# -----------------------------------------------------------------------
if (!exists("SAMPLE_ONLY")) SAMPLE_ONLY <- FALSE
if (!exists("N_SAMPLE"))   N_SAMPLE    <- 10

# --- Load and prepare NDVI raster --------------------------------------
nc_file <- "data_common/raw/ndvi/MOD13A1.061_500m_aid0001.nc"
cat("Loading NDVI raster...\n")
r           <- rast(nc_file, subds = "_500m_16_days_NDVI")
dates_16day <- as.Date(time(r))
cat(sprintf("✓ Loaded %d 16-day layers (%s to %s)\n",
            nlyr(r), min(dates_16day), max(dates_16day)))

# Trim to 2000-2016
keep        <- dates_16day <= as.Date("2016-12-31")
r           <- r[[which(keep)]]
dates_16day <- dates_16day[keep]
cat(sprintf("✓ Trimmed to 2000-2016: %d 16-day layers remaining\n", nlyr(r)))

# --- Sierra Nevada boundary mask ---------------------------------------
cat("Masking to Sierra Nevada boundary...\n")
sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE) |>
  st_transform(crs(r)) |>
  st_union() |>
  vect()

r <- crop(r, sn_boundary)
r <- mask(r, sn_boundary)
cat("✓ Cropped and masked\n")

# --- Aggregate 16-day composites to monthly means ----------------------
cat("Aggregating 16-day layers to monthly means (203 target months)...\n")
month_labels <- as.character(floor_date(dates_16day, "month"))
r_monthly    <- tapp(r, index = month_labels, fun = mean, na.rm = TRUE)
cat(sprintf("✓ %d monthly layers produced\n", nlyr(r_monthly)))

# Rename layers to clean date strings (tapp names them by group label).
# Assumes tapp's groups line up 1:1 and in order with target_months — true for
# the current MODIS data (no missing months), but not guaranteed by tapp() itself,
# so assert it rather than silently mislabeling dates if that ever changes.
target_months       <- seq(as.Date("2000-02-01"), as.Date("2016-12-01"), by = "month")
stopifnot(nlyr(r_monthly) == length(target_months))
names(r_monthly)    <- as.character(target_months)

# --- Identify valid pixels ---------------------------------------------
# Valid cells = exactly the cells mask() kept non-NA. Apply mask() to a
# clean constant raster (no cloud NAs) using the same boundary and same
# spatial rule as the data — gives the true polygon pixel count.
sn_mask     <- mask(init(r_monthly[[1]], fun = 1), sn_boundary)
valid_cells <- which(!is.na(values(sn_mask)[, 1]))
n_pixels    <- length(valid_cells)
cat(sprintf("✓ %d valid pixels within Sierra Nevada boundary\n", n_pixels))

# --- Pixel centroids → WGS84 lon/lat -----------------------------------
cat("Computing pixel centroids (WGS84)...\n")
xy     <- xyFromCell(r_monthly, valid_cells)
pts    <- st_as_sf(as.data.frame(xy), coords = c("x", "y"), crs = crs(r_monthly)) |>
  st_transform(4326)
lonlat <- st_coordinates(pts)

pixel_metadata <- data.frame(
  pixel_id = seq_len(n_pixels),
  lon      = round(lonlat[, 1], 6),
  lat      = round(lonlat[, 2], 6)
)
cat("✓ Pixel metadata ready\n")

# --- SAMPLE mode: print a few rows and stop ----------------------------
if (SAMPLE_ONLY) {
  cat(sprintf("\n--- SAMPLE MODE: first %d pixels ---\n\n", N_SAMPLE))

  sample_idx   <- seq_len(N_SAMPLE)
  sample_cells <- valid_cells[sample_idx]
  n_show       <- min(6, nlyr(r_monthly))
  sample_mat   <- values(r_monthly[[seq_len(n_show)]])[sample_cells, ]
  colnames(sample_mat) <- as.character(target_months)[seq_len(ncol(sample_mat))]

  sample_out <- bind_cols(
    pixel_metadata[sample_idx, ],
    as.data.frame(round(sample_mat, 4))
  )

  cat("pixel_metadata (first rows):\n")
  print(pixel_metadata[sample_idx, ])

  cat("\nNDVI matrix sample (first 6 months shown):\n")
  print(sample_out)

  cat(sprintf(
    "\nExpected full matrix dimensions: %d rows x %d cols\n",
    n_pixels, length(target_months)
  ))
  cat("Set SAMPLE_ONLY <- FALSE to run the full extraction.\n")
  stop("Sample complete — review output above, then set SAMPLE_ONLY <- FALSE.",
       call. = FALSE)
}

# --- FULL run: extract matrix and save ---------------------------------
cat("\nExtracting full pixel matrix (this will take a few minutes)...\n")

# values() returns n_cells x n_layers; subset to valid rows only
ndvi_matrix           <- values(r_monthly)[valid_cells, ]
colnames(ndvi_matrix) <- as.character(target_months)
rownames(ndvi_matrix) <- NULL

cat(sprintf("✓ Matrix dimensions: %d x %d\n", nrow(ndvi_matrix), ncol(ndvi_matrix)))
cat(sprintf("  NDVI range: %.4f to %.4f (NA count: %d)\n",
            min(ndvi_matrix, na.rm = TRUE),
            max(ndvi_matrix, na.rm = TRUE),
            sum(is.na(ndvi_matrix))))

dir.create("pixel_approach/data/processed", showWarnings = FALSE, recursive = TRUE)
write_csv(pixel_metadata, "pixel_approach/data/processed/pixel_metadata.csv")
saveRDS(ndvi_matrix,      "pixel_approach/data/processed/ndvi_matrix.rds")

cat("✓ Saved: pixel_approach/data/pixel_metadata.csv\n")
cat("✓ Saved: pixel_approach/data/ndvi_matrix.rds\n")
cat(sprintf("  pixel_metadata: %d rows x %d cols\n",
            nrow(pixel_metadata), ncol(pixel_metadata)))
cat(sprintf("  ndvi_matrix:    %d rows x %d cols\n",
            nrow(ndvi_matrix), ncol(ndvi_matrix)))
