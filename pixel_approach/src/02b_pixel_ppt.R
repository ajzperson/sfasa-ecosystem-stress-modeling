library(tidyverse)
library(terra)

# Interpolates monthly PRISM precipitation (Feb 2000-Dec 2016) at each pixel centroid.
#
# Output: pixel_approach/data/processed/ppt_matrix.rds

target_dates <- seq(as.Date("2000-02-01"), as.Date("2016-12-01"), by = "month")
target_ym    <- format(target_dates, "%Y%m")

cat("Loading pixel metadata...\n")
pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv",
                       show_col_types = FALSE)
pts <- vect(pixel_meta, geom = c("lon", "lat"), crs = "EPSG:4326")
cat(sprintf("✓ %d pixel centroids\n", nrow(pts)))

# Unzip prismz.zip if .bil files are not already present
ppt_check <- list.files("data_common/raw/prism/prismz",
                        pattern = "prism_ppt_us_25m_[0-9]{6}\\.bil$",
                        recursive = TRUE)
if (length(ppt_check) == 0) {
  zip_path <- "data_common/raw/prism/prismz.zip"
  if (!file.exists(zip_path))
    stop("Neither data_common/raw/prism/prismz/ .bil files nor prismz.zip found.")
  cat("Unzipping prismz.zip — this may take several minutes...\n")
  unzip(zip_path, exdir = "data_common/raw/prism")
  cat("✓ Unzipped\n")
}

ppt_files <- list.files("data_common/raw/prism/prismz",
                        pattern = "prism_ppt_us_25m_[0-9]{6}\\.bil$",
                        recursive = TRUE, full.names = TRUE)

if (length(ppt_files) == 0)
  stop("No PRISM ppt .bil files found after unzip — check prismz.zip contents.")

ppt_ym <- sub(".*_([0-9]{6})\\.bil$", "\\1", basename(ppt_files))

keep_p    <- ppt_ym %in% target_ym
ppt_files <- ppt_files[keep_p]
ppt_ym    <- ppt_ym[keep_p]
ord       <- order(ppt_ym)
ppt_files <- ppt_files[ord]
ppt_ym    <- ppt_ym[ord]
ppt_dates <- as.Date(paste0(ppt_ym, "01"), "%Y%m%d")

cat(sprintf("✓ PPT files: %d (%s to %s)\n", length(ppt_files),
            format(ppt_dates[1]), format(ppt_dates[length(ppt_dates)])))

cat("Stacking PRISM rasters...\n")
ppt_rast <- rast(ppt_files)

# Project pixel points to PRISM CRS (NAD83), crop before extraction
pts_ppt  <- project(pts, crs(ppt_rast))
ppt_rast <- crop(ppt_rast, ext(pts_ppt) + 0.5)
cat("✓ Cropped to Sierra Nevada extent\n")

cat("Extracting PPT at pixel centroids (bilinear interpolation)...\n")
ppt_vals <- terra::extract(ppt_rast, pts_ppt, method = "bilinear", ID = FALSE)
colnames(ppt_vals) <- format(ppt_dates, "%Y-%m-%d")

ppt_mat <- as.matrix(ppt_vals)
cat(sprintf("✓ PPT range: %.1f to %.1f mm | NA count: %d\n",
            min(ppt_mat, na.rm = TRUE), max(ppt_mat, na.rm = TRUE),
            sum(is.na(ppt_mat))))
cat(sprintf("✓ Dimensions: %d pixels × %d months\n", nrow(ppt_mat), ncol(ppt_mat)))

saveRDS(ppt_mat, "pixel_approach/data/processed/ppt_matrix.rds")
cat("✓ Saved: pixel_approach/data/processed/ppt_matrix.rds\n")
