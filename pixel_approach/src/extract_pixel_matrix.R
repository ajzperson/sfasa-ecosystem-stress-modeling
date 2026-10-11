library(terra)
library(sf)
library(dplyr)
library(readr)
library(lubridate)

# Shared by 01_extract_pixel_ndvi.R/ndwi.R: masks an index raster to the
# boundary, averages it to monthly, and saves the pixel x month matrix.
#
# r                    - SpatRaster, one layer per date
# dates                - Date per layer of r
# sn_boundary          - sf polygon to crop/mask to
# metadata_path        - where to write pixel_metadata.csv (NULL to skip)
# verify_metadata_path - existing pixel_metadata.csv this grid must match
# sample_only          - preview the first n_sample pixels and stop
extract_pixel_matrix <- function(r, dates, sn_boundary,
                                  matrix_path, metadata_path = NULL,
                                  verify_metadata_path = NULL,
                                  matrix_label = "index",
                                  sample_only = FALSE, n_sample = 10) {

  sn_vect <- sn_boundary |> st_transform(crs(r)) |> st_union() |> vect()

  cat(sprintf("Masking %s to Sierra Nevada boundary...\n", matrix_label))
  r <- crop(r, sn_vect)
  r <- mask(r, sn_vect)
  cat("✓ Cropped and masked\n")

  cat("Aggregating to monthly means...\n")
  month_labels <- as.character(floor_date(dates, "month"))
  r_monthly    <- tapp(r, index = month_labels, fun = mean, na.rm = TRUE)
  cat(sprintf("✓ %d monthly layers produced\n", nlyr(r_monthly)))

  # Name layers by month; fails if any month is missing
  target_months <- seq(floor_date(min(dates), "month"),
                        floor_date(max(dates), "month"), by = "month")
  stopifnot(nlyr(r_monthly) == length(target_months))
  names(r_monthly) <- as.character(target_months)

  # Pixels inside the boundary, counted on a constant raster so cloud NAs don't drop any
  sn_mask     <- mask(init(r_monthly[[1]], fun = 1), sn_vect)
  valid_cells <- which(!is.na(values(sn_mask)[, 1]))
  n_pixels    <- length(valid_cells)
  cat(sprintf("✓ %d valid pixels within Sierra Nevada boundary\n", n_pixels))

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

  if (!is.null(verify_metadata_path)) {
    ref <- read_csv(verify_metadata_path, show_col_types = FALSE)
    stopifnot(
      "pixel count doesn't match reference metadata" =
        nrow(ref) == nrow(pixel_metadata),
      "lon/lat don't match reference metadata" =
        all(ref$lon == pixel_metadata$lon, ref$lat == pixel_metadata$lat)
    )
    cat(sprintf("✓ Verified %d pixels match %s (shared grid)\n",
                nrow(ref), verify_metadata_path))
  }

  if (sample_only) {
    cat(sprintf("\n--- SAMPLE MODE: first %d pixels ---\n\n", n_sample))

    sample_idx   <- seq_len(n_sample)
    sample_cells <- valid_cells[sample_idx]
    n_show       <- min(6, nlyr(r_monthly))
    sample_mat   <- values(r_monthly[[seq_len(n_show)]])[sample_cells, , drop = FALSE]
    colnames(sample_mat) <- as.character(target_months)[seq_len(ncol(sample_mat))]

    sample_out <- bind_cols(
      pixel_metadata[sample_idx, ],
      as.data.frame(round(sample_mat, 4))
    )

    cat("pixel_metadata (first rows):\n")
    print(pixel_metadata[sample_idx, ])

    cat(sprintf("\n%s matrix sample (first %d months shown):\n", matrix_label, n_show))
    print(sample_out)

    cat(sprintf(
      "\nExpected full matrix dimensions: %d rows x %d cols\n",
      n_pixels, length(target_months)
    ))
    cat("Set SAMPLE_ONLY <- FALSE to run the full extraction.\n")
    stop("Sample complete — review output above, then set SAMPLE_ONLY <- FALSE.",
         call. = FALSE)
  }

  cat(sprintf("\nExtracting full %s pixel matrix (this will take a few minutes)...\n",
              matrix_label))

  out_matrix           <- values(r_monthly)[valid_cells, , drop = FALSE]
  colnames(out_matrix) <- as.character(target_months)
  rownames(out_matrix) <- NULL

  cat(sprintf("✓ Matrix dimensions: %d x %d\n", nrow(out_matrix), ncol(out_matrix)))
  cat(sprintf("  %s range: %.4f to %.4f (NA count: %d)\n", matrix_label,
              min(out_matrix, na.rm = TRUE),
              max(out_matrix, na.rm = TRUE),
              sum(is.na(out_matrix))))

  dir.create(dirname(matrix_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(out_matrix, matrix_path)
  cat(sprintf("✓ Saved: %s\n", matrix_path))

  if (!is.null(metadata_path)) {
    dir.create(dirname(metadata_path), showWarnings = FALSE, recursive = TRUE)
    write_csv(pixel_metadata, metadata_path)
    cat(sprintf("✓ Saved: %s\n", metadata_path))
    cat(sprintf("  pixel_metadata: %d rows x %d cols\n",
                nrow(pixel_metadata), ncol(pixel_metadata)))
  }
  cat(sprintf("  %s_matrix:    %d rows x %d cols\n",
              matrix_label, nrow(out_matrix), ncol(out_matrix)))

  invisible(list(pixel_metadata = pixel_metadata, matrix = out_matrix))
}
