library(terra)
library(tidyverse)
library(lubridate)

zip_file  <- "data_common/raw/prism/prismz.zip"
prism_dir <- "data_common/raw/prism/prismz"

if (!dir.exists(prism_dir)) {
  cat("Unzipping PRISM data (~1 GB, may take a minute)...\n")
  unzip(zip_file, exdir = "data_common/raw/prism")
  unlink("data_common/raw/prism/__MACOSX", recursive = TRUE)
}

cat("Processing PRISM data for Sequoia and Plumas clusters...\n")
write_cluster_files <- FALSE
source("county_approach/src/cluster_polygons.R")

ppt_files   <- sort(list.files(prism_dir, pattern = "prism_ppt_us_25m_",   full.names = TRUE))
tmean_files <- sort(list.files(prism_dir, pattern = "prism_tmean_us_25m_", full.names = TRUE))
tmin_files  <- sort(list.files(prism_dir, pattern = "prism_tmin_us_25m_",  full.names = TRUE))
tmax_files  <- sort(list.files(prism_dir, pattern = "prism_tmax_us_25m_",  full.names = TRUE))

cat(sprintf("Found %d ppt files\n",   length(ppt_files)))
cat(sprintf("Found %d tmean files\n", length(tmean_files)))
cat(sprintf("Found %d tmin files\n",  length(tmin_files)))
cat(sprintf("Found %d tmax files\n",  length(tmax_files)))

# Pre-project cluster polygons to PRISM CRS (done once, reused for all files)
first_bil <- list.files(ppt_files[1], pattern = "\\.bil$", full.names = TRUE)
prism_crs <- crs(rast(first_bil))
clusters_v <- vect(st_transform(clusters_sf, prism_crs))
cat("✓ Cluster polygons projected to PRISM CRS\n\n")

extract_date <- function(foldername) {
  date_str <- str_extract(foldername, "\\d{6}$")
  as.Date(paste0(date_str, "01"), format = "%Y%m%d")
}

process_variable <- function(files, var_name) {
  records <- list()
  for (i in seq_along(files)) {
    folder   <- files[i]
    bil_file <- list.files(folder, pattern = "\\.bil$", full.names = TRUE)
    if (length(bil_file) == 0) next

    date_val <- extract_date(folder)
    tryCatch({
      r    <- rast(bil_file)
      vals <- terra::extract(r, clusters_v, fun = mean, na.rm = TRUE)
      records[[length(records) + 1]] <- data.frame(
        date    = rep(date_val, 2),
        cluster = c("Sequoia", "Plumas"),
        value   = as.numeric(vals[, 2])
      )
      if (i %% 50 == 0) cat(sprintf("  Processed %d/%d\n", i, length(files)))
    }, error = function(e) {
      cat(sprintf("  Error: %s: %s\n", folder, e$message))
    })
  }
  bind_rows(records) |>
    arrange(cluster, date) |>
    rename(!!var_name := value)
}

cat("Extracting precipitation data...\n")
ppt_df <- process_variable(ppt_files, "ppt")
cat(sprintf("✓ %d records\n\n", nrow(ppt_df)))

cat("Extracting mean temperature data...\n")
tmean_df <- process_variable(tmean_files, "tmean")
cat(sprintf("✓ %d records\n\n", nrow(tmean_df)))

cat("Extracting minimum temperature data...\n")
tmin_df <- process_variable(tmin_files, "tmin")
cat(sprintf("✓ %d records\n\n", nrow(tmin_df)))

cat("Extracting maximum temperature data...\n")
tmax_df <- process_variable(tmax_files, "tmax")
cat(sprintf("✓ %d records\n\n", nrow(tmax_df)))

result <- ppt_df |>
  full_join(tmean_df, by = c("date", "cluster")) |>
  full_join(tmin_df,  by = c("date", "cluster")) |>
  full_join(tmax_df,  by = c("date", "cluster")) |>
  arrange(cluster, date) |>
  mutate(across(c(ppt, tmean, tmin, tmax), \(x) round(x, 2)))

write_csv(result, "county_approach/data/processed/Sierra_Nevada_PRISM_2000-2026.csv")

cat("✓ PRISM data saved: county_approach/data/processed/Sierra_Nevada_PRISM_2000-2026.csv\n")
cat(sprintf("  Total records: %d (%d per cluster)\n",
            nrow(result), as.integer(nrow(result) / 2)))
cat(sprintf("  Date range: %s to %s\n", min(result$date), max(result$date)))
cat("  Columns: date, cluster, ppt (mm), tmean, tmin, tmax (°C)\n")
cat("  Spatial method: county polygons intersected with Sierra Nevada boundary\n")
