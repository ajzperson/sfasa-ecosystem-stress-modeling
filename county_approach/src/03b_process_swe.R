library(terra)
library(tidyverse)

cat("Processing TerraClimate SWE for Sequoia and Plumas clusters...\n")
source("county_approach/src/cluster_polygons.R")

raw_file <- "data_common/raw/swe/TerraClimate_swe_2000-2025.nc"
if (!file.exists(raw_file)) {
  stop("Raw SWE file not found. Run 02b_download_swe.R first.")
}
r <- rast(raw_file)
cat(sprintf("✓ Loaded: %d monthly layers\n", nlyr(r)))

clusters_v <- vect(st_transform(clusters_sf, crs(r)))
vals       <- terra::extract(r, clusters_v, fun = mean, na.rm = TRUE)
dates      <- as.Date(terra::time(r))

swe_df <- bind_rows(
  data.frame(date = dates, cluster = "Sequoia",
             swe = round(as.numeric(vals[1, -1]), 1)),
  data.frame(date = dates, cluster = "Plumas",
             swe = round(as.numeric(vals[2, -1]), 1))
) |>
  arrange(cluster, date)

cat(sprintf("\n  Records: %d (%d per cluster)\n",
            nrow(swe_df), as.integer(nrow(swe_df) / 2)))
cat(sprintf("  Date range: %s to %s\n", min(swe_df$date), max(swe_df$date)))
cat(sprintf("  SWE range — Sequoia: %.1f - %.1f mm\n",
            min(swe_df$swe[swe_df$cluster == "Sequoia"]),
            max(swe_df$swe[swe_df$cluster == "Sequoia"])))
cat(sprintf("  SWE range — Plumas:  %.1f - %.1f mm\n",
            min(swe_df$swe[swe_df$cluster == "Plumas"]),
            max(swe_df$swe[swe_df$cluster == "Plumas"])))

write_csv(swe_df, "county_approach/data/processed/Sierra_Nevada_SWE_2000-2025.csv")
cat("\n✓ Saved: county_approach/data/processed/Sierra_Nevada_SWE_2000-2025.csv\n")
cat("  Columns: date, cluster, swe (mm)\n")
