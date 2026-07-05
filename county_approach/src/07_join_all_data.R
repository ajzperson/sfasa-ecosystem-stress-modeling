library(tidyverse)
library(lubridate)

cat("Joining NDVI, VPD, and PRISM data by cluster...\n")

ndvi  <- read_csv("county_approach/data/processed/Sierra_Nevada_NDVI_monthly_2000-2026.csv",
                  show_col_types = FALSE)
vpd   <- read_csv("county_approach/data/processed/Sierra_Nevada_VPD_2000-2025.csv",
                  show_col_types = FALSE)
prism <- read_csv("county_approach/data/processed/Sierra_Nevada_PRISM_2000-2026.csv",
                  show_col_types = FALSE)
cat(sprintf("  NDVI:  %d records\n", nrow(ndvi)))
cat(sprintf("  VPD:   %d records\n", nrow(vpd)))
cat(sprintf("  PRISM: %d records\n", nrow(prism)))

result <- ndvi |>
  full_join(vpd,   by = c("date", "cluster")) |>
  full_join(prism, by = c("date", "cluster")) |>
  filter(!is.na(ndvi), !is.na(vpd)) |>
  arrange(cluster, date)

cat(sprintf("\n  Missing NDVI:  %d\n", sum(is.na(result$ndvi))))
cat(sprintf("  Missing VPD:   %d\n", sum(is.na(result$vpd))))
cat(sprintf("  Missing PRISM: %d\n", sum(is.na(result$ppt))))

write_csv(result, "county_approach/data/processed/Sierra_Nevada_TimeSeries_2000-2025.csv")
cat("\n✓ Saved: county_approach/data/processed/Sierra_Nevada_TimeSeries_2000-2025.csv\n")

result_2014 <- result |> filter(date < as.Date("2015-01-01"))
write_csv(result_2014, "county_approach/data/processed/Sierra_Nevada_TimeSeries_2000-2014.csv")
cat("✓ Saved: county_approach/data/processed/Sierra_Nevada_TimeSeries_2000-2014.csv\n")

cat(sprintf("\n  Full range:     %d records (%s to %s)\n",
            nrow(result), as.character(min(result$date)), as.character(max(result$date))))
cat(sprintf("  Analysis range: %d records (%s to %s)\n",
            nrow(result_2014), as.character(min(result_2014$date)), as.character(max(result_2014$date))))
cat("  Columns: date, cluster, ndvi, vpd, ppt, tmean, tmin, tmax\n")
