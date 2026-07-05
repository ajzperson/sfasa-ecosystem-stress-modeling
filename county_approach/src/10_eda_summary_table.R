library(tidyverse)

ts_data  <- read_csv("county_approach/data/processed/Sierra_Nevada_TimeSeries_2000-2025.csv")
fire_data <- read_csv(
  "county_approach/data/processed/Sierra_Nevada_Fire_Incidents_2013-2016.csv"
)

# aggregate climate data by year and cluster
climate_annual <- ts_data |>
  filter(date <= as.Date("2016-12-31")) |>
  mutate(year = as.integer(format(date, "%Y"))) |>
  group_by(year, cluster) |>
  summarise(
    mean_ndvi  = round(mean(ndvi, na.rm = TRUE), 4),
    mean_vpd   = round(mean(vpd,  na.rm = TRUE), 3),
    total_ppt  = round(sum(ppt,   na.rm = TRUE), 1),
    mean_tmean = round(mean(tmean, na.rm = TRUE), 2),
    .groups = "drop"
  )

# aggregate fire data by year and cluster
fire_annual <- fire_data |>
  group_by(year = ArchiveYear, cluster) |>
  summarise(
    fire_count  = n(),
    total_acres = round(sum(AcresBurned, na.rm = TRUE), 0),
    .groups = "drop"
  )

# join
summary_table <- climate_annual |>
  left_join(fire_annual, by = c("year", "cluster")) |>
  arrange(year, cluster)

print(summary_table, n = Inf)

write_csv(summary_table, "county_approach/data/processed/Sierra_Nevada_AnnualSummary.csv")

cat("Done. Saved to county_approach/data/processed/Sierra_Nevada_AnnualSummary.csv\n")
