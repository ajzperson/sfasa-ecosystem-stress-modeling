library(tidyverse)

ts_data <- read_csv("county_approach/data/processed/Sierra_Nevada_TimeSeries_2000-2025.csv")

collapse_start <- as.Date("2012-01-01")
collapse_end   <- as.Date("2016-12-31")

ts_long <- ts_data |>
  filter(date <= as.Date("2016-12-31")) |>
  pivot_longer(
    cols      = c(vpd, ppt, tmean, tmin, tmax),
    names_to  = "variable",
    values_to = "value"
  ) |>
  mutate(variable = factor(
    variable,
    levels = c("vpd", "ppt", "tmean", "tmin", "tmax"),
    labels = c("VPD (kPa)", "Precipitation (mm)",
               "Mean Temp (°C)", "Min Temp (°C)", "Max Temp (°C)")
  ))

collapse_band <- data.frame(
  xmin = collapse_start, xmax = collapse_end,
  ymin = -Inf, ymax = Inf
)

p <- ggplot(ts_long, aes(x = date, y = value, color = cluster)) +
  geom_rect(
    data = collapse_band,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "red", alpha = 0.12, inherit.aes = FALSE
  ) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") +
  scale_x_date(
    breaks      = seq(as.Date("2000-01-01"), as.Date("2017-01-01"), by = "1 year"),
    date_labels = "%Y"
  ) +
  scale_color_manual(values = c("Sequoia" = "#D55E00", "Plumas" = "#0072B2")) +
  labs(
    title    = "Sierra Nevada Climate Time Series 2000–2016",
    subtitle = "Same drought, same temperatures, same precipitation — both clusters were exposed equally\nShaded region: 2012–2016 collapse window",
    x = "Year", y = NULL, color = "Cluster"
  ) +
  theme_bw() +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 6),
    strip.text    = element_text(face = "bold"),
    plot.title    = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

ggsave(
  "county_approach/charts/climate_timeseries.pdf",
  plot = p, width = 14, height = 18
)

cat("Done. Saved to county_approach/charts/climate_timeseries.pdf\n")
