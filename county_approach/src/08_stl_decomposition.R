library(tidyverse)

ts_data <- read_csv("county_approach/data/processed/Sierra_Nevada_TimeSeries_2000-2014.csv",
                    show_col_types = FALSE)

run_stl <- function(cluster_name) {
  dates <- ts_data |>
    filter(cluster == cluster_name) |>
    arrange(date) |>
    pull(date)

  x <- ts_data |>
    filter(cluster == cluster_name) |>
    arrange(date) |>
    pull(ndvi)

  # monthly data starting Feb 2000
  x_ts <- ts(x, start = c(2000, 2), frequency = 12)

  decomp <- stl(x_ts, s.window = "periodic", robust = TRUE)

  # build tidy data frame of all 4 components
  stl_df <- tibble(
    date      = dates,
    data      = as.numeric(x_ts),
    seasonal  = as.numeric(decomp$time.series[, "seasonal"]),
    trend     = as.numeric(decomp$time.series[, "trend"]),
    remainder = as.numeric(decomp$time.series[, "remainder"])
  ) |>
    pivot_longer(-date, names_to = "component", values_to = "value") |>
    mutate(component = factor(
      component,
      levels = c("data", "seasonal", "trend", "remainder")
    ))

  p <- ggplot(stl_df, aes(x = date, y = value)) +
    geom_line(
      data = ~ filter(.x, component != "remainder"),
      linewidth = 0.4
    ) +
    geom_col(
      data = ~ filter(.x, component == "remainder"),
      width = 20, fill = "black"
    ) +
    facet_wrap(~component, ncol = 1, scales = "free_y") +
    scale_x_date(
      breaks = seq(as.Date("2000-01-01"), as.Date("2015-01-01"), by = "1 year"),
      date_labels = "%Y"
    ) +
    labs(
      title = paste("STL Decomposition —", cluster_name, "NDVI (2000–2014)"),
      x = "Year", y = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
      strip.text     = element_text(face = "bold"),
      plot.title     = element_text(hjust = 0.5)
    )

  ggsave(
    paste0("county_approach/charts/stl_", tolower(cluster_name), ".pdf"),
    plot = p, width = 10, height = 8
  )

  list(
    residual    = round(as.numeric(decomp$time.series[, "remainder"]), 4),
    deseasoned  = as.numeric(x_ts) - as.numeric(decomp$time.series[, "seasonal"]),
    trend       = as.numeric(decomp$time.series[, "trend"]),
    dates       = dates
  )
}

sequoia_out <- run_stl("Sequoia")
plumas_out  <- run_stl("Plumas")

sequoia_residual <- sequoia_out$residual
plumas_residual  <- plumas_out$residual

# overlay chart: STL trend for both clusters
trend_df <- bind_rows(
  tibble(date = sequoia_out$dates, cluster = "Sequoia", trend = sequoia_out$trend),
  tibble(date = plumas_out$dates,  cluster = "Plumas",  trend = plumas_out$trend)
)

p_trend <- ggplot(trend_df, aes(x = date, y = trend, color = cluster)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = as.Date("2012-01-01"), linetype = "dashed", color = "grey40") +
  annotate("text", x = as.Date("2012-03-01"), y = Inf, label = "2012 onset",
           hjust = 0, vjust = 1.5, size = 3, color = "grey40") +
  annotate("text", x = as.Date("2013-01-01"), y = 0.450,
           label = "Sequoia: steady decline,\nno recovery by 2014",
           hjust = 0, size = 3, color = "#D55E00") +
  annotate("text", x = as.Date("2013-01-01"), y = 0.572,
           label = "Plumas: stable,\nno collapse",
           hjust = 0, size = 3, color = "#0072B2") +
  coord_cartesian(ylim = c(0.420, 0.590)) +
  scale_x_date(
    breaks = seq(as.Date("2000-01-01"), as.Date("2015-01-01"), by = "1 year"),
    date_labels = "%Y"
  ) +
  scale_color_manual(values = c("Sequoia" = "#D55E00", "Plumas" = "#0072B2")) +
  labs(
    title    = "STL Trend Component — Sequoia vs Plumas NDVI (2000–2014)",
    subtitle = "Same drought — Sequoia collapsed, Plumas held steady",
    x = "Year", y = "NDVI Trend", color = "Cluster"
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    plot.title       = element_text(hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5),
    legend.position  = "top"
  )

ggsave("county_approach/charts/stl_trend_comparison.pdf", plot = p_trend, width = 10, height = 5)
cat("✓ Saved: county_approach/charts/stl_trend_comparison.pdf\n")

# add residuals back to the 2000-2014 dataset
dates <- ts_data |> filter(cluster == "Sequoia") |> arrange(date) |> pull(date)

residuals_df <- bind_rows(
  tibble(date = dates, cluster = "Sequoia",
         ndvi_residual  = sequoia_residual,
         ndvi_deseasoned = round(sequoia_out$deseasoned, 4)),
  tibble(date = dates, cluster = "Plumas",
         ndvi_residual  = plumas_residual,
         ndvi_deseasoned = round(plumas_out$deseasoned, 4))
)

ts_data <- ts_data |>
  select(-any_of(c("ndvi_residual", "ndvi_residual.x", "ndvi_residual.y",
                   "ndvi_deseasoned", "ndvi_deseasoned.x", "ndvi_deseasoned.y"))) |>
  left_join(residuals_df, by = c("date", "cluster"))

write_csv(ts_data, "county_approach/data/processed/Sierra_Nevada_TimeSeries_2000-2014.csv")

cat("Done. Charts saved to county_approach/charts/. ndvi_residual written to Sierra_Nevada_TimeSeries_2000-2014.csv\n")
