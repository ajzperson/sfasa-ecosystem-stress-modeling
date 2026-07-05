library(tidyverse)

ts_data <- read_csv("county_approach/data/processed/Sierra_Nevada_TimeSeries_2000-2014.csv",
                    show_col_types = FALSE) |>
  filter(!month(date) %in% c(12, 1, 2))  # exclude Dec-Feb: snow masking artifact

WIN <- 54  # rolling window size (~40% of N=135: 9 months/year × 15 years)

# rolling stat functions — each takes a numeric vector of length WIN
#
# roll_r1 needs the window's dates too: filtering out Dec-Feb makes the
# series a concatenation of 9-month blocks, so array-adjacent points are not
# always calendar-adjacent (Nov of year N sits next to Mar of year N+1, a
# 4-month jump). Plain acf() treats array adjacency as time adjacency and
# would silently score that jump as a lag-1 pair, 5 times per 54-point
# window. Compute the lag-1 product only over pairs that are truly one
# calendar month apart; keep the same full-window variance denominator as
# the standard biased ACF estimator (matches earlywarnings convention).
roll_r1 <- function(x, d) {
  xbar <- mean(x)
  midx <- year(d) * 12 + month(d)
  keep <- diff(midx) == 1  # TRUE only for genuine 1-month-apart pairs

  num <- sum((head(x, -1)[keep] - xbar) * (tail(x, -1)[keep] - xbar)) / length(x)
  den <- mean((x - xbar)^2)
  num / den
}

roll_var <- function(x) {
  var(x)
}

roll_skew <- function(x) {
  mu  <- mean(x)
  num <- mean((x - mu)^3)
  den <- (mean((x - mu)^2))^(3/2)
  -(num / den)  # sign flipped by convention: rising gamma = more left-skewed residuals
}

# Climate-correct ndvi_deseasoned: fit regression on pre-drought baseline
# (2000-2011), apply coefficients to full 2000-2014 series. Residuals isolate
# biotic resilience signal from shared VPD/ppt/tmean forcing.
climate_correct <- function(df) {
  train <- df |> filter(date < as.Date("2012-01-01"))
  fit   <- lm(ndvi_deseasoned ~ vpd + ppt, data = train)
  df |> mutate(clim_resid = ndvi_deseasoned - predict(fit, newdata = df))
}

ts_data <- ts_data |>
  group_by(cluster) |>
  group_modify(~ climate_correct(.x)) |>
  ungroup()

compute_csd <- function(cluster_name) {
  df <- ts_data |>
    filter(cluster == cluster_name) |>
    arrange(date)

  dates <- df$date
  n     <- length(dates)
  y     <- df$clim_resid  # ndvi_deseasoned residual after removing VPD + ppt

  # preallocate — first WIN-1 positions are NA (not enough data yet)
  r1   <- rep(NA_real_, n)
  vari <- rep(NA_real_, n)
  skew <- rep(NA_real_, n)

  for (i in WIN:n) {
    window      <- y[(i - WIN + 1):i]
    date_window <- dates[(i - WIN + 1):i]
    r1[i]       <- roll_r1(window, date_window)
    vari[i]     <- roll_var(window)
    skew[i]     <- roll_skew(window)
  }

  tibble(
    date       = dates,
    cluster    = cluster_name,
    clim_resid = round(y,    4),
    r1         = round(r1,   4),
    var        = round(vari, 6),
    skew       = round(skew, 4)
  )
}

csd <- bind_rows(
  compute_csd("Sequoia"),
  compute_csd("Plumas")
)

# Kendall's tau (2000–2014)
kendall_results <- csd |>
  filter(!is.na(r1)) |>
  group_by(cluster) |>
  summarise(
    tau_r1   = round(cor.test(as.numeric(date), r1,   method = "kendall")$estimate, 4),
    p_r1     = round(cor.test(as.numeric(date), r1,   method = "kendall")$p.value,  4),
    tau_var  = round(cor.test(as.numeric(date), var,  method = "kendall")$estimate, 4),
    p_var    = round(cor.test(as.numeric(date), var,  method = "kendall")$p.value,  4),
    tau_skew = round(cor.test(as.numeric(date), skew, method = "kendall")$estimate, 4),
    p_skew   = round(cor.test(as.numeric(date), skew, method = "kendall")$p.value,  4),
    .groups  = "drop"
  )

cat("\nKendall's tau results:\n")
print(kendall_results)
write_csv(kendall_results, "county_approach/data/processed/Sierra_Nevada_KendallTau.csv")
cat("✓ Saved: county_approach/data/processed/Sierra_Nevada_KendallTau.csv\n")

# plot
csd_long <- csd |>
  filter(!is.na(r1)) |>
  pivot_longer(cols = c(r1, var, skew),
               names_to = "stat", values_to = "value") |>
  mutate(stat = factor(stat,
    levels = c("r1", "var", "skew"),
    labels = c("Lag-1 Autocorrelation (r₁)",
               "Variance (σ²)",
               "Skewness (γ)")))

p <- ggplot(csd_long, aes(x = date, y = value, color = cluster)) +
  geom_line(linewidth = 0.6) +
  geom_vline(xintercept = as.Date("2012-01-01"),
             linetype = "dashed", color = "grey40") +
  annotate("text", x = as.Date("2012-03-01"), y = Inf,
           label = "2012 onset", hjust = 0, vjust = 1.5,
           size = 3, color = "grey40") +
  facet_wrap(~stat, ncol = 1, scales = "free_y") +
  scale_x_date(
    breaks = seq(as.Date("2007-01-01"), as.Date("2015-01-01"), by = "1 year"),
    date_labels = "%Y"
  ) +
  scale_color_manual(values = c("Sequoia" = "#D55E00", "Plumas" = "#0072B2")) +
  labs(
    title    = "Critical Slowing Down — Rolling Statistics (window = 54 months, Mar–Nov only)",
    subtitle = "Climate-corrected (VPD + ppt regressed out, baseline 2000–2011); window = 54 months, Mar–Nov",
    x = "Year", y = NULL, color = "Cluster"
  ) +
  theme_bw() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    strip.text      = element_text(face = "bold"),
    plot.title      = element_text(hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5),
    legend.position = "top"
  )

ggsave("county_approach/charts/csd_rolling_stats.pdf", plot = p, width = 10, height = 10)
cat("✓ Saved: county_approach/charts/csd_rolling_stats.pdf\n")

write_csv(csd, "county_approach/data/processed/Sierra_Nevada_CSD_2000-2014.csv")
cat("✓ Saved: county_approach/data/processed/Sierra_Nevada_CSD_2000-2014.csv\n")
