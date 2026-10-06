library(tidyverse)
library(sf)
library(patchwork)

source("pixel_approach/src/clustering_helpers.R")

# Per-pixel probability alternative to 06_clustering.R + 08/09_*_validation.R:
# instead of discretizing tau_r1/tau_skew into a K-means tier and reporting
# a mortality/fire RATE per tier, fits a logistic regression directly on the
# continuous CSD features so every pixel gets its own predicted probability.
#
# Feature set is the combined NDVI+NDWI tau_r1/tau_var/tau_skew (6 features)
# -- tau_var is included here even though 06/06b leave it out of clustering
# (variance's direction near a transition isn't theoretically guaranteed,
# unlike r1/skew, so it was dropped there), since a regression isn't
# limited to 2-3 dimensions the way a K-means map/legend is. Same per-index
# baseline-density residualization as clustering (each index regressed
# against its OWN base_level) to remove the same denser-forest confound
# documented in 06_clustering.R. Same pre-2015 MTBS burned-pixel exclusion,
# for the same reason (a fire crash in the feature window isn't CSD).
#
# Fit twice: once on a 70/30 train/test split, to get an honest held-out
# AUC and calibration curve; once on all eligible pixels, to produce the
# actual per-pixel probability used for the map and CSV output (using only
# the train fold would throw away signal for the test-fold pixels' own
# predictions, and there's no reuse of these probabilities as a training
# label anywhere downstream, so refitting on everything after validating
# generalization is safe).
#
# Output:
#   pixel_approach/data/processed/pixel_regression/pixel_regression_predictions.csv
#   pixel_approach/charts/pixel_regression/pixel_regression_prob_map.pdf
#   pixel_approach/charts/pixel_regression/pixel_regression_calibration.pdf

if (!exists("SEED"))       SEED       <- 42
if (!exists("TRAIN_FRAC")) TRAIN_FRAC <- 0.7

dir.create("pixel_approach/charts/pixel_regression", showWarnings = FALSE, recursive = TRUE)
dir.create("pixel_approach/data/processed/pixel_regression", showWarnings = FALSE, recursive = TRUE)

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

fire_flags <- read_csv("pixel_approach/data/processed/pixel_fire.csv", show_col_types = FALSE) |>
  select(pixel_id, burned_baseline_00_14, burned_2015, burned_2016)
fire_flags$burned_validation_15_16 <- fire_flags$burned_2015 | fire_flags$burned_2016

beetle_flags <- read_csv("pixel_approach/data/processed/pixel_beetle.csv", show_col_types = FALSE)
beetle_flags$any_mortality <- pmax(beetle_flags$tpa_2015, beetle_flags$tpa_2016, na.rm = TRUE) > 0

pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv", show_col_types = FALSE) |>
  select(pixel_id)

# ----------------------------------------------------------------------------
# Feature table: combined NDVI+NDWI tau_r1/tau_var/tau_skew, same inner-join
# + baseline-density residualization as 06b_clustering_combined.R.
# ----------------------------------------------------------------------------
cat("Loading NDVI and NDWI tau outputs...\n")
ndvi_tau <- read_csv("pixel_approach/data/processed/ndvi_pixel_tau.csv", show_col_types = FALSE) |>
  filter(!is.na(tau_r1), !is.na(tau_var), !is.na(tau_skew)) |>
  select(pixel_id, lon, lat, tau_r1, tau_var, tau_skew)

ndwi_tau <- read_csv("pixel_approach/data/processed/ndwi_pixel_tau.csv", show_col_types = FALSE) |>
  filter(!is.na(tau_r1), !is.na(tau_var), !is.na(tau_skew)) |>
  select(pixel_id, tau_r1, tau_var, tau_skew)

d <- ndvi_tau |> inner_join(ndwi_tau, by = "pixel_id", suffix = c("_ndvi", "_ndwi"))
cat(sprintf("✓ %d pixels with valid tau_r1/tau_var/tau_skew for both NDVI and NDWI\n", nrow(d)))

n_before <- nrow(d)
d <- d |>
  left_join(fire_flags |> select(pixel_id, burned_baseline_00_14), by = "pixel_id") |>
  filter(!burned_baseline_00_14) |>
  select(-burned_baseline_00_14)
cat(sprintf("✓ Excluded %d pixels burned (MTBS moderate+) during 2000-2014 baseline (%d remain)\n",
            n_before - nrow(d), nrow(d)))

cat("Computing baseline vegetation density (NDVI and NDWI, each against its own index)...\n")
base_level_ndvi <- compute_base_level("pixel_approach/data/processed/ndvi_matrix.rds",
                                      pixel_meta$pixel_id, "base_level_ndvi")
base_level_ndwi <- compute_base_level("pixel_approach/data/processed/ndwi_matrix.rds",
                                      pixel_meta$pixel_id, "base_level_ndwi")

d <- d |>
  left_join(base_level_ndvi, by = "pixel_id") |>
  left_join(base_level_ndwi, by = "pixel_id")
stopifnot(!anyNA(d$base_level_ndvi), !anyNA(d$base_level_ndwi))

# Single-column version of clustering_helpers.R's residualize_tau() -- that
# helper always residualizes an r1/skew pair together, but tau_var here
# needs the same treatment as a third, independent column.
residualize_col <- function(d, col, base_col, out_col) {
  d[[out_col]] <- resid(lm(reformulate(base_col, col), data = d))
  d
}

d <- d |>
  residualize_col("tau_r1_ndvi",   "base_level_ndvi", "resid_tau_r1_ndvi") |>
  residualize_col("tau_var_ndvi",  "base_level_ndvi", "resid_tau_var_ndvi") |>
  residualize_col("tau_skew_ndvi", "base_level_ndvi", "resid_tau_skew_ndvi") |>
  residualize_col("tau_r1_ndwi",   "base_level_ndwi", "resid_tau_r1_ndwi") |>
  residualize_col("tau_var_ndwi",  "base_level_ndwi", "resid_tau_var_ndwi") |>
  residualize_col("tau_skew_ndwi", "base_level_ndwi", "resid_tau_skew_ndwi")

feature_cols <- c("resid_tau_r1_ndvi", "resid_tau_var_ndvi", "resid_tau_skew_ndvi",
                   "resid_tau_r1_ndwi", "resid_tau_var_ndwi", "resid_tau_skew_ndwi")

d <- d |>
  left_join(beetle_flags |> select(pixel_id, any_mortality), by = "pixel_id") |>
  left_join(fire_flags |> select(pixel_id, burned_validation_15_16), by = "pixel_id")
stopifnot(!anyNA(d$any_mortality), !anyNA(d$burned_validation_15_16))

# ----------------------------------------------------------------------------
# Fit + evaluate one outcome: 70/30 split for a held-out AUC/calibration
# check, then refit on all eligible pixels for the production probability.
# AUC computed by the same rank-and-cumsum trapezoidal method as the
# (commented-out) ROC check in 06_clustering.R, applied here to the
# continuous predicted probability instead of a single stress score.
# ----------------------------------------------------------------------------
fit_and_evaluate <- function(d, outcome_col, feature_cols, label, SEED, TRAIN_FRAC) {
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("  Positive rate: %d / %d (%.1f%%)\n",
              sum(d[[outcome_col]]), nrow(d), 100 * mean(d[[outcome_col]])))

  set.seed(SEED)
  train_idx <- sample(nrow(d), size = floor(TRAIN_FRAC * nrow(d)))
  train <- d[train_idx, ]
  test  <- d[-train_idx, ]

  form      <- reformulate(feature_cols, outcome_col)
  fit_train <- glm(form, data = train, family = binomial())
  cat("\nTrain-fold coefficients:\n")
  print(round(coef(summary(fit_train))[, "Estimate"], 4))

  test$predicted_prob <- predict(fit_train, newdata = test, type = "response")

  roc_df <- test |>
    select(predicted_prob, outcome = all_of(outcome_col)) |>
    arrange(desc(predicted_prob)) |>
    mutate(tpr = cumsum(outcome)  / sum(outcome),
           fpr = cumsum(!outcome) / sum(!outcome))
  auc <- sum(diff(roc_df$fpr) * (head(roc_df$tpr, -1) + tail(roc_df$tpr, -1))) / 2
  cat(sprintf("  Held-out AUC (%.0f%% test split): %.3f\n", 100 * (1 - TRAIN_FRAC), auc))

  calib <- test |>
    mutate(bin = ntile(predicted_prob, 10)) |>
    group_by(bin) |>
    summarise(mean_predicted = mean(predicted_prob),
              observed_rate  = mean(.data[[outcome_col]]),
              n              = n(), .groups = "drop")
  cat("  Calibration (decile bins, held-out):\n")
  print(calib)

  fit_full <- glm(form, data = d, family = binomial())
  d$predicted_prob <- predict(fit_full, newdata = d, type = "response")

  list(fit_full = fit_full, auc = auc, calib = calib, d = d)
}

mortality_res <- fit_and_evaluate(d, "any_mortality", feature_cols,
                                   "Beetle mortality (2015-2016)", SEED, TRAIN_FRAC)
fire_res      <- fit_and_evaluate(d, "burned_validation_15_16", feature_cols,
                                   "Fire (2015-2016)", SEED, TRAIN_FRAC)

# ----------------------------------------------------------------------------
# Combined output table + save
# ----------------------------------------------------------------------------
predictions <- d |>
  select(pixel_id, lon, lat) |>
  mutate(
    predicted_mortality_prob = mortality_res$d$predicted_prob,
    predicted_fire_prob      = fire_res$d$predicted_prob
  )

write_csv(predictions, "pixel_approach/data/processed/pixel_regression/pixel_regression_predictions.csv")
cat("\n✓ Saved: pixel_approach/data/processed/pixel_regression/pixel_regression_predictions.csv\n")

# ----------------------------------------------------------------------------
# Probability map: continuous color in place of the categorical cluster
# map (build_cluster_map_panel() in clustering_helpers.R), same boundary
# overlay and axis-label formatting.
# ----------------------------------------------------------------------------
build_prob_map_panel <- function(d, prob_col, sn_boundary, title) {
  boundary_xy <- as.data.frame(sf::st_coordinates(sn_boundary))
  lon_labels <- function(x) sprintf("%g°%s", abs(x), ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
  lat_labels <- function(y) sprintf("%g°%s", abs(y), ifelse(y < 0, "S", ifelse(y > 0, "N", "")))
  ggplot(d, aes(x = lon, y = lat, color = .data[[prob_col]])) +
    geom_point(size = 0.05, alpha = 0.5) +
    geom_polygon(data = boundary_xy, aes(x = X, y = Y), fill = NA, color = "grey20",
                 linewidth = 0.4, inherit.aes = FALSE) +
    scale_color_viridis_c(name = "Predicted\nprobability", option = "inferno") +
    coord_fixed(ratio = 1 / cos(median(d$lat) * pi / 180),
                xlim = range(d$lon), ylim = range(d$lat)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3), labels = lon_labels) +
    scale_y_continuous(labels = lat_labels) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      plot.title.position = "plot",
      plot.title   = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
      legend.title = element_text(face = "bold", color = "#1B2A1F", size = 9),
      axis.title   = element_text(face = "bold", color = "#1B2A1F"),
      axis.text    = element_text(size = 7, color = "#1B2A1F", face = "bold")
    )
}

p_map <- build_prob_map_panel(predictions, "predicted_mortality_prob", sn_boundary,
                               "Predicted beetle mortality probability") +
  build_prob_map_panel(predictions, "predicted_fire_prob", sn_boundary,
                        "Predicted fire probability") +
  plot_annotation(
    title    = "Per-Pixel Predicted Probability (logistic regression on combined NDVI+NDWI tau_r1/tau_var/tau_skew)",
    caption  = "White gaps: pixels burned (MTBS moderate+ severity) at any point 2000-2014, excluded before fitting",
    theme    = theme(plot.title   = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F", size = 10),
                      plot.caption = element_text(hjust = 0.5, color = "#1B2A1F", size = 8, face = "bold",
                                                   margin = margin(t = 10)))
  )

ggsave("pixel_approach/charts/pixel_regression/pixel_regression_prob_map.pdf",
       plot = p_map, width = 14, height = 7, dpi = 150)
cat("✓ Saved: pixel_approach/charts/pixel_regression/pixel_regression_prob_map.pdf\n")

# ----------------------------------------------------------------------------
# Calibration chart: decile-binned held-out predicted vs observed rate,
# per outcome, against the y=x reference line.
# ----------------------------------------------------------------------------
build_calibration_panel <- function(calib, title, auc) {
  # Bin sizes (n) are dropped as a size aesthetic -- ntile() splits the
  # held-out set into deciles that differ by at most 1 pixel, so there's no
  # real variation for size to encode (ggplot's default continuous-size
  # legend was interpolating 5 "pretty" breaks across that 1-pixel range,
  # producing nonsensical fractional-pixel legend entries).
  ggplot(calib, aes(x = mean_predicted, y = observed_rate)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    geom_line(color = "#D73027", linewidth = 0.5) +
    geom_point(color = "#D73027", size = 3) +
    labs(title = title, subtitle = sprintf("Held-out AUC = %.3f", auc),
         x = "Mean predicted probability (decile bin)", y = "Observed rate") +
    theme_bw() +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
          axis.title    = element_text(face = "bold", color = "#1B2A1F"))
}

p_calib <- build_calibration_panel(mortality_res$calib, "Beetle mortality", mortality_res$auc) +
  build_calibration_panel(fire_res$calib, "Fire", fire_res$auc) +
  plot_annotation(
    title = "Calibration: Predicted Probability vs Observed Rate (held-out 30% test split)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F", size = 10))
  )

ggsave("pixel_approach/charts/pixel_regression/pixel_regression_calibration.pdf",
       plot = p_calib, width = 11, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/pixel_regression/pixel_regression_calibration.pdf\n")
