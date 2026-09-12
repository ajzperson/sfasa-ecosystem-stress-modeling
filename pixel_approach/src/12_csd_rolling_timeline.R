library(tidyverse)
library(patchwork)
library(zoo)

# Timeline of the ACTUAL rolling r1/skew statistics per cluster, 2005-2014
# -- not raw NDVI/NDWI level (which mostly reflects baseline vegetation
# density, not CSD), and not a re-aggregated recompute (averaging pixel
# residuals first, then computing rolling stats on the average, gives a
# DIFFERENT answer than the real per-pixel tau_r1 means -- e.g. NDVI
# cluster 1 flips sign). Instead reads the per-pixel rolling series
# 03_pixel_stl_csd.R saves directly (ndvi_rolling_r1.rds etc.) and averages
# THOSE within each cluster -- the same per-pixel quantity tau_r1_mean
# already reports, just shown as its full time series.
#
# tau_var/rolling variance is excluded -- it isn't part of the clustering
# method (replaced by tau_skew), and its own rolling trend didn't
# distinguish cluster 3 either.
#
# Purpose: tau_r1_mean/tau_skew_mean only say a trend exists somewhere
# across 2005-2014 -- they can't say WHEN. This chart checks whether
# cluster 3's rising autocorrelation was already building years before the
# 2012 drought (a real early-warning signal) or only appears after (a
# post-hoc artifact).
#
# A third row -- rolling MEAN clim_resid level, not just its shape (r1/
# skew) -- is included because skew alone can mislead: NDWI's skew
# recovers after 2012 even though its actual residual level keeps getting
# worse through 2014 (2013-2014 was the driest/highest-VPD stretch of the
# whole record, ruling out a genuine wet-spell recovery). Skew measures
# whether a few extreme months stand out against an otherwise-normal
# majority -- once decline becomes widespread, skew can revert toward
# symmetric even as the mean keeps worsening. Unlike r1/skew (nonlinear,
# so per-pixel-then-averaged differs from averaged-then-computed), a MEAN
# is linear -- averaging clim_resid across a cluster first, then taking a
# rolling mean of that single series, is mathematically identical to the
# reverse order, so this is computed cheaply without redoing 03's
# expensive per-pixel rolling-window computation.
#
# Output: pixel_approach/charts/cluster_csd_rolling_timeline.pdf

cluster_cols <- setNames(
  colorRampPalette(c("#2166AC", "#92C5DE", "#F4A582", "#D73027"))(3),
  c("1", "2", "3")
)
clust_labels <- c("1" = "1 - no CSD", "2" = "2 - partial", "3" = "3 - full CSD")

index_config <- list(
  ndvi = list(clusters_path = "pixel_approach/data/processed/ndvi_pixel_clusters.csv",
              mat_path  = "pixel_approach/data/processed/ndvi_matrix.rds",
              r1_path   = "pixel_approach/data/processed/ndvi_rolling_r1.rds",
              skew_path = "pixel_approach/data/processed/ndvi_rolling_skew.rds",
              label     = "NDVI"),
  ndwi = list(clusters_path = "pixel_approach/data/processed/ndwi_pixel_clusters.csv",
              mat_path  = "pixel_approach/data/processed/ndwi_matrix.rds",
              r1_path   = "pixel_approach/data/processed/ndwi_rolling_r1.rds",
              skew_path = "pixel_approach/data/processed/ndwi_rolling_skew.rds",
              label     = "NDWI")
)

WIN <- 54

# Reconstructs clim_resid for every pixel (same method as 03_pixel_stl_csd.R:
# interpolate NAs -> subtract STL seasonal -> subtract climate prediction,
# using the b0/b1/b2/seas_01..12 already saved per pixel), then averages
# across all pixels in each cluster first (cheap, and valid for a linear
# statistic -- see header), then takes a rolling mean of that single series
# over the same WIN=54 windows as r1/skew, for matching date resolution.
build_level_trajectory <- function(cfg, clusters) {
  mat <- readRDS(cfg$mat_path)
  keep_cols <- colnames(mat) <= "2014-12-01"
  mat <- mat[, keep_cols, drop = FALSE]
  col_dates  <- as.Date(colnames(mat))
  col_months <- as.integer(format(col_dates, "%m"))
  gs_idx     <- which(col_months %in% 3:11)
  gs_dates   <- col_dates[gs_idx]

  vpd_mat <- readRDS("pixel_approach/data/processed/vpd_matrix.rds")[, colnames(mat), drop = FALSE][, gs_idx, drop = FALSE]
  ppt_mat <- readRDS("pixel_approach/data/processed/ppt_matrix.rds")[, colnames(mat), drop = FALSE][, gs_idx, drop = FALSE]

  ids <- clusters$pixel_id
  cat(sprintf("  Reconstructing clim_resid for %d pixels...\n", length(ids)))
  mat_sub <- mat[ids, , drop = FALSE]
  mat_sub <- t(apply(mat_sub, 1, function(y) as.numeric(zoo::na.approx(y, na.rm = FALSE, rule = 2))))

  seas_mat <- as.matrix(clusters[, sprintf("seas_%02d", 1:12)])[, match(col_months, 1:12)]
  y_deseas <- mat_sub - seas_mat
  y_gs     <- y_deseas[, gs_idx, drop = FALSE]
  vpd_gs   <- vpd_mat[ids, , drop = FALSE]
  ppt_gs   <- ppt_mat[ids, , drop = FALSE]

  predicted  <- clusters$b0 + clusters$b1 * vpd_gs + clusters$b2 * ppt_gs
  clim_resid <- y_gs - predicted

  window_idx <- lapply(WIN:length(gs_dates), function(j) (j - WIN + 1):j)
  window_end_dates <- as.Date(vapply(window_idx, function(idx) as.character(gs_dates[tail(idx, 1)]), character(1)))

  out <- vector("list", 3)
  for (k in 1:3) {
    rows <- which(clusters$cluster == k)
    # Cluster mean first (cheap), then rolling mean of that single series --
    # identical result to per-pixel-rolling-then-averaging, since mean is linear.
    cluster_mean_series <- colMeans(clim_resid[rows, , drop = FALSE], na.rm = TRUE)
    rolled <- vapply(window_idx, function(idx) mean(cluster_mean_series[idx]), numeric(1))
    out[[k]] <- tibble(date = window_end_dates, cluster = factor(k, levels = 1:3),
                        statistic = "level", value = rolled)
  }
  bind_rows(out)
}

build_cluster_trajectories <- function(cfg) {
  cat(sprintf("\n=== %s ===\n", cfg$label))
  clusters <- read_csv(cfg$clusters_path, show_col_types = FALSE) |>
    select(pixel_id, cluster, b0, b1, b2, starts_with("seas_"))
  cat(sprintf("✓ %d pixels with cluster assignment\n", nrow(clusters)))

  r1_mat   <- readRDS(cfg$r1_path)
  skew_mat <- readRDS(cfg$skew_path)
  dates    <- as.Date(colnames(r1_mat))

  # Row i of each matrix corresponds to pixel_id == i -- 03_pixel_stl_csd.R
  # builds these via future_lapply(seq_len(n_pixels), ...), same order as
  # pixel_meta$pixel_id (which is 1:n_pixels).
  build_one <- function(mat, stat_name) {
    out <- vector("list", 3)
    for (k in 1:3) {
      rows <- clusters$pixel_id[clusters$cluster == k]
      out[[k]] <- tibble(
        date      = dates,
        cluster   = factor(k, levels = 1:3),
        statistic = stat_name,
        value     = colMeans(mat[rows, , drop = FALSE], na.rm = TRUE)
      )
    }
    bind_rows(out)
  }

  cat("  Building level trajectory (reconstructing clim_resid)...\n")
  level_traj <- build_level_trajectory(cfg, clusters)

  traj <- bind_rows(
    build_one(r1_mat,   "r1"),
    build_one(skew_mat, "skew"),
    level_traj
  )

  cat("✓ Cluster trajectories built (level/r1/skew x 3 clusters x 82 dates)\n")
  traj
}

build_panel <- function(traj, stat_name, ylab, label) {
  d <- traj |> filter(statistic == stat_name)
  ggplot(d, aes(x = date, y = value, color = cluster)) +
    geom_vline(xintercept = as.numeric(as.Date("2012-01-01")),
               linetype = "dashed", color = "grey40", linewidth = 0.4) +
    annotate("text", x = as.Date("2012-02-01"), y = Inf, label = "Drought onset",
             hjust = 0, vjust = 1.5, size = 2.6, fontface = "bold", color = "#1B2A1F") +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = cluster_cols, name = "Stress level", labels = clust_labels) +
    labs(title = label, x = NULL, y = ylab) +
    theme_bw() +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
          axis.title  = element_text(face = "bold", color = "#1B2A1F"),
          axis.text   = element_text(face = "bold", color = "#1B2A1F"),
          legend.title = element_text(face = "bold", color = "#1B2A1F"),
          legend.text  = element_text(face = "bold", color = "#1B2A1F"))
}

results <- lapply(index_config, build_cluster_trajectories)

panels <- list()
for (idx in names(index_config)) {
  cfg  <- index_config[[idx]]
  traj <- results[[idx]]
  panels[[paste0(idx, "_level")]] <- build_panel(traj, "level", "Rolling mean clim_resid (level)", cfg$label)
  panels[[paste0(idx, "_r1")]]    <- build_panel(traj, "r1",    "Rolling r1 (autocorrelation)", cfg$label)
  panels[[paste0(idx, "_skew")]]  <- build_panel(traj, "skew",  "Rolling skewness", cfg$label)
}

p <- (panels$ndvi_level | panels$ndwi_level) /
     (panels$ndvi_r1 | panels$ndwi_r1) /
     (panels$ndvi_skew | panels$ndwi_skew) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Rolling CSD Statistics by Cluster, 2005-2014",
    subtitle = "Per-pixel rolling level/r1/skew (WIN=54mo) averaged within each cluster -- level is the real check on whether r1/skew's shape changes reflect actual worsening",
    theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
                      plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"))
  )

out_file <- "pixel_approach/charts/cluster_csd_rolling_timeline.pdf"
ggsave(out_file, plot = p, width = 12, height = 12, dpi = 150)
cat(sprintf("\n✓ Saved: %s\n", out_file))
