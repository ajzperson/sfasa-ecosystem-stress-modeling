library(tidyverse)
library(sf)
library(cluster)

# Shared logic between 06_clustering.R (per-index) and
# 06b_clustering_combined.R (combined NDVI+NDWI) -- the two scripts differ
# in how many features they cluster on (2 vs 4) and relabel clusters by
# stress rank the same way (unweighted composite sum of the residualized
# features), but share everything else: baseline-density residualization,
# the K-sweep loop, and the two chart builders. Sourced by both, not run
# standalone.

# --- Cluster color/label style, shared by every K-means chart in this repo ---
build_cluster_style <- function(K) {
  cols <- setNames(
    colorRampPalette(c("#2166AC", "#92C5DE", "#F4A582", "#D73027"))(K),
    as.character(seq_len(K))
  )
  # K=1 is degenerate (one cluster can't be both least and most stressed),
  # so it gets its own generic label. K=4's two middle tiers get distinct
  # labels (low-/high-moderate); any other K just repeats "Moderate".
  mid_labels <- if (K == 4) {
    c(" - Low-moderate", " - High-moderate")
  } else {
    rep(" - Moderate", max(K - 2, 0))
  }
  labels <- if (K == 1) {
    "1 - all"
  } else {
    paste0(seq_len(K), c(" - Least-stressed", mid_labels, " - Most-stressed"))
  }
  list(cols = cols, labels = labels)
}

# --- Baseline vegetation density (mean pre-2015 level), per pixel ---
# tau_r1/tau_skew both correlate with how much forest a pixel has to begin
# with, which otherwise leaks into cluster assignment -- this is what gets
# residualized out below.
compute_base_level <- function(mat_path, pixel_ids, col_name) {
  mat       <- readRDS(mat_path)
  mat_dates <- as.Date(colnames(mat))
  out <- tibble(
    pixel_id = pixel_ids,
    !!col_name := rowMeans(mat[, mat_dates < as.Date("2015-01-01")], na.rm = TRUE)
  )
  rm(mat); gc(verbose = FALSE)
  out
}

# Adds resid_<tau_r1_col>/resid_<tau_skew_col>-equivalent columns to `d` by
# regressing each against `base_level_col`. Caller supplies the output
# column names so 06 (unsuffixed: resid_tau_r1) and 06b (suffixed per
# index: resid_tau_r1_ndvi) can both use this unchanged.
residualize_tau <- function(d, tau_r1_col, tau_skew_col, base_level_col,
                            out_r1_col, out_skew_col) {
  d[[out_r1_col]]   <- resid(lm(reformulate(base_level_col, tau_r1_col),   data = d))
  d[[out_skew_col]] <- resid(lm(reformulate(base_level_col, tau_skew_col), data = d))
  d
}

# --- K-means sweep (K=2-8) for the elbow/silhouette charts, plus the fit
# at the chosen K (reused instead of re-running kmeans() a second time) ---
run_kmeans_sweep <- function(feat, K, SEED, k_range = 2:8,
                             nstart = 25, iter.max = 100) {
  wcss_vec <- numeric(length(k_range))
  sil_vec  <- numeric(length(k_range))

  # Subsample for silhouette (full distance matrix on ~250k pixels is too large)
  set.seed(SEED)
  sil_idx  <- sample(nrow(feat), min(5000, nrow(feat)))
  feat_sil <- feat[sil_idx, ]

  km_final <- NULL
  for (ki in seq_along(k_range)) {
    k_i <- k_range[ki]
    set.seed(SEED)
    km_i <- kmeans(feat, centers = k_i, nstart = nstart, iter.max = iter.max)
    wcss_vec[ki] <- km_i$tot.withinss
    sil_vec[ki]  <- mean(silhouette(km_i$cluster[sil_idx], dist(feat_sil))[, 3])
    cat(sprintf("  k=%d  WCSS=%.1f  Silhouette=%.5f\n", k_i, wcss_vec[ki], sil_vec[ki]))
    if (k_i == K) km_final <- km_i
  }
  sweep_df <- data.frame(k = k_range, wcss = wcss_vec, silhouette = sil_vec)

  if (is.null(km_final)) {
    # K fell outside the swept range, so it wasn't computed above
    set.seed(SEED)
    km_final <- kmeans(feat, centers = K, nstart = nstart, iter.max = iter.max)
  }

  list(sweep_df = sweep_df, km_final = km_final)
}

# --- Elbow + silhouette panel (one "row"), titled per caller ---
build_k_sweep_panel <- function(sweep_df, K, title_prefix) {
  k_range <- sweep_df$k

  p_elbow <- ggplot(sweep_df, aes(x = k, y = wcss)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    geom_point(data = sweep_df[sweep_df$k == K, ],
               size = 5, shape = 21, fill = "#D73027", color = "black") +
    scale_x_continuous(breaks = k_range) +
    labs(title = sprintf("%s - Elbow Method", title_prefix),
         x = "Number of clusters (K)", y = "Total within-cluster SS") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
      axis.title = element_text(face = "bold", color = "#1B2A1F"),
      axis.text  = element_text(face = "bold", size = 11, color = "#1B2A1F")
    )

  p_sil <- ggplot(sweep_df, aes(x = k, y = silhouette)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    geom_point(data = sweep_df[sweep_df$k == K, ],
               size = 5, shape = 21, fill = "#D73027", color = "black") +
    scale_x_continuous(breaks = k_range) +
    labs(title = sprintf("%s - Silhouette Score", title_prefix),
         x = "Number of clusters (K)", y = "Mean silhouette width") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
      axis.title = element_text(face = "bold", color = "#1B2A1F"),
      axis.text  = element_text(face = "bold", size = 11, color = "#1B2A1F")
    )

  p_elbow + p_sil
}

# --- Cluster map panel, titled per caller. `d` must have lon/lat/cluster. ---
# Plain polygon on coord_fixed(), not coord_sf() -- its graticule labels
# ignore theme's bold face, so degree labels are recreated manually.
build_cluster_map_panel <- function(d, sn_boundary, cluster_cols, clust_labels, title) {
  boundary_xy <- as.data.frame(sf::st_coordinates(sn_boundary))
  lon_labels <- function(x) sprintf("%g°%s", abs(x), ifelse(x < 0, "W", ifelse(x > 0, "E", "")))
  lat_labels <- function(y) sprintf("%g°%s", abs(y), ifelse(y < 0, "S", ifelse(y > 0, "N", "")))
  ggplot(d, aes(x = lon, y = lat, color = cluster)) +
    geom_point(size = 0.05, alpha = 0.5) +
    geom_polygon(data = boundary_xy, aes(x = X, y = Y), fill = NA, color = "grey20",
                 linewidth = 0.4, inherit.aes = FALSE) +
    scale_color_manual(values = cluster_cols, name = "Stress level",
                       labels = clust_labels) +
    coord_fixed(ratio = 1 / cos(median(d$lat) * pi / 180),
                xlim = range(d$lon), ylim = range(d$lat)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3), labels = lon_labels) +
    scale_y_continuous(labels = lat_labels) +
    labs(title = title, x = "Longitude", y = "Latitude",
         caption = "White gaps: pixels burned (MTBS moderate+ severity) at any point\n2000-2014, excluded before clustering") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_bw() +
    theme(
      # "plot" (not the default "panel") so a long title centers over the
      # full device width, including the legend's horizontal space -- at
      # "panel" a long title centers over just the (legend-narrowed) panel
      # and overflows past the left edge of the device, where it's clipped
      # (no such room on the right, which just extends harmlessly into the
      # legend's own whitespace).
      plot.title.position = "plot",
      plot.title      = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
      plot.caption    = element_text(hjust = 0.5, color = "#1B2A1F", size = 8, face = "bold",
                                      margin = margin(t = 10)),
      legend.position = "right",
      legend.title    = element_text(face = "bold", color = "#1B2A1F"),
      legend.text     = element_text(face = "bold", color = "#1B2A1F"),
      axis.title      = element_text(face = "bold", color = "#1B2A1F"),
      axis.text       = element_text(size = 7, color = "#1B2A1F", face = "bold")
    )
}
