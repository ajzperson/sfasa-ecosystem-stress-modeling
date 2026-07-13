library(tidyverse)
library(sf)
library(patchwork)
library(multcompView)
library(ggsignif)
library(cluster)

# K-means clustering on tau_r1 + tau_var, run independently for NDVI and
# NDWI (not combined into one 4-feature clustering) — the point is to
# cross-validate two independent stress signals against each other, which
# only works if each keeps its own cluster assignment to compare. tau_skew
# stays excluded from clustering for both: it was already the weakest
# NDVI predictor, and NDVI's and NDWI's tau_skew values are actually
# negatively correlated (r=-0.12) — not just noisy, but inconsistent
# across indices.
#
# Each index is validated against its own actual-drop ground truth
# (Actual_Drop_15_16 for NDVI, Actual_NDWI_Drop_15_16 for NDWI) rather than
# cross-validated against the other's, so tau is always tested against that
# same variable's own future collapse — matching standard CSD-literature
# methodology and keeping the two analyses symmetric/comparable.
#
# Outputs:
#   pixel_approach/data/processed/ndvi_pixel_clusters.csv
#   pixel_approach/data/processed/ndwi_pixel_clusters.csv
#   pixel_approach/charts/cluster_k_sweep.pdf         (NDVI/NDWI stacked)
#   pixel_approach/charts/cluster_map.pdf              (NDVI/NDWI side by side)
#   pixel_approach/charts/cluster_scatter_pairwise.pdf (NDVI/NDWI side by side)
#   pixel_approach/charts/cluster_timelines.pdf   (NDVI/NDWI side by side, independent y-scales)
#   pixel_approach/charts/cluster_validation.pdf       (NDVI/NDWI side by side, independent y-scales)
#   pixel_approach/charts/cluster_roc.pdf              (NDVI vs NDWI overlaid — both axes already unitless [0,1])

if (!exists("K"))       K       <- 3
if (!exists("SEED"))    SEED    <- 42
if (!exists("INDICES")) INDICES <- c("ndvi", "ndwi")

dir.create("pixel_approach/charts", showWarnings = FALSE, recursive = TRUE)

index_config <- list(
  ndvi = list(
    tau_path     = "pixel_approach/data/processed/ndvi_pixel_tau.csv",
    fire_path    = "pixel_approach/data/processed/pixel_fire_ndvi.csv",
    mat_path     = "pixel_approach/data/processed/ndvi_matrix.rds",
    drop_col     = "Actual_Drop_15_16",
    clusters_out = "pixel_approach/data/processed/ndvi_pixel_clusters.csv",
    label        = "NDVI"
  ),
  ndwi = list(
    tau_path     = "pixel_approach/data/processed/ndwi_pixel_tau.csv",
    fire_path    = "pixel_approach/data/processed/pixel_fire_ndwi.csv",
    mat_path     = "pixel_approach/data/processed/ndwi_matrix.rds",
    drop_col     = "Actual_NDWI_Drop_15_16",
    clusters_out = "pixel_approach/data/processed/ndwi_pixel_clusters.csv",
    label        = "NDWI"
  )
)

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

# Colors match heatmap convention: blue (rank 1, least stressed) -> red
# (rank K, most stressed). Shared across both indices' charts since the
# cluster labels (1..K, "no CSD".."full CSD") mean the same thing in both.
cluster_cols <- setNames(
  colorRampPalette(c("#2166AC", "#92C5DE", "#F4A582", "#D73027"))(K),
  as.character(seq_len(K))
)
# rep(" - partial", K - 2) errors for K=1 (negative times); K=1 is
# degenerate anyway (one cluster can't be both least and most stressed),
# so it gets its own generic label rather than trying to force the
# no-CSD/partial/full-CSD framing onto it.
clust_labels <- if (K == 1) {
  "1 - all"
} else {
  paste0(seq_len(K), c(" - no CSD", rep(" - partial", max(K - 2, 0)),
                       " - full CSD"))
}

# ----------------------------------------------------------------------------
# Core clustering step: K-sweep, k-means at K, stress-rank relabeling, save
# CSV. Generic over cfg (tau_path/drop_col/etc.) so the same code produces
# both ndvi_pixel_clusters.csv and ndwi_pixel_clusters.csv.
# ----------------------------------------------------------------------------
run_clustering <- function(cfg, K, SEED) {
  cat(sprintf("\n=== %s ===\n", cfg$label))
  cat(sprintf("Loading %s...\n", cfg$tau_path))
  d <- read_csv(cfg$tau_path, show_col_types = FALSE) |>
    filter(!is.na(tau_r1), !is.na(tau_var), !is.na(tau_skew))
  cat(sprintf("✓ %d pixels with valid tau_r1, tau_var, tau_skew\n", nrow(d)))

  # Exclude pixels burned at any point during the 2000-2014 baseline window --
  # a fire causes a sudden, discrete NDVI/NDWI crash unrelated to gradual
  # drought stress, which would contaminate the tau_r1/tau_var trend
  # statistics computed over that same window. Verified against real data
  # (09_fire_validation.R) that even the earliest possible fires (2000-2003)
  # still show a persistent NDVI deficit through 2014 -- no recovery-time
  # cutoff is safe within this dataset's window, so any fire excludes.
  fire_flags <- read_csv(cfg$fire_path, show_col_types = FALSE) |>
    select(pixel_id, burned_baseline_00_14)
  n_before <- nrow(d)
  d <- d |>
    left_join(fire_flags, by = "pixel_id") |>
    filter(!burned_baseline_00_14) |>
    select(-burned_baseline_00_14)
  cat(sprintf("✓ Excluded %d pixels burned during 2000-2014 baseline (%d remain)\n",
              n_before - nrow(d), nrow(d)))

  # Scale features before clustering (Euclidean distance is scale-sensitive)
  feat <- d |> select(tau_r1, tau_var) |> scale()

  # --- K sweep - elbow + silhouette (K = 2 to 8) ---
  cat("Running K sweep (K=2 to 8) for elbow and silhouette plots...\n")
  k_range  <- 2:8
  wcss_vec <- numeric(length(k_range))
  sil_vec  <- numeric(length(k_range))

  # Subsample for silhouette (full distance matrix on 248k pixels is too large)
  set.seed(SEED)
  sil_idx  <- sample(nrow(feat), min(5000, nrow(feat)))
  feat_sil <- feat[sil_idx, ]

  # Captures the sweep's own fit at k_i == K (same feat/seed/nstart/iter.max
  # as the "final" fit below, so it's deterministically identical) to avoid
  # running kmeans() a second time for the same K.
  km_final <- NULL
  for (ki in seq_along(k_range)) {
    k_i <- k_range[ki]
    set.seed(SEED)
    km_i <- kmeans(feat, centers = k_i, nstart = 25, iter.max = 100)
    wcss_vec[ki] <- km_i$tot.withinss
    sil_vec[ki]  <- mean(silhouette(km_i$cluster[sil_idx],
                                    dist(feat_sil))[, 3])
    cat(sprintf("  k=%d  WCSS=%.1f  Silhouette=%.3f\n",
                k_i, wcss_vec[ki], sil_vec[ki]))
    if (k_i == K) km_final <- km_i
  }
  sweep_df <- data.frame(k = k_range, wcss = wcss_vec, silhouette = sil_vec)

  # --- K-means at chosen K ---
  cat(sprintf("Running k-means (k=%d, seed=%d)...\n", K, SEED))
  if (is.null(km_final)) {
    # K fell outside the swept 2:8 range, so it wasn't computed above
    set.seed(SEED)
    km_final <- kmeans(feat, centers = K, nstart = 25, iter.max = 100)
  }
  km <- km_final
  d$cluster_raw <- km$cluster

  # Re-label clusters by composite stress rank so colors are semantically
  # consistent with the heatmaps (red = more stressed, blue = less stressed).
  stress_rank <- d |>
    group_by(cluster_raw) |>
    summarise(composite = mean(tau_r1) + mean(tau_var), .groups = "drop") |>
    arrange(composite) |>
    mutate(stress_rank = row_number())   # 1 = least stressed, K = most stressed

  d <- d |>
    left_join(stress_rank |> select(cluster_raw, stress_rank),
              by = "cluster_raw") |>
    mutate(cluster = factor(stress_rank, levels = seq_len(K)))

  cat(sprintf("✓ Cluster sizes (stress rank 1=least → %d=most): %s\n",
              K, paste(table(d$cluster), collapse = " / ")))

  # Print cluster centroids in original (unscaled) units
  drop_col <- cfg$drop_col
  centers_orig <- d |>
    group_by(cluster) |>
    summarise(
      n           = n(),
      tau_r1      = round(mean(tau_r1),   4),
      tau_var     = round(mean(tau_var),  4),
      tau_skew    = round(mean(tau_skew), 4),
      composite   = round(mean(tau_r1) + mean(tau_var), 4),
      actual_drop = round(mean(.data[[drop_col]], na.rm = TRUE), 4),
      .groups = "drop"
    ) |>
    arrange(as.integer(cluster))
  cat("\nCluster profile (1=least stressed, K=most stressed):\n")
  print(centers_orig)

  write_csv(d, cfg$clusters_out)
  cat(sprintf("\n✓ Saved: %s\n", cfg$clusters_out))

  list(d = d, sweep_df = sweep_df, cfg = cfg)
}

results <- setNames(
  lapply(INDICES, function(idx) run_clustering(index_config[[idx]], K, SEED)),
  INDICES
)

# Combines one panel per index into a single patchwork figure — side by
# side (or stacked, for combine_op = `/`) when both indices are present,
# returned unchanged when only one is. Every chart below is built from
# this so it genuinely follows INDICES instead of assuming both NDVI and
# NDWI always ran (an earlier version hardcoded results$ndvi/results$ndwi
# directly, which crashed if INDICES had only one entry).
combine_panels <- function(panels, guides = "collect", combine_op = `+`) {
  if (length(panels) == 1) return(panels[[1]])
  Reduce(combine_op, panels) + plot_layout(guides = guides)
}

# Chart dimensions below were tuned for 2 side-by-side panels — scale by
# how many indices are actually being plotted so a single-index run
# doesn't produce an oddly stretched-wide (or, for the stacked k-sweep
# chart, stretched-tall) single panel.
n_idx <- length(INDICES)

# ----------------------------------------------------------------------------
# K-sweep chart: elbow + silhouette per index, stacked (NDVI row, NDWI row)
# ----------------------------------------------------------------------------
build_k_sweep_row <- function(res, K) {
  sweep_df <- res$sweep_df
  k_range  <- sweep_df$k

  p_elbow <- ggplot(sweep_df, aes(x = k, y = wcss)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    geom_point(data = sweep_df[sweep_df$k == K, ],
               size = 5, shape = 21, fill = "#D73027", color = "black") +
    scale_x_continuous(breaks = k_range) +
    labs(title = sprintf("%s - Elbow Method", res$cfg$label),
         x = "Number of clusters (K)", y = "Total within-cluster SS") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p_sil <- ggplot(sweep_df, aes(x = k, y = silhouette)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    geom_point(data = sweep_df[sweep_df$k == K, ],
               size = 5, shape = 21, fill = "#D73027", color = "black") +
    scale_x_continuous(breaks = k_range) +
    labs(title = sprintf("%s - Silhouette Score", res$cfg$label),
         x = "Number of clusters (K)", y = "Mean silhouette width") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  p_elbow + p_sil
}

p_sweep <- combine_panels(
  lapply(INDICES, function(idx) build_k_sweep_row(results[[idx]], K)),
  combine_op = `/`
) +
  plot_annotation(
    title = sprintf("K-means Hyperparameter Sweep (K = 2-8); selected K = %d shown in red", K),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

ggsave("pixel_approach/charts/cluster_k_sweep.pdf",
       plot = p_sweep, width = 10, height = 4 * n_idx, dpi = 150)
cat("\n✓ Saved: pixel_approach/charts/cluster_k_sweep.pdf\n")

# ----------------------------------------------------------------------------
# Cluster map: NDVI / NDWI side by side, shared legend (same cluster labels)
# ----------------------------------------------------------------------------
build_map_panel <- function(res) {
  d <- res$d
  ggplot(d, aes(x = lon, y = lat, color = cluster)) +
    geom_point(size = 0.05, alpha = 0.5) +
    geom_sf(data = sn_boundary, fill = NA, color = "grey20",
            linewidth = 0.4, inherit.aes = FALSE) +
    scale_color_manual(values = cluster_cols, name = "Stress level",
                       labels = clust_labels) +
    coord_sf(xlim = range(d$lon), ylim = range(d$lat)) +
    labs(title = res$cfg$label, x = "Longitude", y = "Latitude") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_bw() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      axis.text       = element_text(size = 7)
    )
}

p_map <- combine_panels(lapply(INDICES, function(idx) build_map_panel(results[[idx]]))) +
  plot_annotation(
    title    = sprintf("K-means Clusters (k=%d) on tau_r1 + tau_var", K),
    subtitle = "Color: blue = least stressed -> red = most stressed (consistent with heatmaps)",
    theme    = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

ggsave("pixel_approach/charts/cluster_map.pdf",
       plot = p_map, width = 7 * n_idx, height = 7, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_map.pdf\n")

# ----------------------------------------------------------------------------
# Scatter: tau_r1 vs tau_var, NDVI / NDWI side by side, shared legend
# ----------------------------------------------------------------------------
build_scatter_panel <- function(res) {
  # Every other random draw in this script (K-sweep silhouette sampling,
  # kmeans, spatial thinning) is seeded via SEED for reproducibility —
  # this one was missing it, so the specific 30k points plotted differed
  # run to run even though nothing else did.
  set.seed(SEED)
  samp <- res$d |> slice_sample(n = min(30000, nrow(res$d)))
  ggplot(samp, aes(x = tau_r1, y = tau_var, color = cluster)) +
    geom_point(size = 0.3, alpha = 0.4) +
    scale_color_manual(values = cluster_cols, name = "Stress level",
                       labels = clust_labels) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = res$cfg$label, x = "tau_r1 (autocorrelation)", y = "tau_var") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

p_scatter <- combine_panels(lapply(INDICES, function(idx) build_scatter_panel(results[[idx]]))) +
  plot_annotation(
    title = "Cluster Separation - tau_r1 vs tau_var",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

ggsave("pixel_approach/charts/cluster_scatter_pairwise.pdf",
       plot = p_scatter, width = 6 * n_idx, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_scatter_pairwise.pdf\n")

# ----------------------------------------------------------------------------
# Per-cluster timelines: NDVI / NDWI side by side. Independent y-scales —
# NDVI and NDWI are different physical units, same reasoning as the
# Actual_Drop heatmap in 05 (what matters is each index's own trajectory
# shape, not whether the raw values line up numerically).
# ----------------------------------------------------------------------------
build_timeline_panel <- function(res) {
  cat(sprintf("Loading %s matrix for per-cluster timelines...\n", res$cfg$label))
  mat <- readRDS(res$cfg$mat_path)
  d   <- res$d
  col_dates <- as.Date(colnames(mat))

  timeline_list <- vector("list", K)
  for (k in seq_len(K)) {
    idx <- d$pixel_id[d$cluster == k]
    sub <- mat[idx, , drop = FALSE]
    monthly <- data.frame(date = col_dates,
                          value = colMeans(sub, na.rm = TRUE),
                          cluster = factor(k))
    timeline_list[[k]] <- monthly |>
      mutate(year  = as.integer(format(date, "%Y")),
             month = as.integer(format(date, "%m"))) |>
      filter(month %in% 3:11) |>       # growing season only (Mar–Nov)
      group_by(cluster, year) |>
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  }
  timelines <- bind_rows(timeline_list)

  ggplot(timelines, aes(x = year, y = value, color = cluster)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.5) +
    scale_color_manual(values = cluster_cols, name = "Stress level",
                       labels = clust_labels) +
    geom_vline(xintercept = 2012, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 2015, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    annotate("text", x = 2012.2, y = Inf, label = "Drought onset",
             hjust = 0, vjust = 1.5, size = 2.8, color = "grey40") +
    annotate("text", x = 2015.2, y = Inf, label = "Peak mortality",
             hjust = 0, vjust = 1.5, size = 2.8, color = "grey40") +
    labs(
      title = res$cfg$label,
      x = NULL, y = sprintf("Annual Mean %s", res$cfg$label)
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

p_ts <- combine_panels(lapply(INDICES, function(idx) build_timeline_panel(results[[idx]]))) +
  plot_annotation(
    title = "Annual Mean Trajectory by Cluster (2000-2016)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

ggsave("pixel_approach/charts/cluster_timelines.pdf",
       plot = p_ts, width = 7 * n_idx, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_timelines.pdf\n")

# ----------------------------------------------------------------------------
# ANOVA + Tukey HSD on each index's own actual-drop column by cluster.
# Spatial thinning (1 pixel per 0.04° x 0.04° cell, ~4km spacing) avoids
# inflating the F-statistic via autocorrelation between neighboring 500m
# pixels. Boxplot shows all pixels; ANOVA uses the thinned subset.
# ----------------------------------------------------------------------------
run_validation <- function(res, SEED) {
  d        <- res$d
  drop_col <- res$cfg$drop_col
  val_data <- d |> filter(!is.na(.data[[drop_col]]))

  set.seed(SEED)
  val_thin <- val_data |>
    mutate(lon_bin = round(lon / 0.04) * 0.04,
           lat_bin = round(lat / 0.04) * 0.04) |>
    group_by(lon_bin, lat_bin) |>
    slice_sample(n = 1) |>
    ungroup()

  cat(sprintf("\n=== %s: spatial thinning: %d → %d pixels (%.0f%%, >=4km spacing) ===\n",
              res$cfg$label, nrow(val_data), nrow(val_thin),
              100 * nrow(val_thin) / nrow(val_data)))

  # bquote()/.() substitutes the real column name in as a symbol (not just
  # a string that evaluates to one), so aov_fit$call — and therefore
  # TukeyHSD's printed "Fit:" line — shows a clean `Actual_Drop_15_16 ~
  # cluster` instead of a literal val_thin[[drop_col]] indexing expression
  # or an as.formula(paste(...)) call. as.formula(paste(drop_col, "~
  # cluster")) fits correctly too but still prints the paste() expression
  # itself, since R's stored call reflects the calling code, not its value.
  aov_fit <- eval(bquote(aov(.(as.name(drop_col)) ~ cluster, data = val_thin)))
  cat(sprintf("=== One-way ANOVA (spatially thinned): %s ~ cluster ===\n", drop_col))
  print(summary(aov_fit))

  tukey_result <- TukeyHSD(aov_fit)
  cat("\n=== Tukey HSD pairwise comparisons ===\n")
  print(tukey_result)

  mean_df <- val_data |>
    group_by(cluster) |>
    summarise(mean_val = mean(.data[[drop_col]], na.rm = TRUE), .groups = "drop") |>
    mutate(label = sprintf("mean = %.3f", mean_val))

  format_p <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
  tk <- tukey_result$cluster
  # geom_signif comparisons order (1,2)/(2,3)/(1,3) vs. TukeyHSD row order (2-1)/(3-1)/(3-2)
  sig_labels <- c(format_p(tk["2-1", "p adj"]),
                   format_p(tk["3-2", "p adj"]),
                   format_p(tk["3-1", "p adj"]))
  overall_p <- summary(aov_fit)[[1]][["Pr(>F)"]][1]

  p_val <- ggplot(val_data, aes(x = cluster, y = .data[[drop_col]], fill = cluster)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, width = 0.5) +
    geom_text(data = mean_df, aes(x = cluster, y = mean_val, label = label),
              inherit.aes = FALSE, size = 3.2, vjust = -0.6, fontface = "italic") +
    geom_signif(
      comparisons  = list(c("1", "2"), c("2", "3"), c("1", "3")),
      annotations  = sig_labels,
      step_increase = 0.12,
      tip_length   = 0.01,
      textsize     = 3.2
    ) +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    labs(
      title    = res$cfg$label,
      subtitle = sprintf("ANOVA F = %.1f, %s (n = %s, >=4km spacing) - Tukey HSD shown",
                         summary(aov_fit)[[1]][["F value"]][1],
                         format_p(overall_p),
                         format(nrow(val_thin), big.mark = ",")),
      x = "Cluster", y = sprintf("%s (baseline mean - 2015-2016 min)", drop_col)
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9)
    )

  list(plot = p_val, aov_fit = aov_fit, val_thin = val_thin)
}

val_results <- setNames(
  lapply(INDICES, function(idx) run_validation(results[[idx]], SEED)),
  INDICES
)

p_val <- combine_panels(lapply(INDICES, function(idx) val_results[[idx]]$plot)) +
  plot_annotation(
    title = "Actual Drop 2015-2016 by Cluster",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

ggsave("pixel_approach/charts/cluster_validation.pdf",
       plot = p_val, width = 7 * n_idx, height = 6, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_validation.pdf\n")

# ----------------------------------------------------------------------------
# Binary classification validation - confusion matrix + ROC-AUC on each
# index's own actual-drop ground truth. ROC axes (TPR/FPR) are already
# unitless [0,1] regardless of source index, so NDVI and NDWI are overlaid
# on one plot rather than shown as separate panels.
# ----------------------------------------------------------------------------
compute_roc <- function(res, K) {
  drop_col <- res$cfg$drop_col
  clf_data <- res$d |>
    filter(!is.na(.data[[drop_col]])) |>
    mutate(
      # Ground truth: top 20% of actual drop = "severe collapse"
      severe    = .data[[drop_col]] >= quantile(.data[[drop_col]], 0.80, na.rm = TRUE),
      # Prediction: cluster K (full CSD) = "predicted severe"
      predicted = as.integer(cluster) == K
    )

  TP <- sum( clf_data$severe &  clf_data$predicted)
  FP <- sum(!clf_data$severe &  clf_data$predicted)
  TN <- sum(!clf_data$severe & !clf_data$predicted)
  FN <- sum( clf_data$severe & !clf_data$predicted)

  sensitivity <- TP / (TP + FN)
  precision   <- TP / (TP + FP)
  specificity <- TN / (TN + FP)
  fpr_binary  <- FP / (FP + TN)
  accuracy    <- (TP + TN) / (TP + FP + TN + FN)
  lift        <- precision / 0.20   # baseline rate = 20% by definition

  cat(sprintf("\n=== %s: Binary classification, Cluster %d vs top-20%% %s ===\n",
              res$cfg$label, K, drop_col))
  cat(sprintf("  Threshold: top 20%% (>= %.4f)\n",
              quantile(clf_data[[drop_col]], 0.80, na.rm = TRUE)))
  cat(sprintf("  TP=%d  FP=%d  TN=%d  FN=%d\n", TP, FP, TN, FN))
  cat(sprintf("  Sensitivity (Recall): %.3f\n", sensitivity))
  cat(sprintf("  Precision:            %.3f\n", precision))
  cat(sprintf("  Specificity:          %.3f\n", specificity))
  cat(sprintf("  False Positive Rate:  %.3f\n", fpr_binary))
  cat(sprintf("  Accuracy:             %.3f\n", accuracy))
  cat(sprintf("  Lift factor:          %.2fx  (vs 1.00x random)\n", lift))

  # ROC-AUC via trapezoidal rule on composite stress score
  roc_df <- clf_data |>
    mutate(stress_score = tau_r1 + tau_var) |>
    select(stress_score, severe) |>
    arrange(desc(stress_score)) |>
    mutate(
      tpr   = cumsum(severe)  / sum(severe),
      fpr   = cumsum(!severe) / sum(!severe),
      index = res$cfg$label
    )

  auc <- sum(diff(roc_df$fpr) * (head(roc_df$tpr, -1) + tail(roc_df$tpr, -1))) / 2
  cat(sprintf("  ROC-AUC:              %.3f\n", auc))

  list(roc_df = roc_df, auc = auc)
}

roc_results <- setNames(
  lapply(INDICES, function(idx) compute_roc(results[[idx]], K)),
  INDICES
)

roc_combined <- bind_rows(lapply(roc_results, function(r) r$roc_df))
roc_colors   <- c(NDVI = "#D73027", NDWI = "#2166AC")
auc_label    <- paste(sprintf("%s AUC = %.3f",
                              sapply(INDICES, function(idx) index_config[[idx]]$label),
                              sapply(roc_results, function(r) r$auc)),
                      collapse = "\n")

p_roc <- ggplot(roc_combined, aes(x = fpr, y = tpr, color = index)) +
  geom_line(linewidth = 0.8) +
  geom_abline(linetype = "dashed", color = "grey50") +
  scale_color_manual(values = roc_colors, name = NULL) +
  annotate("text", x = 0.55, y = 0.15,
           label = auc_label,
           size = 4, fontface = "bold", hjust = 0) +
  labs(
    title    = "ROC Curve: CSD Stress Score vs Severe Canopy Collapse",
    subtitle = "Ground truth = top 20% actual drop (own index); predictor = tau_r1 + tau_var",
    x = "False Positive Rate", y = "True Positive Rate (Sensitivity)"
  ) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))

ggsave("pixel_approach/charts/cluster_roc.pdf",
       plot = p_roc, width = 6.5, height = 5.5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_roc.pdf\n")

# ----------------------------------------------------------------------------
# Summary table, per index
# ----------------------------------------------------------------------------
for (idx in INDICES) {
  res <- results[[idx]]
  drop_col <- res$cfg$drop_col
  cat(sprintf("\n=== %s: final cluster profile (1=least stressed, K=most stressed) ===\n",
              res$cfg$label))
  res$d |>
    group_by(cluster) |>
    summarise(
      n              = n(),
      tau_r1_mean    = round(mean(tau_r1),   3),
      tau_var_mean   = round(mean(tau_var),  3),
      tau_skew_mean  = round(mean(tau_skew), 3),
      actual_mean    = round(mean(.data[[drop_col]],   na.rm = TRUE), 3),
      actual_median  = round(median(.data[[drop_col]], na.rm = TRUE), 3),
      .groups = "drop"
    ) |>
    arrange(as.integer(cluster)) |>
    print()
}
