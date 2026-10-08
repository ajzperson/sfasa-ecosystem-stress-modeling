library(tidyverse)
library(sf)
library(patchwork)
library(cluster)

source("pixel_approach/src/clustering_helpers.R")

# K-means clustering on density-residualized tau_r1 + tau_skew, run
# independently for NDVI and NDWI (not combined into one 4-feature
# clustering) so each keeps its own cluster assignment and the two can be
# cross-validated against each other. tau_r1 and tau_skew both correlate
# with a pixel's own baseline vegetation density (denser forest has more
# fuel to burn, independent of actual CSD dynamics), so both are
# residualized against it before clustering.
#
# The ANOVA/Tukey/ROC validation below (against each index's own
# actual-drop ground truth) is commented out -- MTBS fire + USFS ADS
# beetle mortality (checked separately in 09/08) is the current
# validation instead.
#
# Outputs:
#   pixel_approach/data/processed/ndvi_pixel_clusters.csv
#   pixel_approach/data/processed/ndwi_pixel_clusters.csv
#   pixel_approach/data/processed/cluster_profile_ndvi.csv (per-cluster mean tau/resid values)
#   pixel_approach/data/processed/cluster_profile_ndwi.csv
#   pixel_approach/charts/cluster_k_sweep.pdf         (NDVI/NDWI stacked)
#   pixel_approach/charts/cluster_map.pdf              (NDVI/NDWI side by side)

if (!exists("K"))       K       <- 3
if (!exists("SEED"))    SEED    <- 42
if (!exists("INDICES")) INDICES <- c("ndvi", "ndwi")

dir.create("pixel_approach/charts", showWarnings = FALSE, recursive = TRUE)

if (!exists("NDVI_TAU_FILE")) NDVI_TAU_FILE <- "pixel_approach/data/processed/ndvi_pixel_tau.csv"
if (!exists("NDWI_TAU_FILE")) NDWI_TAU_FILE <- "pixel_approach/data/processed/ndwi_pixel_tau.csv"

index_config <- list(
  ndvi = list(
    tau_path     = NDVI_TAU_FILE,
    mat_path     = "pixel_approach/data/processed/ndvi_matrix.rds",
    # drop_col   = "Actual_Drop_15_16",  # deprioritized (see header)
    clusters_out = "pixel_approach/data/processed/ndvi_pixel_clusters.csv",
    label        = "NDVI"
  ),
  ndwi = list(
    tau_path     = NDWI_TAU_FILE,
    mat_path     = "pixel_approach/data/processed/ndwi_matrix.rds",
    # drop_col   = "Actual_NDWI_Drop_15_16",  # deprioritized (see header)
    clusters_out = "pixel_approach/data/processed/ndwi_pixel_clusters.csv",
    label        = "NDWI"
  )
)

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

# MTBS-based fire flags (02c_pixel_fire.R) -- shared across both indices,
# since fire exposure doesn't depend on which vegetation index is used.
fire_flags <- read_csv("pixel_approach/data/processed/pixel_fire.csv", show_col_types = FALSE) |>
  select(pixel_id, burned_baseline_00_14)

# Shared pixel_id key for joining baseline-density residualization onto
# either index's matrix, which has no pixel_id column of its own -- row i
# of ndvi_matrix.rds / ndwi_matrix.rds corresponds to pixel_meta$pixel_id[i].
pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv", show_col_types = FALSE) |>
  select(pixel_id)

# Colors: blue (rank 1, least stressed) -> red (rank K, most stressed).
# Shared across both indices' charts since the cluster labels (1..K,
# "Least-stressed".."Most-stressed") mean the same thing in both.
style        <- build_cluster_style(K)
cluster_cols <- style$cols
clust_labels <- style$labels

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

  # Exclude pixels burned (MTBS moderate+ severity) at any point during the
  # 2000-2014 baseline window -- a fire causes a sudden, discrete NDVI/NDWI
  # crash unrelated to gradual drought stress, which would contaminate the
  # tau_r1/tau_skew trend statistics computed over that same window.
  n_before <- nrow(d)
  d <- d |>
    left_join(fire_flags, by = "pixel_id") |>
    filter(!burned_baseline_00_14) |>
    select(-burned_baseline_00_14)
  cat(sprintf("✓ Excluded %d pixels burned (MTBS moderate+) during 2000-2014 baseline (%d remain)\n",
              n_before - nrow(d), nrow(d)))

  # Baseline vegetation density (mean pre-2015 level) -- tau_r1 and
  # tau_skew both correlate with how much forest a pixel has to begin
  # with, which otherwise leaks into cluster assignment (see header).
  base_level <- compute_base_level(cfg$mat_path, pixel_meta$pixel_id, "base_level")
  d <- d |> left_join(base_level, by = "pixel_id")
  stopifnot(!anyNA(d$base_level))

  d <- residualize_tau(d, "tau_r1", "tau_skew", "base_level",
                       "resid_tau_r1", "resid_tau_skew")

  # Scale features before clustering (Euclidean distance is scale-sensitive).
  # resid_tau_skew is used unflipped -- its sign isn't forced to match
  # resid_tau_r1's "higher = more stressed" convention, since skewness's
  # direction near a transition depends on which state the system is
  # heading toward (not universal, unlike r1). Left raw so k-means can
  # separate pixels by that direction rather than having it collapsed away.
  feat <- d |> transmute(resid_tau_r1, resid_tau_skew) |> scale()

  # --- K sweep - elbow + silhouette (K = 2 to 8), then K-means at chosen K ---
  cat("Running K sweep (K=2 to 8) for elbow and silhouette plots...\n")
  sweep_res <- run_kmeans_sweep(feat, K, SEED)
  sweep_df  <- sweep_res$sweep_df

  cat(sprintf("Running k-means (k=%d, seed=%d)...\n", K, SEED))
  km <- sweep_res$km_final
  d$cluster_raw <- km$cluster

  # Re-label clusters by composite stress rank (red = more stressed,
  # blue = less stressed).
  stress_rank <- d |>
    group_by(cluster_raw) |>
    summarise(composite = mean(resid_tau_r1) + mean(resid_tau_skew), .groups = "drop") |>
    arrange(composite) |>
    mutate(stress_rank = row_number())   # 1 = least stressed, K = most stressed

  d <- d |>
    left_join(stress_rank |> select(cluster_raw, stress_rank),
              by = "cluster_raw") |>
    mutate(cluster = factor(stress_rank, levels = seq_len(K)))

  cat(sprintf("✓ Cluster sizes (stress rank 1=least → %d=most): %s\n",
              K, paste(table(d$cluster), collapse = " / ")))

  # Print cluster centroids in original (unscaled) units
  # drop_col <- cfg$drop_col  # deprioritized (see header)
  centers_orig <- d |>
    group_by(cluster) |>
    summarise(
      n           = n(),
      tau_r1        = round(mean(tau_r1),   4),
      tau_var       = round(mean(tau_var),  4),
      tau_skew      = round(mean(tau_skew), 4),
      resid_tau_r1  = round(mean(resid_tau_r1),   4),
      resid_tau_skew = round(mean(resid_tau_skew), 4),
      composite     = round(mean(resid_tau_r1) + mean(resid_tau_skew), 4),
      # actual_drop = round(mean(.data[[drop_col]], na.rm = TRUE), 4),  # deprioritized
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
# returned unchanged when only one is, so every chart below genuinely
# follows INDICES instead of assuming both NDVI and NDWI always ran.
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
  build_k_sweep_panel(res$sweep_df, K, res$cfg$label)
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
  build_cluster_map_panel(res$d, sn_boundary, cluster_cols, clust_labels, res$cfg$label)
}

map_panels <- setNames(lapply(INDICES, function(idx) build_map_panel(results[[idx]])), INDICES)

# The caption is generic boilerplate (not panel-specific), but
# build_cluster_map_panel puts it on every panel -- with panels side by
# side, two copies of the same text sit close enough to visually collide.
# Blanked on every panel here and re-added once below via plot_annotation
# instead, so it's centered under the full combined width (shared by both
# panels) rather than attached to just one side.
map_panels <- lapply(map_panels, function(pnl) pnl + labs(caption = NULL))

p_map <- combine_panels(map_panels) +
  plot_annotation(
    title    = sprintf("K-means Clusters (k=%d) on density-residualized tau_r1 + tau_skew", K),
    subtitle = "Color: blue = least stressed -> red = most stressed",
    caption  = "White gaps: pixels burned (MTBS moderate+ severity) at any point\n2000-2014, excluded before clustering",
    theme    = theme(plot.title    = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"),
                      plot.subtitle = element_text(hjust = 0.5, color = "#1B2A1F", size = 9, face = "bold"),
                      plot.caption  = element_text(hjust = 0.5, color = "#1B2A1F", size = 8, face = "bold",
                                                    margin = margin(t = 10)))
  )

ggsave("pixel_approach/charts/cluster_map.pdf",
       plot = p_map, width = 7 * n_idx, height = 7, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_map.pdf\n")

# Individual per-index map panels are NOT saved as .rds here -- a ggplot
# object embeds its full underlying data, and this map plots ~220k raw
# points, making each serialized panel ~175MB. 11_dashboard.R instead
# rebuilds this panel directly from ndvi_pixel_clusters.csv/
# ndwi_pixel_clusters.csv (cheap, and that data already persists).

# # ----------------------------------------------------------------------------
# # ANOVA + Tukey HSD on each index's own actual-drop column by cluster.
# # Spatial thinning (1 pixel per 0.04° x 0.04° cell, ~4km spacing) avoids
# # inflating the F-statistic via autocorrelation between neighboring 500m
# # pixels. Boxplot shows all pixels; ANOVA uses the thinned subset.
# # ----------------------------------------------------------------------------
# run_validation <- function(res, SEED) {
#   d        <- res$d
#   drop_col <- res$cfg$drop_col
#   val_data <- d |> filter(!is.na(.data[[drop_col]]))
# 
#   set.seed(SEED)
#   val_thin <- val_data |>
#     mutate(lon_bin = round(lon / 0.04) * 0.04,
#            lat_bin = round(lat / 0.04) * 0.04) |>
#     group_by(lon_bin, lat_bin) |>
#     slice_sample(n = 1) |>
#     ungroup()
# 
#   cat(sprintf("\n=== %s: spatial thinning: %d → %d pixels (%.0f%%, >=4km spacing) ===\n",
#               res$cfg$label, nrow(val_data), nrow(val_thin),
#               100 * nrow(val_thin) / nrow(val_data)))
# 
#   # bquote()/.() substitutes the real column name in as a symbol (not just
#   # a string that evaluates to one), so aov_fit$call — and therefore
#   # TukeyHSD's printed "Fit:" line — shows a clean `Actual_Drop_15_16 ~
#   # cluster` instead of a literal val_thin[[drop_col]] indexing expression
#   # or an as.formula(paste(...)) call. as.formula(paste(drop_col, "~
#   # cluster")) fits correctly too but still prints the paste() expression
#   # itself, since R's stored call reflects the calling code, not its value.
#   aov_fit <- eval(bquote(aov(.(as.name(drop_col)) ~ cluster, data = val_thin)))
#   cat(sprintf("=== One-way ANOVA (spatially thinned): %s ~ cluster ===\n", drop_col))
#   print(summary(aov_fit))
# 
#   tukey_result <- TukeyHSD(aov_fit)
#   cat("\n=== Tukey HSD pairwise comparisons ===\n")
#   print(tukey_result)
# 
#   mean_df <- val_data |>
#     group_by(cluster) |>
#     summarise(mean_val = mean(.data[[drop_col]], na.rm = TRUE), .groups = "drop") |>
#     mutate(label = sprintf("mean = %.3f", mean_val))
# 
#   format_p <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
#   tk <- tukey_result$cluster
#   # geom_signif comparisons order (1,2)/(2,3)/(1,3) vs. TukeyHSD row order (2-1)/(3-1)/(3-2)
#   sig_labels <- c(format_p(tk["2-1", "p adj"]),
#                    format_p(tk["3-2", "p adj"]),
#                    format_p(tk["3-1", "p adj"]))
#   overall_p <- summary(aov_fit)[[1]][["Pr(>F)"]][1]
# 
#   p_val <- ggplot(val_data, aes(x = cluster, y = .data[[drop_col]], fill = cluster)) +
#     geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, width = 0.5) +
#     geom_text(data = mean_df, aes(x = cluster, y = mean_val, label = label),
#               inherit.aes = FALSE, size = 3.2, vjust = -0.6, fontface = "italic") +
#     geom_signif(
#       comparisons  = list(c("1", "2"), c("2", "3"), c("1", "3")),
#       annotations  = sig_labels,
#       step_increase = 0.12,
#       tip_length   = 0.01,
#       textsize     = 3.2
#     ) +
#     scale_fill_manual(values = cluster_cols, guide = "none") +
#     labs(
#       title    = res$cfg$label,
#       subtitle = sprintf("ANOVA F = %.1f, %s (n = %s, >=4km spacing) - Tukey HSD shown",
#                          summary(aov_fit)[[1]][["F value"]][1],
#                          format_p(overall_p),
#                          format(nrow(val_thin), big.mark = ",")),
#       x = "Cluster", y = sprintf("%s (baseline mean - 2015-2016 min)", drop_col)
#     ) +
#     theme_bw() +
#     theme(
#       plot.title    = element_text(hjust = 0.5, face = "bold"),
#       plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9)
#     )
# 
#   list(plot = p_val, aov_fit = aov_fit, val_thin = val_thin)
# }
# 
# val_results <- setNames(
#   lapply(INDICES, function(idx) run_validation(results[[idx]], SEED)),
#   INDICES
# )
# 
# p_val <- combine_panels(lapply(INDICES, function(idx) val_results[[idx]]$plot)) +
#   plot_annotation(
#     title = "Actual Drop 2015-2016 by Cluster",
#     theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#   )
# 
# ggsave("pixel_approach/charts/cluster_validation.pdf",
#        plot = p_val, width = 7 * n_idx, height = 6, dpi = 150)
# cat("✓ Saved: pixel_approach/charts/cluster_validation.pdf\n")
# 
# # ----------------------------------------------------------------------------
# # Binary classification validation - confusion matrix + ROC-AUC on each
# # index's own actual-drop ground truth. ROC axes (TPR/FPR) are already
# # unitless [0,1] regardless of source index, so NDVI and NDWI are overlaid
# # on one plot rather than shown as separate panels.
# # ----------------------------------------------------------------------------
# compute_roc <- function(res, K) {
#   drop_col <- res$cfg$drop_col
#   clf_data <- res$d |>
#     filter(!is.na(.data[[drop_col]])) |>
#     mutate(
#       # Ground truth: top 20% of actual drop = "severe collapse"
#       severe    = .data[[drop_col]] >= quantile(.data[[drop_col]], 0.80, na.rm = TRUE),
#       # Prediction: cluster K (full CSD) = "predicted severe"
#       predicted = as.integer(cluster) == K
#     )
# 
#   TP <- sum( clf_data$severe &  clf_data$predicted)
#   FP <- sum(!clf_data$severe &  clf_data$predicted)
#   TN <- sum(!clf_data$severe & !clf_data$predicted)
#   FN <- sum( clf_data$severe & !clf_data$predicted)
# 
#   sensitivity <- TP / (TP + FN)
#   precision   <- TP / (TP + FP)
#   specificity <- TN / (TN + FP)
#   fpr_binary  <- FP / (FP + TN)
#   accuracy    <- (TP + TN) / (TP + FP + TN + FN)
#   lift        <- precision / 0.20   # baseline rate = 20% by definition
# 
#   cat(sprintf("\n=== %s: Binary classification, Cluster %d vs top-20%% %s ===\n",
#               res$cfg$label, K, drop_col))
#   cat(sprintf("  Threshold: top 20%% (>= %.4f)\n",
#               quantile(clf_data[[drop_col]], 0.80, na.rm = TRUE)))
#   cat(sprintf("  TP=%d  FP=%d  TN=%d  FN=%d\n", TP, FP, TN, FN))
#   cat(sprintf("  Sensitivity (Recall): %.3f\n", sensitivity))
#   cat(sprintf("  Precision:            %.3f\n", precision))
#   cat(sprintf("  Specificity:          %.3f\n", specificity))
#   cat(sprintf("  False Positive Rate:  %.3f\n", fpr_binary))
#   cat(sprintf("  Accuracy:             %.3f\n", accuracy))
#   cat(sprintf("  Lift factor:          %.2fx  (vs 1.00x random)\n", lift))
# 
#   # ROC-AUC via trapezoidal rule on composite stress score (same
#   # density-residualized tau_r1 - tau_skew composite used for clustering)
#   roc_df <- clf_data |>
#     mutate(stress_score = resid_tau_r1 - resid_tau_skew) |>
#     select(stress_score, severe) |>
#     arrange(desc(stress_score)) |>
#     mutate(
#       tpr   = cumsum(severe)  / sum(severe),
#       fpr   = cumsum(!severe) / sum(!severe),
#       index = res$cfg$label
#     )
# 
#   auc <- sum(diff(roc_df$fpr) * (head(roc_df$tpr, -1) + tail(roc_df$tpr, -1))) / 2
#   cat(sprintf("  ROC-AUC:              %.3f\n", auc))
# 
#   list(roc_df = roc_df, auc = auc)
# }
# 
# roc_results <- setNames(
#   lapply(INDICES, function(idx) compute_roc(results[[idx]], K)),
#   INDICES
# )
# 
# roc_combined <- bind_rows(lapply(roc_results, function(r) r$roc_df))
# roc_colors   <- c(NDVI = "#D73027", NDWI = "#2166AC")
# auc_label    <- paste(sprintf("%s AUC = %.3f",
#                               sapply(INDICES, function(idx) index_config[[idx]]$label),
#                               sapply(roc_results, function(r) r$auc)),
#                       collapse = "\n")
# 
# p_roc <- ggplot(roc_combined, aes(x = fpr, y = tpr, color = index)) +
#   geom_line(linewidth = 0.8) +
#   geom_abline(linetype = "dashed", color = "grey50") +
#   scale_color_manual(values = roc_colors, name = NULL) +
#   annotate("text", x = 0.55, y = 0.15,
#            label = auc_label,
#            size = 4, fontface = "bold", hjust = 0) +
#   labs(
#     title    = "ROC Curve: CSD Stress Score vs Severe Canopy Collapse",
#     subtitle = "Ground truth = top 20% actual drop (own index); predictor = density-residualized tau_r1 - tau_skew",
#     x = "False Positive Rate", y = "True Positive Rate (Sensitivity)"
#   ) +
#   theme_bw() +
#   theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
#         plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))
# 
# ggsave("pixel_approach/charts/cluster_roc.pdf",
#        plot = p_roc, width = 6.5, height = 5.5, dpi = 150)
# cat("✓ Saved: pixel_approach/charts/cluster_roc.pdf\n")

# ----------------------------------------------------------------------------
# Summary table, per index
# ----------------------------------------------------------------------------
for (idx in INDICES) {
  res <- results[[idx]]
  cat(sprintf("\n=== %s: final cluster profile (1=least stressed, K=most stressed) ===\n",
              res$cfg$label))
  profile <- res$d |>
    group_by(cluster) |>
    summarise(
      n                   = n(),
      tau_r1_mean         = round(mean(tau_r1),   3),
      tau_var_mean        = round(mean(tau_var),  3),
      tau_skew_mean       = round(mean(tau_skew), 3),
      resid_tau_r1_mean   = round(mean(resid_tau_r1),   3),
      resid_tau_skew_mean = round(mean(resid_tau_skew), 3),
      .groups = "drop"
    ) |>
    arrange(as.integer(cluster))
  print(profile)

  # Persisted to disk (not just console output) so these numbers are
  # inspectable later without re-running the script.
  profile_path <- sprintf("pixel_approach/data/processed/cluster_profile_%s.csv", idx)
  write_csv(profile, profile_path)
  cat(sprintf("✓ Saved: %s\n", profile_path))
}
