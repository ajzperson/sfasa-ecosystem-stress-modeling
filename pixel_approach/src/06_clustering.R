library(tidyverse)
library(sf)
library(patchwork)
library(multcompView)
library(ggsignif)
library(cluster)

# K-means clustering on tau_r1 + tau_var + tau_skew (all three CSD signals).
# tau_skew included after VIF confirmed no multicollinearity (all VIF ~1.01).
# Actual_Drop_15_16 is NOT a clustering feature — post-hoc validation only.
#
# Outputs:
#   pixel_approach/data/processed/pixel_clusters.csv
#   pixel_approach/charts/cluster_map.pdf
#   pixel_approach/charts/cluster_ndvi_timelines.pdf
#   pixel_approach/charts/cluster_validation.pdf  (Actual_Drop by cluster)

if (!exists("K"))    K    <- 3
if (!exists("SEED")) SEED <- 42

dir.create("pixel_approach/charts", showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------
# 1. Load and prepare
# --------------------------------------------------------------------------
cat("Loading pixel_tau.csv...\n")
d <- read_csv("pixel_approach/data/processed/pixel_tau.csv",
              show_col_types = FALSE) |>
  filter(!is.na(tau_r1), !is.na(tau_var), !is.na(tau_skew))

cat(sprintf("✓ %d pixels with valid tau_r1, tau_var, tau_skew\n", nrow(d)))

# Scale features before clustering (Euclidean distance is scale-sensitive)
feat <- d |> select(tau_r1, tau_var, tau_skew) |> scale()

# --------------------------------------------------------------------------
# 2. K sweep — elbow + silhouette (K = 2 to 8)
# --------------------------------------------------------------------------
cat("Running K sweep (K=2 to 8) for elbow and silhouette plots...\n")

k_range  <- 2:8
wcss_vec <- numeric(length(k_range))
sil_vec  <- numeric(length(k_range))

# Subsample for silhouette (full distance matrix on 248k pixels is too large)
set.seed(SEED)
sil_idx  <- sample(nrow(feat), min(5000, nrow(feat)))
feat_sil <- feat[sil_idx, ]

for (ki in seq_along(k_range)) {
  k_i <- k_range[ki]
  set.seed(SEED)
  km_i <- kmeans(feat, centers = k_i, nstart = 25, iter.max = 100)
  wcss_vec[ki] <- km_i$tot.withinss
  sil_vec[ki]  <- mean(silhouette(km_i$cluster[sil_idx],
                                  dist(feat_sil))[, 3])
  cat(sprintf("  k=%d  WCSS=%.1f  Silhouette=%.3f\n",
              k_i, wcss_vec[ki], sil_vec[ki]))
}

sweep_df <- data.frame(k = k_range, wcss = wcss_vec, silhouette = sil_vec)

p_elbow <- ggplot(sweep_df, aes(x = k, y = wcss)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_point(data = sweep_df[sweep_df$k == K, ],
             size = 5, shape = 21, fill = "#D73027", color = "black") +
  scale_x_continuous(breaks = k_range) +
  labs(title = "Elbow Method",
       x = "Number of clusters (K)", y = "Total within-cluster SS") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_sil <- ggplot(sweep_df, aes(x = k, y = silhouette)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_point(data = sweep_df[sweep_df$k == K, ],
             size = 5, shape = 21, fill = "#D73027", color = "black") +
  scale_x_continuous(breaks = k_range) +
  labs(title = "Silhouette Score",
       x = "Number of clusters (K)", y = "Mean silhouette width") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p_sweep <- (p_elbow + p_sil) +
  plot_annotation(
    title    = sprintf("K-means Hyperparameter Sweep (K = 2–8); selected K = %d shown in red", K),
    theme    = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

ggsave("pixel_approach/charts/cluster_k_sweep.pdf",
       plot = p_sweep, width = 10, height = 4, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_k_sweep.pdf\n")

# --------------------------------------------------------------------------
# 3. K-means
# --------------------------------------------------------------------------
cat(sprintf("Running k-means (k=%d, seed=%d)...\n", K, SEED))
set.seed(SEED)
km <- kmeans(feat, centers = K, nstart = 25, iter.max = 100)

d$cluster_raw <- km$cluster

# Re-label clusters by composite stress rank so colors are semantically
# consistent with the heatmaps (red = more stressed, blue = less stressed).
# Composite = tau_r1 + tau_var - tau_skew (tau_skew negative = worse, so negate).
stress_rank <- d |>
  group_by(cluster_raw) |>
  summarise(composite = mean(tau_r1) + mean(tau_var) - mean(tau_skew),
            .groups = "drop") |>
  arrange(composite) |>
  mutate(stress_rank = row_number())   # 1 = least stressed, K = most stressed

d <- d |>
  left_join(stress_rank |> select(cluster_raw, stress_rank),
            by = "cluster_raw") |>
  mutate(cluster = factor(stress_rank, levels = seq_len(K)))

cat(sprintf("✓ Cluster sizes (stress rank 1=least → %d=most): %s\n",
            K, paste(table(d$cluster), collapse = " / ")))

# Print cluster centroids in original (unscaled) units
centers_orig <- d |>
  group_by(cluster) |>
  summarise(
    n           = n(),
    tau_r1      = round(mean(tau_r1),   4),
    tau_var     = round(mean(tau_var),  4),
    tau_skew    = round(mean(tau_skew), 4),
    composite   = round(mean(tau_r1) + mean(tau_var) - mean(tau_skew), 4),
    Actual_Drop = round(mean(Actual_Drop_15_16, na.rm = TRUE), 4),
    .groups = "drop"
  ) |>
  arrange(as.integer(cluster))

cat("\nCluster profile (1=least stressed, K=most stressed):\n")
print(centers_orig)

# --------------------------------------------------------------------------
# 4. Save clustered pixel table
# --------------------------------------------------------------------------
out_csv <- "pixel_approach/data/processed/pixel_clusters.csv"
write_csv(d, out_csv)
cat(sprintf("\n✓ Saved: %s\n", out_csv))

# --------------------------------------------------------------------------
# 5. Cluster map
# --------------------------------------------------------------------------
sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

# Colors match heatmap convention: blue (rank 1, least stressed) → red (rank K, most stressed).
# Generated dynamically so rank K is always red regardless of K.
cluster_cols <- setNames(
  colorRampPalette(c("#2166AC", "#92C5DE", "#F4A582", "#D73027"))(K),
  as.character(seq_len(K))
)

p_map <- ggplot(d, aes(x = lon, y = lat, color = cluster)) +
  geom_point(size = 0.05, alpha = 0.5) +
  geom_sf(data = sn_boundary, fill = NA, color = "grey20",
          linewidth = 0.4, inherit.aes = FALSE) +
  scale_color_manual(
    values = cluster_cols,
    name   = "Stress level",
    labels = paste0(seq_len(K), c(" — no CSD", rep(" — partial", K - 2), " — full CSD"))
  ) +
  coord_sf(xlim = range(d$lon), ylim = range(d$lat)) +
  labs(
    title    = sprintf("K-means Clusters (k=%d) on τ_r1 + τ_var + τ_skew", K),
    subtitle = "Color: blue = least stressed → red = most stressed (consistent with heatmaps)",
    x = "Longitude", y = "Latitude"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    axis.text       = element_text(size = 7)
  )

ggsave("pixel_approach/charts/cluster_map.pdf",
       plot = p_map, width = 8, height = 7, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_map.pdf\n")

# --------------------------------------------------------------------------
# 6. Pairwise scatter: all three 2D projections of the 3D feature space
# --------------------------------------------------------------------------
samp <- d |> slice_sample(n = min(30000, nrow(d)))
clust_labels <- paste0(seq_len(K), c(" — no CSD", rep(" — partial", K - 2),
                                     " — full CSD"))

make_scatter <- function(xvar, yvar, xlabel, ylabel) {
  ggplot(samp, aes(x = .data[[xvar]], y = .data[[yvar]], color = cluster)) +
    geom_point(size = 0.3, alpha = 0.4) +
    scale_color_manual(values = cluster_cols, name = "Stress level",
                       labels = clust_labels) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = xlabel, y = ylabel) +
    theme_bw()
}

p1 <- make_scatter("tau_r1",  "tau_var",  "τ r₁ (autocorrelation)", "τ variance")
p2 <- make_scatter("tau_r1",  "tau_skew", "τ r₁ (autocorrelation)", "τ skewness")
p3 <- make_scatter("tau_var", "tau_skew", "τ variance",              "τ skewness")

p_scatter <- (p1 + p2 + p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Cluster Separation — Pairwise CSD Feature Projections",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  )

ggsave("pixel_approach/charts/cluster_scatter_pairwise.pdf",
       plot = p_scatter, width = 14, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_scatter_pairwise.pdf\n")

# --------------------------------------------------------------------------
# 7. Per-cluster NDVI timelines
# --------------------------------------------------------------------------
cat("Loading NDVI matrix for per-cluster timelines...\n")
ndvi_mat  <- readRDS("pixel_approach/data/processed/ndvi_matrix.rds")

# Align pixel_clusters rows to ndvi_matrix rows via pixel_id
# pixel_id is 1-indexed row position in the original ndvi_matrix
pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv",
                       show_col_types = FALSE)

# d$pixel_id indexes into ndvi_mat rows
col_dates <- as.Date(colnames(ndvi_mat))

# Aggregate to annual means to remove seasonal oscillations
timeline_list <- vector("list", K)
for (k in seq_len(K)) {
  idx <- d$pixel_id[d$cluster == k]
  sub <- ndvi_mat[idx, , drop = FALSE]
  monthly <- data.frame(date = col_dates,
                        ndvi = colMeans(sub, na.rm = TRUE),
                        cluster = factor(k))
  timeline_list[[k]] <- monthly |>
    mutate(year  = as.integer(format(date, "%Y")),
           month = as.integer(format(date, "%m"))) |>
    filter(month %in% 3:11) |>       # growing season only (Mar–Nov)
    group_by(cluster, year) |>
    summarise(ndvi = mean(ndvi, na.rm = TRUE), .groups = "drop")
}
timelines <- bind_rows(timeline_list)

p_ts <- ggplot(timelines, aes(x = year, y = ndvi, color = cluster)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = cluster_cols, name = "Stress level",
    labels = paste0(seq_len(K), c(" — no CSD", rep(" — partial", K - 2),
                                  " — full CSD"))) +
  geom_vline(xintercept = 2012, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 2015, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = 2012.2, y = Inf, label = "Drought onset",
           hjust = 0, vjust = 1.5, size = 2.8, color = "grey40") +
  annotate("text", x = 2015.2, y = Inf, label = "Peak mortality",
           hjust = 0, vjust = 1.5, size = 2.8, color = "grey40") +
  labs(
    title = "Annual Mean NDVI Trajectory by Cluster (2000–2016)",
    x = NULL, y = "Annual Mean NDVI"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("pixel_approach/charts/cluster_ndvi_timelines.pdf",
       plot = p_ts, width = 9, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_ndvi_timelines.pdf\n")

# --------------------------------------------------------------------------
# 8. ANOVA + Tukey HSD — run first so letters can go on the boxplot
# --------------------------------------------------------------------------
val_data <- d |> filter(!is.na(Actual_Drop_15_16))

# Spatial thinning: sample one pixel per 0.04° × 0.04° cell (~4km spacing)
# to avoid inflating F-statistic due to spatial autocorrelation between
# neighboring 500m pixels. Boxplot shows all pixels; ANOVA uses thinned subset.
set.seed(SEED)
val_thin <- val_data |>
  mutate(lon_bin = round(lon / 0.04) * 0.04,
         lat_bin = round(lat / 0.04) * 0.04) |>
  group_by(lon_bin, lat_bin) |>
  slice_sample(n = 1) |>
  ungroup()

cat(sprintf("\n=== Spatial thinning: %d → %d pixels (%.0f%%, ≥4km spacing) ===\n",
            nrow(val_data), nrow(val_thin),
            100 * nrow(val_thin) / nrow(val_data)))

aov_fit  <- aov(Actual_Drop_15_16 ~ cluster, data = val_thin)

cat("=== One-way ANOVA (spatially thinned): Actual_Drop_15_16 ~ cluster ===\n")
print(summary(aov_fit))

tukey_result <- TukeyHSD(aov_fit)
cat("\n=== Tukey HSD pairwise comparisons ===\n")
print(tukey_result)

# Mean value label per cluster (printed inside each box)
mean_df <- val_data |>
  group_by(cluster) |>
  summarise(mean_drop = mean(Actual_Drop_15_16, na.rm = TRUE), .groups = "drop") |>
  mutate(label = sprintf("mean = %.3f", mean_drop))

# --------------------------------------------------------------------------
# 9. Validation boxplot with mean labels + significance brackets
# --------------------------------------------------------------------------
p_val <- ggplot(val_data, aes(x = cluster, y = Actual_Drop_15_16, fill = cluster)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, width = 0.5) +
  geom_text(data = mean_df, aes(x = cluster, y = mean_drop, label = label),
            inherit.aes = FALSE, size = 3.2, vjust = -0.6, fontface = "italic") +
  geom_signif(
    comparisons  = list(c("1", "2"), c("2", "3"), c("1", "3")),
    annotations  = c("p < 0.001", "p < 0.001", "p < 0.001"),
    step_increase = 0.12,
    tip_length   = 0.01,
    textsize     = 3.2
  ) +
  scale_fill_manual(values = cluster_cols, guide = "none") +
  labs(
    title    = "Post-hoc Validation: Actual NDVI Drop 2015–2016 by Cluster",
    subtitle = sprintf("ANOVA F = %.0f, p < 0.001 (spatially thinned n = %s, ≥4km spacing) — Tukey HSD shown",
                       summary(aov_fit)[[1]][["F value"]][1],
                       format(nrow(val_thin), big.mark = ",")),
    x = "Cluster", y = "Actual Drop (baseline mean − 2015–2016 min)"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9)
  )

ggsave("pixel_approach/charts/cluster_validation.pdf",
       plot = p_val, width = 7, height = 6, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_validation.pdf\n")

# --------------------------------------------------------------------------
# 10. Binary classification validation — confusion matrix + ROC-AUC
# --------------------------------------------------------------------------
# Reload from saved CSV so Actual_Drop reflects the latest 04_actual_drop.R run
clf_data <- read_csv("pixel_approach/data/processed/pixel_clusters.csv",
                     show_col_types = FALSE) |>
  filter(!is.na(Actual_Drop_15_16)) |>
  mutate(
    # Ground truth: top 20% of actual drop = "severe collapse"
    actual_severe    = Actual_Drop_15_16 >= quantile(Actual_Drop_15_16, 0.80, na.rm = TRUE),
    # Prediction: cluster K (full CSD) = "predicted severe"
    predicted_severe = as.integer(cluster) == K,
    # Continuous score for ROC: composite stress (higher = more stressed)
    stress_score     = tau_r1 + tau_var - tau_skew
  )

TP <- sum( clf_data$actual_severe &  clf_data$predicted_severe)
FP <- sum(!clf_data$actual_severe &  clf_data$predicted_severe)
TN <- sum(!clf_data$actual_severe & !clf_data$predicted_severe)
FN <- sum( clf_data$actual_severe & !clf_data$predicted_severe)

sensitivity <- TP / (TP + FN)
precision   <- TP / (TP + FP)
specificity <- TN / (TN + FP)
fpr_binary  <- FP / (FP + TN)
accuracy    <- (TP + TN) / (TP + FP + TN + FN)

cat("\n=== Binary classification: Cluster 3 vs top-20% Actual_Drop ===\n")
cat(sprintf("  Threshold: top %.0f%% of Actual_Drop (≥ %.4f)\n",
            20, quantile(clf_data$Actual_Drop_15_16, 0.80, na.rm = TRUE)))
cat(sprintf("  TP=%d  FP=%d  TN=%d  FN=%d\n", TP, FP, TN, FN))
lift <- precision / 0.20   # baseline rate = 20% by definition
cat(sprintf("  Sensitivity (Recall): %.3f\n", sensitivity))
cat(sprintf("  Precision:            %.3f\n", precision))
cat(sprintf("  Specificity:          %.3f\n", specificity))
cat(sprintf("  False Positive Rate:  %.3f\n", fpr_binary))
cat(sprintf("  Accuracy:             %.3f\n", accuracy))
cat(sprintf("  Lift factor:          %.2fx  (vs 1.00x random)\n", lift))

# ROC-AUC via trapezoidal rule on composite stress score
roc_df <- clf_data |>
  select(stress_score, actual_severe) |>
  arrange(desc(stress_score)) |>
  mutate(
    tpr = cumsum(actual_severe)  / sum(actual_severe),
    fpr = cumsum(!actual_severe) / sum(!actual_severe)
  )

auc <- sum(diff(roc_df$fpr) * (head(roc_df$tpr, -1) + tail(roc_df$tpr, -1))) / 2
cat(sprintf("  ROC-AUC:              %.3f\n", auc))

# ROC curve plot
p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr)) +
  geom_line(linewidth = 0.8, color = "#D73027") +
  geom_abline(linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.6, y = 0.2,
           label = sprintf("AUC = %.3f", auc),
           size = 4.5, fontface = "bold") +
  labs(
    title    = "ROC Curve: CSD Stress Score vs Severe Canopy Collapse",
    subtitle = "Ground truth = top 20% Actual_Drop; predictor = τ_r1 + τ_var − τ_skew",
    x = "False Positive Rate", y = "True Positive Rate (Sensitivity)"
  ) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))

ggsave("pixel_approach/charts/cluster_roc.pdf",
       plot = p_roc, width = 6, height = 5, dpi = 150)
cat("✓ Saved: pixel_approach/charts/cluster_roc.pdf\n")

# --------------------------------------------------------------------------
# 11. Summary table
# --------------------------------------------------------------------------
cat("\n=== Final cluster profile (1=least stressed, K=most stressed) ===\n")
d |>
  group_by(cluster) |>
  summarise(
    n             = n(),
    tau_r1_mean   = round(mean(tau_r1),   3),
    tau_var_mean  = round(mean(tau_var),  3),
    tau_skew_mean = round(mean(tau_skew), 3),
    drop_mean     = round(mean(Actual_Drop_15_16,   na.rm = TRUE), 3),
    drop_median   = round(median(Actual_Drop_15_16, na.rm = TRUE), 3),
    .groups = "drop"
  ) |>
  arrange(as.integer(cluster)) |>
  print()
