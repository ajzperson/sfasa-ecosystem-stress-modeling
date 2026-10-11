library(tidyverse)
library(sf)
library(patchwork)
library(cluster)

source("pixel_approach/src/clustering_helpers.R")

# K-means stress tiers from density-residualized tau_r1 + tau_skew, run
# separately for NDVI and NDWI.
#
# Output: pixel_approach/data/processed/{ndvi,ndwi}_pixel_clusters.csv
#         pixel_approach/data/processed/cluster_profile_{ndvi,ndwi}.csv
#         pixel_approach/charts/cluster_k_sweep.pdf
#         pixel_approach/charts/cluster_map.pdf

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
    clusters_out = "pixel_approach/data/processed/ndvi_pixel_clusters.csv",
    label        = "NDVI"
  ),
  ndwi = list(
    tau_path     = NDWI_TAU_FILE,
    mat_path     = "pixel_approach/data/processed/ndwi_matrix.rds",
    clusters_out = "pixel_approach/data/processed/ndwi_pixel_clusters.csv",
    label        = "NDWI"
  )
)

sn_boundary <- st_read("data_common/processed/sierra_nevada_boundary.geojson", quiet = TRUE)

fire_flags <- read_csv("pixel_approach/data/processed/pixel_fire.csv", show_col_types = FALSE) |>
  select(pixel_id, burned_baseline_00_14)

# Row i of each index matrix is pixel_meta$pixel_id[i]
pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv", show_col_types = FALSE) |>
  select(pixel_id)

style        <- build_cluster_style(K)
cluster_cols <- style$cols
clust_labels <- style$labels

run_clustering <- function(cfg, K, SEED) {
  cat(sprintf("\n=== %s ===\n", cfg$label))
  cat(sprintf("Loading %s...\n", cfg$tau_path))
  d <- read_csv(cfg$tau_path, show_col_types = FALSE) |>
    filter(!is.na(tau_r1), !is.na(tau_var), !is.na(tau_skew))
  cat(sprintf("✓ %d pixels with valid tau_r1, tau_var, tau_skew\n", nrow(d)))

  # Drop pixels burned in 2000-2014: a fire crash isn't gradual stress
  n_before <- nrow(d)
  d <- d |>
    left_join(fire_flags, by = "pixel_id") |>
    filter(!burned_baseline_00_14) |>
    select(-burned_baseline_00_14)
  cat(sprintf("✓ Excluded %d pixels burned (MTBS moderate+) during 2000-2014 baseline (%d remain)\n",
              n_before - nrow(d), nrow(d)))

  # Remove the part of tau explained by baseline vegetation density
  base_level <- compute_base_level(cfg$mat_path, pixel_meta$pixel_id, "base_level")
  d <- d |> left_join(base_level, by = "pixel_id")
  stopifnot(!anyNA(d$base_level))

  d <- residualize_tau(d, "tau_r1", "tau_skew", "base_level",
                       "resid_tau_r1", "resid_tau_skew")

  # Skew is left unflipped: its direction near a transition isn't fixed
  feat <- d |> transmute(resid_tau_r1, resid_tau_skew) |> scale()

  cat("Running K sweep (K=2 to 8) for elbow and silhouette plots...\n")
  sweep_res <- run_kmeans_sweep(feat, K, SEED)
  sweep_df  <- sweep_res$sweep_df

  cat(sprintf("Running k-means (k=%d, seed=%d)...\n", K, SEED))
  km <- sweep_res$km_final
  d$cluster_raw <- km$cluster

  # Relabel clusters 1 (least stressed) .. K (most stressed)
  stress_rank <- d |>
    group_by(cluster_raw) |>
    summarise(composite = mean(resid_tau_r1) + mean(resid_tau_skew), .groups = "drop") |>
    arrange(composite) |>
    mutate(stress_rank = row_number())

  d <- d |>
    left_join(stress_rank |> select(cluster_raw, stress_rank),
              by = "cluster_raw") |>
    mutate(cluster = factor(stress_rank, levels = seq_len(K)))

  cat(sprintf("✓ Cluster sizes (stress rank 1=least → %d=most): %s\n",
              K, paste(table(d$cluster), collapse = " / ")))

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

# Joins one panel per index; returns a single panel unchanged
combine_panels <- function(panels, guides = "collect", combine_op = `+`) {
  if (length(panels) == 1) return(panels[[1]])
  Reduce(combine_op, panels) + plot_layout(guides = guides)
}

# Chart sizes scale with the number of indices
n_idx <- length(INDICES)

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

build_map_panel <- function(res) {
  build_cluster_map_panel(res$d, sn_boundary, cluster_cols, clust_labels, res$cfg$label)
}

map_panels <- setNames(lapply(INDICES, function(idx) build_map_panel(results[[idx]])), INDICES)

# Show the caption once under both panels instead of on each
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

  profile_path <- sprintf("pixel_approach/data/processed/cluster_profile_%s.csv", idx)
  write_csv(profile, profile_path)
  cat(sprintf("✓ Saved: %s\n", profile_path))
}
