library(tidyverse)
library(sf)
library(patchwork)
library(cluster)

source("pixel_approach/src/clustering_helpers.R")

# Exploratory alternative to 06_clustering.R: instead of clustering NDVI and
# NDWI independently (2 features each: density-residualized tau_r1 +
# tau_skew), combines both indices' features into one 4-feature K-means
# clustering. Does not modify 06_clustering.R or its outputs -- this is a
# separate approach, run and compared alongside it.
#
# Output:
#   pixel_approach/data/processed/combined_pixel_clusters.csv
#   pixel_approach/data/processed/combined_cluster_profile.csv
#   pixel_approach/charts/combined_cluster_k_sweep.pdf
#   pixel_approach/charts/combined_cluster_map.pdf

if (!exists("K"))              K              <- 4
if (!exists("SEED"))           SEED           <- 42
if (!exists("NDVI_TAU_FILE"))  NDVI_TAU_FILE  <- "pixel_approach/data/processed/ndvi_pixel_tau.csv"
if (!exists("NDWI_TAU_FILE"))  NDWI_TAU_FILE  <- "pixel_approach/data/processed/ndwi_pixel_tau.csv"

dir.create("pixel_approach/charts", showWarnings = FALSE, recursive = TRUE)

sn_boundary <- st_read("docs_common/sierra_nevada_boundary.geojson", quiet = TRUE)

fire_flags <- read_csv("pixel_approach/data/processed/pixel_fire.csv", show_col_types = FALSE) |>
  select(pixel_id, burned_baseline_00_14)

pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv", show_col_types = FALSE) |>
  select(pixel_id)

style        <- build_cluster_style(K)
cluster_cols <- style$cols
clust_labels <- style$labels

cat(sprintf("Loading NDVI and NDWI tau outputs (%s, %s)...\n", NDVI_TAU_FILE, NDWI_TAU_FILE))
ndvi_tau <- read_csv(NDVI_TAU_FILE, show_col_types = FALSE) |>
  filter(!is.na(tau_r1), !is.na(tau_skew)) |>
  select(pixel_id, lon, lat, tau_r1, tau_skew)

ndwi_tau <- read_csv(NDWI_TAU_FILE, show_col_types = FALSE) |>
  filter(!is.na(tau_r1), !is.na(tau_skew)) |>
  select(pixel_id, tau_r1, tau_skew)

# Inner join -- NDVI and NDWI don't drop exactly the same pixels in 03 (each
# has its own missing-data pattern), so only pixels with a valid tau for
# BOTH indices can go into a combined clustering.
d <- ndvi_tau |>
  inner_join(ndwi_tau, by = "pixel_id", suffix = c("_ndvi", "_ndwi"))
cat(sprintf("✓ %d pixels with valid tau_r1/tau_skew for both NDVI and NDWI\n", nrow(d)))

n_before <- nrow(d)
d <- d |>
  left_join(fire_flags, by = "pixel_id") |>
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

# Each index's tau_r1/tau_skew residualized against its OWN baseline density
# (not the other index's) -- same reasoning as 06_clustering.R, applied twice.
d <- residualize_tau(d, "tau_r1_ndvi", "tau_skew_ndvi", "base_level_ndvi",
                     "resid_tau_r1_ndvi", "resid_tau_skew_ndvi")
d <- residualize_tau(d, "tau_r1_ndwi", "tau_skew_ndwi", "base_level_ndwi",
                     "resid_tau_r1_ndwi", "resid_tau_skew_ndwi")

# Each index's skew residual is used unflipped, same reasoning as
# 06_clustering.R: skewness's direction near a transition isn't forced to
# match resid_tau_r1's "higher = more stressed" convention, since it depends
# on which state the system is heading toward (not universal, unlike r1).
feat <- d |>
  transmute(resid_tau_r1_ndvi,
            resid_tau_skew_ndvi,
            resid_tau_r1_ndwi,
            resid_tau_skew_ndwi) |>
  scale()

cat("Running K sweep (K=2 to 8) for elbow and silhouette plots...\n")
sweep_res <- run_kmeans_sweep(feat, K, SEED)
sweep_df  <- sweep_res$sweep_df

cat(sprintf("Running k-means (k=%d, seed=%d) on combined 4-feature space...\n", K, SEED))
km <- sweep_res$km_final
d$cluster_raw <- km$cluster

# Composite/stress-rank score: unweighted sum of the same 4 z-scored
# features used for clustering. An earlier version used PC1 (first
# principal component) instead, reasoned that NDVI/NDWI's weak cross-index
# correlation (r1: r=0.20; skew: r=-0.15) meant equal weights weren't
# supported by the data. Tested directly: PC1 only explains 30.8% of
# variance (barely above the 25% no-structure baseline for 4 features) and
# assigns tau_skew_ndvi a counterintuitive *negative* loading (-0.293) --
# despite that, re-ranking the same raw K-means clusters with PC1 vs. a
# plain unweighted sum produces 100% identical cluster assignments (verified
# directly: both orderings are cluster 1 < 4 < 3 < 2 by mean score). Since
# the two methods are provably equivalent for this data, the simpler,
# easier-to-explain unweighted sum is used -- no PCA-specific caveat needed.
d$stress_score <- rowSums(feat)

stress_rank <- d |>
  group_by(cluster_raw) |>
  summarise(composite = mean(stress_score), .groups = "drop") |>
  arrange(composite) |>
  mutate(stress_rank = row_number())

d <- d |>
  left_join(stress_rank |> select(cluster_raw, stress_rank), by = "cluster_raw") |>
  mutate(cluster = factor(stress_rank, levels = seq_len(K)))

cat(sprintf("✓ Cluster sizes (stress rank 1=least → %d=most): %s\n",
            K, paste(table(d$cluster), collapse = " / ")))

write_csv(d, "pixel_approach/data/processed/combined_pixel_clusters.csv")
cat("✓ Saved: pixel_approach/data/processed/combined_pixel_clusters.csv\n")

# --- K-sweep chart ---
p_sweep <- build_k_sweep_panel(sweep_df, K, "NDVI+NDWI-combined") +
  plot_annotation(
    title = sprintf("K-means Hyperparameter Sweep (K = 2-8), Combined NDVI+NDWI Features; selected K = %d shown in red", K),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "#1B2A1F"))
  )

ggsave("pixel_approach/charts/combined_cluster_k_sweep.pdf", plot = p_sweep, width = 10, height = 4, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_cluster_k_sweep.pdf\n")
# --- Cluster map ---
p_map <- build_cluster_map_panel(
  d, sn_boundary, cluster_cols, clust_labels,
  sprintf("Combined NDVI+NDWI Clustering (k=%d) on 4 density-residualized features", K)
)

ggsave("pixel_approach/charts/combined_cluster_map.pdf", plot = p_map, width = 7, height = 7, dpi = 150)
cat("✓ Saved: pixel_approach/charts/combined_cluster_map.pdf\n")

# --- Cluster profile ---
profile <- d |>
  group_by(cluster) |>
  summarise(
    n                        = n(),
    tau_r1_ndvi_mean         = round(mean(tau_r1_ndvi), 3),
    tau_skew_ndvi_mean       = round(mean(tau_skew_ndvi), 3),
    tau_r1_ndwi_mean         = round(mean(tau_r1_ndwi), 3),
    tau_skew_ndwi_mean       = round(mean(tau_skew_ndwi), 3),
    resid_tau_r1_ndvi_mean   = round(mean(resid_tau_r1_ndvi), 3),
    resid_tau_skew_ndvi_mean = round(mean(resid_tau_skew_ndvi), 3),
    resid_tau_r1_ndwi_mean   = round(mean(resid_tau_r1_ndwi), 3),
    resid_tau_skew_ndwi_mean = round(mean(resid_tau_skew_ndwi), 3),
    .groups = "drop"
  ) |>
  arrange(as.integer(cluster))

cat("\n=== Combined cluster profile (1=least stressed, K=most stressed) ===\n")
print(profile)

write_csv(profile, "pixel_approach/data/processed/combined_cluster_profile.csv")
cat("✓ Saved: pixel_approach/data/processed/combined_cluster_profile.csv\n")