library(tidyverse)

# Adds an "actual drop" ground-truth column to each index's pixel_tau file:
#   Actual_*_Drop_15_16 = mean(Mar-Nov, 2000-2011) - min(Mar-Nov, 2015-2016)
# Growing season only (Mar-Nov) to match the 9-month CSD computation and
# avoid winter/snow contamination. Positive = collapsed more.
#
# NDWI gets its own Actual_NDWI_Drop_15_16 (computed identically, from its
# own matrix) rather than reusing NDVI's — the two clustering pipelines are
# independent, so each index's tau is validated against its own future
# collapse, matching the standard CSD-literature approach (does a variable's
# own destabilization signal forecast its own future collapse).
#
# DEPRIORITIZED: real-world validation (MTBS fire, LCMAP mortality, USFS
# ADS beetle -- 09/10/08) is now the trusted ground truth; Actual_Drop
# doesn't agree with it on cluster severity ordering and is only
# growing-season, unlike the other three. 06_clustering.R's ANOVA/Tukey/
# ROC block that validated against this column is commented out. This
# script is not auto-invoked by anything (no master pipeline script calls
# it -- only pixel_approach/tests/test-04-actual-drop.R does), so nothing
# breaks by leaving it alone. It's still needed, though: 05_heatmaps.R and
# 07_regression.R both still read Actual_Drop_15_16 from the tau CSVs this
# script produces -- if 03_pixel_stl_csd.R is ever re-run (which rebuilds
# those CSVs from scratch without this column), re-run this script
# afterward or those two will break.

# --------------------------------------------------------------------------
# Parameters — set from console before source()ing to override defaults:
#   INDICES <- c("ndvi", "ndwi")  # set to c("ndvi") or c("ndwi") to run only one
# --------------------------------------------------------------------------
if (!exists("INDICES")) INDICES <- c("ndvi", "ndwi")

index_config <- list(
  ndvi = list(
    tau_path = "pixel_approach/data/processed/ndvi_pixel_tau.csv",
    mat_path = "pixel_approach/data/processed/ndvi_matrix.rds",
    drop_col = "Actual_Drop_15_16",
    label    = "NDVI"
  ),
  ndwi = list(
    tau_path = "pixel_approach/data/processed/ndwi_pixel_tau.csv",
    mat_path = "pixel_approach/data/processed/ndwi_matrix.rds",
    drop_col = "Actual_NDWI_Drop_15_16",
    label    = "NDWI"
  )
)

# min(x, na.rm=TRUE) on an all-NA row returns Inf with a warning instead of NA.
# safe_min() avoids that — needed since some pixels have zero valid months in
# the 2015-2016 window despite having enough 2000-2014 data to pass 03's checks.
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)

add_actual_drop <- function(tau_path, mat_path, drop_col, var_label) {
  cat(sprintf("\n=== %s ===\n", var_label))
  cat("Loading inputs...\n")
  pixel_tau <- read_csv(tau_path, show_col_types = FALSE)
  mat       <- readRDS(mat_path)

  col_dates  <- colnames(mat)
  col_months <- as.integer(format(as.Date(col_dates), "%m"))
  gs_months  <- col_months %in% 3:11   # growing season: Mar–Nov

  # Baseline: 2000-2011, growing season only
  baseline_cols <- col_dates >= "2000-02-01" & col_dates <= "2011-12-01" & gs_months
  cat(sprintf("✓ Baseline columns (2000-2011, Mar-Nov): %d months\n", sum(baseline_cols)))

  # Collapse window: 2015-2016, growing season only
  drop_cols <- col_dates >= "2015-01-01" & col_dates <= "2016-12-01" & gs_months
  cat(sprintf("✓ Collapse columns (2015-2016, Mar-Nov): %d months\n", sum(drop_cols)))

  # Compute per pixel
  cat("Computing baseline mean and 2015-2016 min...\n")
  baseline_mean <- rowMeans(mat[, baseline_cols, drop = FALSE], na.rm = TRUE)
  min_2015_2016 <- apply(mat[, drop_cols, drop = FALSE], 1, safe_min)

  # pixel_tau[[drop_col]] below is a positional assignment (no join) —
  # relies entirely on pixel_tau and mat sharing the same pixel_id row
  # order. That's true today (both come from the same 1-indexed pipeline),
  # but nothing else here checks it; a silent misalignment here would
  # assign every drop value to the wrong pixel with no error at all.
  # (Uses == rather than identical(): pixel_id round-trips through CSV as
  # double, not integer, so identical() against seq_len()'s integer output
  # would false-positive on a type mismatch even when values match.)
  stopifnot(all(pixel_tau$pixel_id == seq_len(nrow(mat))))

  pixel_tau[[drop_col]] <- round(baseline_mean - min_2015_2016, 4)

  cat(sprintf("✓ %s range: %.4f to %.4f (mean: %.4f)\n", drop_col,
              min(pixel_tau[[drop_col]], na.rm = TRUE),
              max(pixel_tau[[drop_col]], na.rm = TRUE),
              mean(pixel_tau[[drop_col]], na.rm = TRUE)))

  write_csv(pixel_tau, tau_path)
  cat(sprintf("✓ Saved: %s\n", tau_path))
}

for (idx in INDICES) {
  cfg <- index_config[[idx]]
  add_actual_drop(cfg$tau_path, cfg$mat_path, cfg$drop_col, cfg$label)
}
