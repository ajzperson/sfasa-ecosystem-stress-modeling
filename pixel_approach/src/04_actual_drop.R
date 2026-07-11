library(tidyverse)

# Adds Actual_Drop_15_16 to pixel_tau.csv:
#   Actual_Drop_15_16 = mean_NDVI(Mar-Nov, 2000-2011) - min_NDVI(Mar-Nov, 2015-2016)
# Growing season only (Mar-Nov) to match the 9-month CSD computation and
# avoid winter/snow contamination. Positive = collapsed more.

# min(x, na.rm=TRUE) on an all-NA row returns Inf with a warning instead of NA.
# safe_min() avoids that — needed since some pixels have zero valid NDVI months
# in the 2015-2016 window despite having enough 2000-2014 data to pass 03's checks.
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)

cat("Loading inputs...\n")
pixel_tau <- read_csv("pixel_approach/data/processed/pixel_tau.csv",
                      show_col_types = FALSE) |>
  select(-any_of("Corrected_Drop_15_16"))  # drop stale column from older runs
ndvi_mat  <- readRDS("pixel_approach/data/processed/ndvi_matrix.rds")

col_dates  <- colnames(ndvi_mat)
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
baseline_mean  <- rowMeans(ndvi_mat[, baseline_cols], na.rm = TRUE)
min_2015_2016  <- apply(ndvi_mat[, drop_cols], 1, safe_min)

pixel_tau$Actual_Drop_15_16 <- round(baseline_mean - min_2015_2016, 4)

cat(sprintf("✓ Actual_Drop_15_16 range: %.4f to %.4f (mean: %.4f)\n",
            min(pixel_tau$Actual_Drop_15_16, na.rm = TRUE),
            max(pixel_tau$Actual_Drop_15_16, na.rm = TRUE),
            mean(pixel_tau$Actual_Drop_15_16, na.rm = TRUE)))

write_csv(pixel_tau, "pixel_approach/data/processed/pixel_tau.csv")
cat("✓ Saved: pixel_approach/data/processed/pixel_tau.csv\n")
