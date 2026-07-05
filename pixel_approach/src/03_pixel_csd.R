library(tidyverse)
library(zoo)

# --------------------------------------------------------------------------
# Parameters — set from console before source()ing to override defaults:
#   WIN          <- 54    # rolling window (months)
#   NA_THRESHOLD <- 0.30  # skip pixel if fraction of NAs >= this
# --------------------------------------------------------------------------
if (!exists("WIN"))          WIN          <- 54
if (!exists("NA_THRESHOLD")) NA_THRESHOLD <- 0.30

# --- Load inputs ----------------------------------------------------------
cat("Loading inputs...\n")
pixel_meta <- read_csv("pixel_approach/data/processed/pixel_metadata.csv",
                       show_col_types = FALSE)
ndvi_mat   <- readRDS("pixel_approach/data/processed/ndvi_matrix.rds")
vpd_mat    <- readRDS("pixel_approach/data/processed/vpd_matrix.rds")
ppt_mat    <- readRDS("pixel_approach/data/processed/ppt_matrix.rds")

# Trim ndvi to 2000-2014; align vpd/ppt to same columns by name
# (vpd/ppt matrices are already 2000-2014, so keep_cols would be wrong length)
keep_cols <- colnames(ndvi_mat) <= "2014-12-01"
ndvi_mat  <- ndvi_mat[, keep_cols]
vpd_mat   <- vpd_mat[, colnames(ndvi_mat)]
ppt_mat   <- ppt_mat[, colnames(ndvi_mat)]
n_months  <- ncol(ndvi_mat)   # 179 (Feb 2000 – Dec 2014)
n_pixels  <- nrow(ndvi_mat)   # 251,971

cat(sprintf("✓ Analysis window: %s to %s (%d months)\n",
            colnames(ndvi_mat)[1], colnames(ndvi_mat)[n_months], n_months))
cat(sprintf("✓ Pixels: %d | WIN: %d months | NA threshold: %.0f%%\n",
            n_pixels, WIN, NA_THRESHOLD * 100))

# --- Rolling stat functions -----------------------------------------------
roll_r1 <- function(x, d) {
  xbar <- mean(x)
  midx <- as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m"))
  keep <- diff(midx) == 1L   # only genuine 1-month-apart pairs (excludes Nov→Mar gap)
  num  <- sum((head(x, -1)[keep] - xbar) * (tail(x, -1)[keep] - xbar)) / length(x)
  den  <- mean((x - xbar)^2)
  num / den
}

roll_var <- function(x) var(x)

roll_skew <- function(x) {
  mu  <- mean(x)
  num <- mean((x - mu)^3)
  den <- (mean((x - mu)^2))^(3/2)
  -(num / den)
}

# Growing-season index — filter to Mar–Nov after STL (same as cluster-level
# CSD in src/11_1_csd_manual.R). STL still runs on all 12 months.
col_dates  <- as.Date(colnames(ndvi_mat))
col_months <- as.integer(format(col_dates, "%m"))
gs_idx     <- which(col_months %in% 3:11)
gs_dates   <- col_dates[gs_idx]

# Baseline mask for climate regression: 2000–2011 pre-drought growing season
baseline_gs <- gs_dates <= as.Date("2011-12-31")

# --- Per-pixel CSD --------------------------------------------------------
na_limit <- NA_THRESHOLD * n_months

results <- data.frame(
  pixel_id = pixel_meta$pixel_id,
  lon      = pixel_meta$lon,
  lat      = pixel_meta$lat,
  tau_r1   = NA_real_,  p_r1   = NA_real_,
  tau_var  = NA_real_,  p_var  = NA_real_,
  tau_skew = NA_real_,  p_skew = NA_real_
)

n_dropped <- sum(apply(ndvi_mat, 1, function(y) sum(is.na(y))) >= na_limit)
cat(sprintf("  Pixels dropped (>= %.0f%% NA): %d / %d\n",
            NA_THRESHOLD * 100, n_dropped, n_pixels))

cat(sprintf("\nRunning STL + climate correction + CSD for %d pixels...\n", n_pixels))
t_start <- proc.time()["elapsed"]

for (i in seq_len(n_pixels)) {

  if (i %% 10000 == 0) {
    elapsed <- proc.time()["elapsed"] - t_start
    cat(sprintf("  [%d / %d] %.0f%% — %.1f min elapsed\n",
                i, n_pixels, 100 * i / n_pixels, elapsed / 60))
  }

  y <- ndvi_mat[i, ]

  # skip pixels with too many NAs — interpolation would be mostly fabricated
  if (sum(is.na(y)) >= na_limit) next

  # linear interpolation for remaining NAs; rule=2 fills edge NAs by extension
  y <- as.numeric(na.approx(y, na.rm = FALSE, rule = 2))

  # STL — monthly series starting Feb 2000, frequency 12
  ts_y   <- ts(y, start = c(2000, 2), frequency = 12)
  decomp <- tryCatch(
    stl(ts_y, s.window = "periodic", robust = TRUE),
    error = function(e) NULL
  )
  if (is.null(decomp)) next

  # ndvi_deseasoned = ndvi - seasonal component (keeps trend + remainder).
  # Same approach as cluster-level src/11_1_csd_manual.R.
  seasonal    <- as.numeric(decomp$time.series[, "seasonal"])
  ndvi_deseas <- y - seasonal

  # Filter to growing season
  ndvi_gs <- ndvi_deseas[gs_idx]
  vpd_gs  <- vpd_mat[i, gs_idx]
  ppt_gs  <- ppt_mat[i, gs_idx]

  # Climate correction: fit lm on 2000–2011 baseline, apply to full series
  lm_data    <- data.frame(y = ndvi_gs, vpd = vpd_gs, ppt = ppt_gs)
  clim_fit   <- tryCatch(
    lm(y ~ vpd + ppt, data = lm_data[baseline_gs, ]),
    error = function(e) NULL
  )
  if (is.null(clim_fit)) next
  clim_resid <- ndvi_gs - predict(clim_fit, newdata = lm_data)

  n <- length(clim_resid)

  # rolling stats over WIN-month sliding window
  r1_vec   <- rep(NA_real_, n)
  var_vec  <- rep(NA_real_, n)
  skew_vec <- rep(NA_real_, n)

  for (j in WIN:n) {
    idx_w        <- (j - WIN + 1):j
    w            <- clim_resid[idx_w]
    r1_vec[j]   <- roll_r1(w, gs_dates[idx_w])
    var_vec[j]  <- roll_var(w)
    skew_vec[j] <- roll_skew(w)
  }

  valid <- which(!is.na(r1_vec))
  if (length(valid) < 10) next

  idx    <- seq_along(valid)
  t_r1   <- cor.test(idx, r1_vec[valid],   method = "kendall")
  t_var  <- cor.test(idx, var_vec[valid],  method = "kendall")
  t_skew <- cor.test(idx, skew_vec[valid], method = "kendall")

  results$tau_r1[i]   <- round(t_r1$estimate,   4)
  results$p_r1[i]     <- round(t_r1$p.value,    4)
  results$tau_var[i]  <- round(t_var$estimate,  4)
  results$p_var[i]    <- round(t_var$p.value,   4)
  results$tau_skew[i] <- round(t_skew$estimate, 4)
  results$p_skew[i]   <- round(t_skew$p.value,  4)
}

elapsed_total <- (proc.time()["elapsed"] - t_start) / 60
cat(sprintf("\n✓ Completed in %.1f minutes\n", elapsed_total))
cat(sprintf("  Pixels with valid tau: %d / %d\n",
            sum(!is.na(results$tau_r1)), n_pixels))

out_file <- "pixel_approach/data/processed/pixel_tau.csv"
write_csv(results, out_file)
cat(sprintf("✓ Saved: %s\n", out_file))
