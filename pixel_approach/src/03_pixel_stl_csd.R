library(tidyverse)
library(zoo)
library(future.apply)
library(progressr)

# progressr defaults to disabled when interactive() is FALSE (e.g. Rscript),
# so force it on to guarantee the progress bar shows regardless of how this
# script is run.
options(progressr.enable = TRUE)
handlers("txtprogressbar")

# --------------------------------------------------------------------------
# Parameters — set from console before source()ing to override defaults:
#   WIN           <- 54    # rolling window (months)
#   NA_THRESHOLD  <- 0.30  # skip pixel if fraction of NAs >= this
#   EDGE_NA_LIMIT <- 6     # skip pixel if leading/trailing NA run exceeds this
#                            (guards against na.approx(rule=2) flatlining a long
#                            edge gap — verified not an issue in the current
#                            MODIS NDVI data (worst case: 2-month leading /
#                            1-month trailing run), but a real hazard if this
#                            pipeline is reused with a sensor that has larger
#                            cloud-masking gaps, e.g. Landsat or Sentinel-2)
#   N_WORKERS     <- 5     # parallel worker processes (default: 5)
#   FUTURE_MAX_SIZE_MB <- 2048  # cap on globals shipped to each worker (MiB);
#                            must exceed ndvi_mat + vpd_mat + ppt_mat combined
#                            (~864 MiB currently) — future's own default cap
#                            is 500 MiB, meant to catch accidental large
#                            exports, not this intentional one
# --------------------------------------------------------------------------
if (!exists("WIN"))                WIN                <- 54
if (!exists("NA_THRESHOLD"))       NA_THRESHOLD       <- 0.30
if (!exists("EDGE_NA_LIMIT"))      EDGE_NA_LIMIT      <- 6
if (!exists("N_WORKERS"))          N_WORKERS          <- 5
if (!exists("FUTURE_MAX_SIZE_MB")) FUTURE_MAX_SIZE_MB <- 2048
options(future.globals.maxSize = FUTURE_MAX_SIZE_MB * 1024^2)

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
# roll_r1 is date-aware because Dec-Feb are already removed from the series
# (see gs_idx below), so Nov (year N) sits array-adjacent to Mar (year N+1) —
# 4 calendar months apart, not 1. Without this check that pair would be
# wrongly scored as a normal 1-month lag.
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

# --- Data guardrail: contiguous edge-NA run length -------------------------
# na.approx(rule=2) fills leading/trailing NAs by copying the nearest real
# value, which flatlines a long edge run. Checked empirically against the
# current MODIS NDVI data: worst case across all 251,971 pixels is a 2-month
# leading / 1-month trailing run — not a live problem here. Kept as a guard
# in case this pipeline is reused with a sensor that has larger cloud-masking
# gaps (e.g. Landsat, Sentinel-2).
edge_na_run <- function(x) {
  na_idx <- is.na(x)
  if (!na_idx[1]) return(0L)
  rle(na_idx)$lengths[1]
}

# Growing-season index — filter to Mar–Nov after STL (same as cluster-level
# CSD in src/11_1_csd_manual.R). STL still runs on all 12 months.
col_dates  <- as.Date(colnames(ndvi_mat))
col_months <- as.integer(format(col_dates, "%m"))
gs_idx     <- which(col_months %in% 3:11)
gs_dates   <- col_dates[gs_idx]

# vpd_mat/ppt_mat are only ever read via [i, gs_idx] below — trim to growing
# season columns once here so each parallel worker gets ~25% less data.
vpd_mat <- vpd_mat[, gs_idx]
ppt_mat <- ppt_mat[, gs_idx]

# Baseline mask for climate regression: 2000–2011 pre-drought growing season
baseline_gs <- gs_dates <= as.Date("2011-12-31")

# --- Per-pixel CSD --------------------------------------------------------
na_limit <- NA_THRESHOLD * n_months
n_gs     <- length(gs_idx)   # fixed for every pixel (135 months)

# Rolling-window row indices and their dates are identical for every pixel
# (WIN and n_gs are both fixed), so compute them once instead of rebuilding
# the same (j-WIN+1):j slice and gs_dates lookup 251,971 x 82 times.
window_idx   <- lapply(WIN:n_gs, function(j) (j - WIN + 1):j)
window_dates <- lapply(window_idx, function(idx) gs_dates[idx])

# seas_01..seas_12 (periodic monthly seasonal offsets) and b0/b1/b2 (fitted
# climate-model coefficients) are exported for downstream reuse without
# re-running STL/lm from scratch.
na_row <- c(tau_r1 = NA_real_, p_r1 = NA_real_, tau_var = NA_real_,
            p_var = NA_real_, tau_skew = NA_real_, p_skew = NA_real_,
            setNames(rep(NA_real_, 12), sprintf("seas_%02d", 1:12)),
            b0 = NA_real_, b1 = NA_real_, b2 = NA_real_)

# Every pixel's 7-step process (interpolate -> STL -> deseasonalize -> filter
# to growing season -> climate correction -> rolling CSD stats -> Kendall's
# tau) is fully independent of every other pixel, so this is embarrassingly
# parallel. Each call returns a status string plus its result row (NA row if
# dropped); both are tabulated together after future_lapply() finishes.
process_pixel <- function(i) {
  p()
  y <- ndvi_mat[i, ]

  # skip pixels with too many NAs — interpolation would be mostly fabricated
  if (sum(is.na(y)) >= na_limit)
    return(list(status = "dropped_high_na", row = na_row))

  # skip pixels with a long leading/trailing NA run — na.approx(rule=2) would
  # flatline it, distorting rolling var/r1 at the flat-to-real transition
  if (edge_na_run(y) > EDGE_NA_LIMIT || edge_na_run(rev(y)) > EDGE_NA_LIMIT)
    return(list(status = "dropped_edge_na", row = na_row))

  # linear interpolation for remaining NAs; rule=2 fills edge NAs by extension
  y <- as.numeric(zoo::na.approx(y, na.rm = FALSE, rule = 2))

  # STL — monthly series starting Feb 2000, frequency 12
  ts_y   <- ts(y, start = c(2000, 2), frequency = 12)
  decomp <- tryCatch(
    stl(ts_y, s.window = "periodic", robust = TRUE),
    error = function(e) NULL
  )
  if (is.null(decomp))
    return(list(status = "dropped_stl_fail", row = na_row))

  # ndvi_deseasoned = ndvi - seasonal component (keeps trend + remainder).
  # Same approach as cluster-level src/11_1_csd_manual.R.
  seasonal    <- as.numeric(decomp$time.series[, "seasonal"])
  ndvi_deseas <- y - seasonal

  # 12 monthly seasonal offsets (s.window="periodic" -> exactly one value per
  # calendar month, verified empirically). col_months is the global calendar-
  # month vector for every column of ndvi_mat (same for every pixel).
  seas_by_month <- seasonal[match(1:12, col_months)]
  names(seas_by_month) <- sprintf("seas_%02d", 1:12)

  # Filter to growing season (vpd_mat/ppt_mat pre-trimmed to gs_idx columns)
  ndvi_gs <- ndvi_deseas[gs_idx]
  vpd_gs  <- vpd_mat[i, ]
  ppt_gs  <- ppt_mat[i, ]

  # Climate correction: fit lm on 2000–2011 baseline, apply to full series
  lm_data    <- data.frame(y = ndvi_gs, vpd = vpd_gs, ppt = ppt_gs)
  clim_fit   <- tryCatch(
    lm(y ~ vpd + ppt, data = lm_data[baseline_gs, ]),
    error = function(e) NULL
  )
  if (is.null(clim_fit))
    return(list(status = "dropped_clim_fail", row = na_row))
  clim_resid <- ndvi_gs - predict(clim_fit, newdata = lm_data)

  n <- length(clim_resid)

  # rolling stats over WIN-month sliding window (window_idx/window_dates
  # precomputed once outside this function — same for every pixel)
  r1_vec   <- rep(NA_real_, n)
  var_vec  <- rep(NA_real_, n)
  skew_vec <- rep(NA_real_, n)

  for (k in seq_along(window_idx)) {
    j           <- WIN + k - 1
    w           <- clim_resid[window_idx[[k]]]
    r1_vec[j]   <- roll_r1(w, window_dates[[k]])
    var_vec[j]  <- roll_var(w)
    skew_vec[j] <- roll_skew(w)
  }

  valid <- which(!is.na(r1_vec))
  if (length(valid) < 10)
    return(list(status = "dropped_insufficient_valid", row = na_row))

  idx    <- seq_along(valid)
  t_r1   <- cor.test(idx, r1_vec[valid],   method = "kendall")
  t_var  <- cor.test(idx, var_vec[valid],  method = "kendall")
  t_skew <- cor.test(idx, skew_vec[valid], method = "kendall")

  # unname() strips cor.test()$estimate's own "tau" name — otherwise c(tau_r1 = ...)
  # concatenates it into "tau_r1.tau" instead of the plain column name "tau_r1"
  # (same reason coef(clim_fit)[...] below is also unname()'d)
  coefs <- coef(clim_fit)
  list(
    status = "ok",
    row = c(
      tau_r1   = unname(round(t_r1$estimate,   4)), p_r1   = round(t_r1$p.value,   4),
      tau_var  = unname(round(t_var$estimate,  4)), p_var  = round(t_var$p.value,  4),
      tau_skew = unname(round(t_skew$estimate, 4)), p_skew = round(t_skew$p.value, 4),
      round(seas_by_month, 4),
      b0 = unname(round(coefs["(Intercept)"], 4)),
      b1 = unname(round(coefs["vpd"], 4)),
      b2 = unname(round(coefs["ppt"], 4))
    )
  )
}

# R here is linked against Apple's Accelerate/vecLib BLAS, which multi-threads
# lm()'s linear algebra by default. Left uncapped, each of the N_WORKERS
# processes could spawn its own BLAS threads on top of the per-pixel
# parallelism already provided by future_lapply(), oversubscribing the
# machine's cores instead of speeding things up. Force each worker's BLAS
# to single-threaded; the parallelism we want already comes from splitting
# pixels across workers. Set before plan() so workers inherit it at spawn.
Sys.setenv(VECLIB_MAXIMUM_THREADS = "1")

future::plan(future::multisession, workers = N_WORKERS)

cat(sprintf("\nRunning STL + climate correction + CSD for %d pixels across %d workers...\n",
            n_pixels, N_WORKERS))
t_start <- proc.time()["elapsed"]

with_progress({
  p <- progressor(along = seq_len(n_pixels))
  # future.seed = NULL disables future's RNG-misuse check. process_pixel() itself
  # never touches R's RNG (verified empirically: ran all steps on 20,000 real
  # pixels sequentially, checked .Random.seed before/after each — zero changes),
  # so the "UNRELIABLE VALUE" warning some runs produce is a false positive from
  # the multisession/progressr worker machinery itself, not from this function.
  pixel_results <- future_lapply(seq_len(n_pixels), process_pixel,
                                  future.seed = NULL,
                                  future.packages = "zoo")
})

future::plan(future::sequential)  # shut down worker processes

elapsed_total <- (proc.time()["elapsed"] - t_start) / 60

pixel_status <- vapply(pixel_results, `[[`, character(1), "status")
result_rows  <- do.call(rbind, lapply(pixel_results, `[[`, "row"))

results <- data.frame(
  pixel_id = pixel_meta$pixel_id,
  lon      = pixel_meta$lon,
  lat      = pixel_meta$lat,
  result_rows
)

cat(sprintf("\n✓ Completed in %.1f minutes\n", elapsed_total))
cat("  --- Summary ---\n")
cat(sprintf("  Total pixels evaluated:               %d\n", n_pixels))
cat(sprintf("  Dropped (missing > %.0f%%):             %d\n",
            NA_THRESHOLD * 100, sum(pixel_status == "dropped_high_na")))
cat(sprintf("  Dropped (edge NA run > %d months):     %d\n",
            EDGE_NA_LIMIT, sum(pixel_status == "dropped_edge_na")))
cat(sprintf("  Dropped (STL/climate fit failed):      %d\n",
            sum(pixel_status %in% c("dropped_stl_fail", "dropped_clim_fail",
                                     "dropped_insufficient_valid"))))
cat(sprintf("  Final pixels processed:                %d\n",
            sum(pixel_status == "ok")))

out_file <- "pixel_approach/data/processed/pixel_tau.csv"
write_csv(results, out_file)
cat(sprintf("✓ Saved: %s\n", out_file))
