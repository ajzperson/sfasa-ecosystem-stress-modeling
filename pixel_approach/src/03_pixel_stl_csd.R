library(tidyverse)
library(zoo)
library(future.apply)
library(progressr)

# Computes CSD (critical slowing down) trend statistics per pixel: STL-
# deseasonalize -> climate-correct (VPD+PPT, fit on the 2000-2011 baseline)
# -> rolling WIN-month r1/variance/skewness over the growing season (Mar-
# Nov) -> Kendall's tau trend in each. Runs once per index in INDICES
# (NDVI and/or NDWI), parallelized across pixels with future.apply.
#
# Output (per index):
#   pixel_approach/data/processed/<index>_pixel_tau.csv
#     (pixel_id, lon, lat, tau_r1, p_r1, tau_var, p_var, tau_skew, p_skew,
#      seas_01..seas_12, b0, b1, b2)
#   pixel_approach/data/processed/<index>_rolling_r1.rds
#   pixel_approach/data/processed/<index>_rolling_var.rds
#   pixel_approach/data/processed/<index>_rolling_skew.rds

# progressr defaults to disabled when interactive() is FALSE (e.g. Rscript),
# so force it on to guarantee the progress bar shows regardless of how this
# script is run.
options(progressr.enable = TRUE)
handlers("txtprogressbar")

# --------------------------------------------------------------------------
# Parameters — set from console before source()ing to override defaults:
#   INDICES       <- c("ndvi", "ndwi")  # which index/indices to process;
#                            set to c("ndvi") or c("ndwi") to run only one
#   WIN           <- 54    # rolling window (months)
#   NA_THRESHOLD  <- 0.30  # skip pixel if fraction of NAs >= this
#   EDGE_NA_LIMIT <- 6     # skip pixel if leading/trailing NA run exceeds this
#                            (guards against na.approx(rule=2) flatlining a
#                            long edge gap -- see edge_na_run() below)
#   N_WORKERS     <- 5     # parallel worker processes (default: 5)
#   FUTURE_MAX_SIZE_MB <- 2048  # cap on globals shipped to each worker (MiB);
#                            must exceed mat + vpd_mat + ppt_mat combined
#                            (~864 MiB currently) — future's own default cap
#                            is 500 MiB, meant to catch accidental large
#                            exports, not this intentional one
# --------------------------------------------------------------------------
if (!exists("INDICES"))            INDICES            <- c("ndvi", "ndwi")
if (!exists("WIN"))                WIN                <- 54
if (!exists("NA_THRESHOLD"))       NA_THRESHOLD       <- 0.30
if (!exists("EDGE_NA_LIMIT"))      EDGE_NA_LIMIT      <- 6
if (!exists("N_WORKERS"))          N_WORKERS          <- 5
if (!exists("FUTURE_MAX_SIZE_MB")) FUTURE_MAX_SIZE_MB <- 2048
options(future.globals.maxSize = FUTURE_MAX_SIZE_MB * 1024^2)

index_config <- list(
  ndvi = list(
    matrix_path      = "pixel_approach/data/processed/ndvi_matrix.rds",
    out_file         = "pixel_approach/data/processed/ndvi_pixel_tau.csv",
    rolling_r1_file  = "pixel_approach/data/processed/ndvi_rolling_r1.rds",
    rolling_var_file = "pixel_approach/data/processed/ndvi_rolling_var.rds",
    rolling_skew_file = "pixel_approach/data/processed/ndvi_rolling_skew.rds",
    label       = "NDVI"
  ),
  ndwi = list(
    matrix_path      = "pixel_approach/data/processed/ndwi_matrix.rds",
    out_file         = "pixel_approach/data/processed/ndwi_pixel_tau.csv",
    rolling_r1_file  = "pixel_approach/data/processed/ndwi_rolling_r1.rds",
    rolling_var_file = "pixel_approach/data/processed/ndwi_rolling_var.rds",
    rolling_skew_file = "pixel_approach/data/processed/ndwi_rolling_skew.rds",
    label       = "NDWI"
  )
)

# --- Load shared inputs (same pixel grid/climate data for every index) ----
cat("Loading pixel metadata and climate inputs...\n")
pixel_meta   <- read_csv("pixel_approach/data/processed/pixel_metadata.csv",
                         show_col_types = FALSE)
vpd_mat_full <- readRDS("pixel_approach/data/processed/vpd_matrix.rds")
ppt_mat_full <- readRDS("pixel_approach/data/processed/ppt_matrix.rds")

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
  num / den
}

# --- Data guardrail: contiguous edge-NA run length -------------------------
# na.approx(rule=2) fills leading/trailing NAs by copying the nearest real
# value, which flatlines a long edge run. Not a live problem in the current
# MODIS data (worst case: 2-month leading / 1-month trailing run), but kept
# as a guard for a sensor with larger cloud-masking gaps (e.g. Landsat,
# Sentinel-2).
edge_na_run <- function(x) {
  na_idx <- is.na(x)
  if (!na_idx[1]) return(0L)
  rle(na_idx)$lengths[1]
}

# Runs the full STL + climate-correction + rolling-CSD pipeline for one
# vegetation index. Generic over `mat` — nothing below assumes NDVI
# specifically, so the same code produces pixel_tau.csv for NDVI and
# ndwi_pixel_tau.csv for NDWI.
run_pixel_csd <- function(mat, var_label, out_file,
                           rolling_r1_file, rolling_var_file, rolling_skew_file,
                           vpd_mat_full, ppt_mat_full, pixel_meta,
                           WIN, NA_THRESHOLD, EDGE_NA_LIMIT, N_WORKERS) {

  cat(sprintf("\n=== %s ===\n", var_label))

  # Trim to 2000-2014; align vpd/ppt to the same columns by name (vpd/ppt
  # cover Feb 2000-Dec 2016, wider than mat after this trim, so align by
  # name rather than assuming a shared positional range).
  keep_cols <- colnames(mat) <= "2014-12-01"
  mat       <- mat[, keep_cols, drop = FALSE]
  vpd_mat   <- vpd_mat_full[, colnames(mat), drop = FALSE]
  ppt_mat   <- ppt_mat_full[, colnames(mat), drop = FALSE]
  n_months  <- ncol(mat)   # 179 (Feb 2000 – Dec 2014)
  n_pixels  <- nrow(mat)   # 251,971

  # process_pixel() is defined inside this function, so its closure
  # captures this whole frame by reference — including arguments it never
  # touches. vpd_mat_full/ppt_mat_full (~409 MiB each, full Feb 2000-Dec
  # 2016 range) are superseded by the trimmed vpd_mat/ppt_mat above but
  # would otherwise still be reachable from process_pixel's environment,
  # roughly doubling what future_lapply() ships to every worker. Drop them
  # before process_pixel is defined so they're no longer reachable.
  rm(vpd_mat_full, ppt_mat_full)

  cat(sprintf("✓ Analysis window: %s to %s (%d months)\n",
              colnames(mat)[1], colnames(mat)[n_months], n_months))
  cat(sprintf("✓ Pixels: %d | WIN: %d months | NA threshold: %.0f%%\n",
              n_pixels, WIN, NA_THRESHOLD * 100))

  # Growing-season index — filter to Mar–Nov after STL (same as cluster-level
  # CSD in county_approach/src/11_1_csd_manual.R). STL still runs on all 12 months.
  col_dates  <- as.Date(colnames(mat))
  col_months <- as.integer(format(col_dates, "%m"))
  gs_idx     <- which(col_months %in% 3:11)
  gs_dates   <- col_dates[gs_idx]

  # ts() needs a numeric (year, month) start — derived from the data's own
  # first column rather than hardcoded, so this doesn't silently mislabel
  # the series if this pipeline is ever reused with data starting a
  # different month (both NDVI and NDWI currently start Feb 2000, but
  # nothing enforces that going forward).
  series_start <- c(as.integer(format(col_dates[1], "%Y")),
                     as.integer(format(col_dates[1], "%m")))

  # vpd_mat/ppt_mat are only ever read via [i, gs_idx] below — trim to growing
  # season columns once here so each parallel worker gets ~25% less data.
  vpd_mat <- vpd_mat[, gs_idx, drop = FALSE]
  ppt_mat <- ppt_mat[, gs_idx, drop = FALSE]

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

  # Placeholder for dropped pixels' rolling r1/var/skew series (below) --
  # same length as every "ok" pixel's series (length(window_idx), fixed for
  # every pixel), so rbind-ing all pixels' series into a matrix afterward
  # produces a consistent shape regardless of drop status.
  na_series <- rep(NA_real_, length(window_idx))

  # Every pixel's 7-step process (interpolate -> STL -> deseasonalize -> filter
  # to growing season -> climate correction -> rolling CSD stats -> Kendall's
  # tau) is fully independent of every other pixel, so this is embarrassingly
  # parallel. Each call returns a status string plus its result row (NA row if
  # dropped); both are tabulated together after future_lapply() finishes.
  # process_pixel() is assigned into .GlobalEnv below (see note before
  # future_lapply) so future/globals doesn't serialize a whole closure
  # frame. Everything it reads from outside its own body is therefore a
  # ".pp_"-prefixed global, not the bare local names above — bare names
  # here (mat, WIN, EDGE_NA_LIMIT, ...) would collide with the driver
  # loop's own top-level variables of the same name.
  process_pixel <- function(i) {
    .pp_p()
    y <- .pp_mat[i, ]

    # skip pixels with too many NAs — interpolation would be mostly fabricated
    if (sum(is.na(y)) >= .pp_na_limit)
      return(list(status = "dropped_high_na", row = .pp_na_row,
                  r1_series = .pp_na_series, var_series = .pp_na_series, skew_series = .pp_na_series))

    # skip pixels with a long leading/trailing NA run — na.approx(rule=2) would
    # flatline it, distorting rolling var/r1 at the flat-to-real transition
    if (edge_na_run(y) > .pp_EDGE_NA_LIMIT || edge_na_run(rev(y)) > .pp_EDGE_NA_LIMIT)
      return(list(status = "dropped_edge_na", row = .pp_na_row,
                  r1_series = .pp_na_series, var_series = .pp_na_series, skew_series = .pp_na_series))

    # linear interpolation for remaining NAs; rule=2 fills edge NAs by extension
    y <- as.numeric(zoo::na.approx(y, na.rm = FALSE, rule = 2))

    # STL — monthly series, frequency 12
    ts_y   <- ts(y, start = .pp_series_start, frequency = 12)
    decomp <- tryCatch(
      stl(ts_y, s.window = "periodic", robust = TRUE),
      error = function(e) NULL
    )
    if (is.null(decomp))
      return(list(status = "dropped_stl_fail", row = .pp_na_row,
                  r1_series = .pp_na_series, var_series = .pp_na_series, skew_series = .pp_na_series))

    # y_deseas = y - seasonal component (keeps trend + remainder). Same
    # approach as cluster-level county_approach/src/11_1_csd_manual.R.
    seasonal <- as.numeric(decomp$time.series[, "seasonal"])
    y_deseas <- y - seasonal

    # 12 monthly seasonal offsets (s.window="periodic" -> exactly one value
    # per calendar month). .pp_col_months is the calendar-month vector for
    # every column of .pp_mat (same for every pixel).
    seas_by_month <- seasonal[match(1:12, .pp_col_months)]
    names(seas_by_month) <- sprintf("seas_%02d", 1:12)

    # Filter to growing season (vpd/ppt pre-trimmed to gs_idx columns)
    y_gs   <- y_deseas[.pp_gs_idx]
    vpd_gs <- .pp_vpd_mat[i, ]
    ppt_gs <- .pp_ppt_mat[i, ]

    # Climate correction: fit lm on 2000–2011 baseline, apply to full series
    lm_data  <- data.frame(y = y_gs, vpd = vpd_gs, ppt = ppt_gs)
    clim_fit <- tryCatch(
      lm(y ~ vpd + ppt, data = lm_data[.pp_baseline_gs, ]),
      error = function(e) NULL
    )
    if (is.null(clim_fit))
      return(list(status = "dropped_clim_fail", row = .pp_na_row,
                  r1_series = .pp_na_series, var_series = .pp_na_series, skew_series = .pp_na_series))
    clim_resid <- y_gs - predict(clim_fit, newdata = lm_data)

    n <- length(clim_resid)

    # rolling stats over WIN-month sliding window (window_idx/window_dates
    # precomputed once outside this function — same for every pixel)
    r1_vec   <- rep(NA_real_, n)
    var_vec  <- rep(NA_real_, n)
    skew_vec <- rep(NA_real_, n)

    for (k in seq_along(.pp_window_idx)) {
      j           <- .pp_WIN + k - 1
      w           <- clim_resid[.pp_window_idx[[k]]]
      r1_vec[j]   <- roll_r1(w, .pp_window_dates[[k]])
      var_vec[j]  <- roll_var(w)
      skew_vec[j] <- roll_skew(w)
    }

    valid <- which(!is.na(r1_vec))
    if (length(valid) < 10)
      return(list(status = "dropped_insufficient_valid", row = .pp_na_row,
                  r1_series = .pp_na_series, var_series = .pp_na_series, skew_series = .pp_na_series))

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
      ),
      # Rolling r1/var/skew series (length(valid) == length(window_idx) == 82
      # for every "ok" pixel, since clim_resid has no internal NAs -- these
      # are the actual month-by-month values Kendall's tau above collapses
      # into a single number. Saved so a downstream timeline chart can show
      # the real trajectory instead of just the summary statistic.
      r1_series   = unname(round(r1_vec[valid],   4)),
      var_series  = unname(round(var_vec[valid],  4)),
      skew_series = unname(round(skew_vec[valid], 4))
    )
  }

  # future/globals only skips serializing a closure environment when that
  # environment is .GlobalEnv itself -- any other environment still gets
  # shipped to workers as an opaque blob on top of its own bindings being
  # separately identified by name, doubling the payload. So process_pixel()
  # must live in .GlobalEnv.
  #
  # Its ".pp_"-prefixed names (process_pixel's body reads only those) avoid
  # colliding with this driver loop's own top-level variables of the same
  # bare name (mat, WIN, EDGE_NA_LIMIT, ...).
  needed <- list(.pp_mat = mat, .pp_vpd_mat = vpd_mat, .pp_ppt_mat = ppt_mat,
                 .pp_gs_idx = gs_idx, .pp_window_idx = window_idx,
                 .pp_window_dates = window_dates, .pp_na_row = na_row,
                 .pp_na_series = na_series,
                 .pp_WIN = WIN, .pp_EDGE_NA_LIMIT = EDGE_NA_LIMIT,
                 .pp_na_limit = na_limit, .pp_col_months = col_months,
                 .pp_baseline_gs = baseline_gs, .pp_series_start = series_start)
  list2env(needed, envir = globalenv())
  environment(process_pixel) <- globalenv()
  # process_pixel itself stays local to this frame (only its environment
  # attribute points at .GlobalEnv) — future_lapply() receives it by value,
  # so it never needs to be a named global too.
  on.exit(rm(list = c(names(needed), ".pp_p"), envir = globalenv()),
          add = TRUE)

  cat(sprintf("\nRunning STL + climate correction + CSD for %d pixels across %d workers (%s)...\n",
              n_pixels, N_WORKERS, var_label))
  t_start <- proc.time()["elapsed"]

  with_progress({
    assign(".pp_p", progressor(along = seq_len(n_pixels)), envir = globalenv())
    # future.seed = NULL disables future's RNG-misuse check. process_pixel()
    # itself never touches R's RNG, so the "UNRELIABLE VALUE" warning some
    # runs produce is a false positive from the multisession/progressr
    # worker machinery, not this function.
    pixel_results <- future_lapply(seq_len(n_pixels), process_pixel,
                                    future.seed = NULL,
                                    future.packages = "zoo")
  })

  elapsed_total <- (proc.time()["elapsed"] - t_start) / 60

  pixel_status <- vapply(pixel_results, `[[`, character(1), "status")
  result_rows  <- do.call(rbind, lapply(pixel_results, `[[`, "row"))

  results <- data.frame(
    pixel_id = pixel_meta$pixel_id,
    lon      = pixel_meta$lon,
    lat      = pixel_meta$lat,
    result_rows
  )

  cat(sprintf("\n✓ %s completed in %.1f minutes\n", var_label, elapsed_total))
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

  write_csv(results, out_file)
  cat(sprintf("✓ Saved: %s\n", out_file))

  # Rolling r1/var/skew series, one row per pixel (same order as results
  # above) and one column per rolling-window end date -- the month-by-month
  # trajectory that Kendall's tau collapses into tau_r1/tau_var/tau_skew.
  # Column names are the window end dates so downstream scripts (e.g. a
  # timeline chart) can read dates directly off the matrix.
  window_end_dates <- as.character(vapply(window_dates, function(d) as.character(tail(d, 1)), character(1)))

  r1_matrix   <- do.call(rbind, lapply(pixel_results, `[[`, "r1_series"))
  var_matrix  <- do.call(rbind, lapply(pixel_results, `[[`, "var_series"))
  skew_matrix <- do.call(rbind, lapply(pixel_results, `[[`, "skew_series"))
  colnames(r1_matrix) <- colnames(var_matrix) <- colnames(skew_matrix) <- window_end_dates
  rownames(r1_matrix) <- rownames(var_matrix) <- rownames(skew_matrix) <- NULL

  saveRDS(r1_matrix, rolling_r1_file)
  cat(sprintf("✓ Saved: %s\n", rolling_r1_file))
  saveRDS(var_matrix, rolling_var_file)
  cat(sprintf("✓ Saved: %s\n", rolling_var_file))
  saveRDS(skew_matrix, rolling_skew_file)
  cat(sprintf("✓ Saved: %s\n", rolling_skew_file))
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

for (idx in INDICES) {
  cfg <- index_config[[idx]]
  cat(sprintf("\nLoading %s matrix (%s)...\n", cfg$label, cfg$matrix_path))
  mat <- readRDS(cfg$matrix_path)
  run_pixel_csd(mat, cfg$label, cfg$out_file,
                cfg$rolling_r1_file, cfg$rolling_var_file, cfg$rolling_skew_file,
                vpd_mat_full, ppt_mat_full, pixel_meta,
                WIN, NA_THRESHOLD, EDGE_NA_LIMIT, N_WORKERS)
  rm(mat)
  gc()
}

future::plan(future::sequential)  # shut down worker processes
