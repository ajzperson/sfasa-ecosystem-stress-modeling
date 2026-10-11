library(tidyverse)
library(zoo)
library(future.apply)
library(progressr)

# Per-pixel CSD trends: STL-deseasonalize, climate-correct (VPD + PPT), then
# Kendall's tau of rolling r1/variance/skewness over the growing season.
#
# Output (per index): pixel_approach/data/processed/<index>_pixel_tau.csv
#                     pixel_approach/data/processed/<index>_rolling_{r1,var,skew}.rds

# Show the progress bar under Rscript too
options(progressr.enable = TRUE)
handlers("txtprogressbar")

# Parameters (set before source() to override)
if (!exists("INDICES"))            INDICES            <- c("ndvi", "ndwi")
if (!exists("WIN"))                WIN                <- 54     # rolling window (growing-season months)
if (!exists("NA_THRESHOLD"))       NA_THRESHOLD       <- 0.30   # skip pixel if NA fraction >= this
if (!exists("EDGE_NA_LIMIT"))      EDGE_NA_LIMIT      <- 6      # skip pixel if leading/trailing NA run > this
if (!exists("N_WORKERS"))          N_WORKERS          <- 5
if (!exists("FUTURE_MAX_SIZE_MB")) FUTURE_MAX_SIZE_MB <- 2048   # above future's 500 MiB default; workers need the matrices
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

cat("Loading pixel metadata and climate inputs...\n")
pixel_meta   <- read_csv("pixel_approach/data/processed/pixel_metadata.csv",
                         show_col_types = FALSE)
vpd_mat_full <- readRDS("pixel_approach/data/processed/vpd_matrix.rds")
ppt_mat_full <- readRDS("pixel_approach/data/processed/ppt_matrix.rds")

# Lag-1 autocorrelation over 1-month-apart pairs only (skips the Nov->Mar gap)
roll_r1 <- function(x, d) {
  xbar <- mean(x)
  midx <- as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m"))
  keep <- diff(midx) == 1L
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

# Length of the leading NA run (na.approx would flatline a long one)
edge_na_run <- function(x) {
  na_idx <- is.na(x)
  if (!na_idx[1]) return(0L)
  rle(na_idx)$lengths[1]
}

# Runs the full per-pixel pipeline for one index matrix
run_pixel_csd <- function(mat, var_label, out_file,
                           rolling_r1_file, rolling_var_file, rolling_skew_file,
                           vpd_mat_full, ppt_mat_full, pixel_meta,
                           WIN, NA_THRESHOLD, EDGE_NA_LIMIT, N_WORKERS) {

  cat(sprintf("\n=== %s ===\n", var_label))

  # Trim to 2000-2014 and align vpd/ppt columns by name
  keep_cols <- colnames(mat) <= "2014-12-01"
  mat       <- mat[, keep_cols, drop = FALSE]
  vpd_mat   <- vpd_mat_full[, colnames(mat), drop = FALSE]
  ppt_mat   <- ppt_mat_full[, colnames(mat), drop = FALSE]
  n_months  <- ncol(mat)
  n_pixels  <- nrow(mat)

  # Free the untrimmed matrices so they aren't shipped to workers
  rm(vpd_mat_full, ppt_mat_full)

  cat(sprintf("✓ Analysis window: %s to %s (%d months)\n",
              colnames(mat)[1], colnames(mat)[n_months], n_months))
  cat(sprintf("✓ Pixels: %d | WIN: %d months | NA threshold: %.0f%%\n",
              n_pixels, WIN, NA_THRESHOLD * 100))

  # Growing season (Mar-Nov), applied after STL
  col_dates  <- as.Date(colnames(mat))
  col_months <- as.integer(format(col_dates, "%m"))
  gs_idx     <- which(col_months %in% 3:11)
  gs_dates   <- col_dates[gs_idx]

  series_start <- c(as.integer(format(col_dates[1], "%Y")),
                     as.integer(format(col_dates[1], "%m")))

  vpd_mat <- vpd_mat[, gs_idx, drop = FALSE]
  ppt_mat <- ppt_mat[, gs_idx, drop = FALSE]

  # Climate model is fit on the pre-drought 2000-2011 growing seasons
  baseline_gs <- gs_dates <= as.Date("2011-12-31")

  na_limit <- NA_THRESHOLD * n_months
  n_gs     <- length(gs_idx)

  # Rolling windows are the same for every pixel, so build them once
  window_idx   <- lapply(WIN:n_gs, function(j) (j - WIN + 1):j)
  window_dates <- lapply(window_idx, function(idx) gs_dates[idx])

  # Seasonal offsets and climate coefficients are saved for reuse in 12
  na_row <- c(tau_r1 = NA_real_, p_r1 = NA_real_, tau_var = NA_real_,
              p_var = NA_real_, tau_skew = NA_real_, p_skew = NA_real_,
              setNames(rep(NA_real_, 12), sprintf("seas_%02d", 1:12)),
              b0 = NA_real_, b1 = NA_real_, b2 = NA_real_)

  # Rolling series placeholder for dropped pixels
  na_series <- rep(NA_real_, length(window_idx))

  # Reads its inputs from ".pp_"-prefixed globals (see list2env below)
  process_pixel <- function(i) {
    .pp_p()
    y <- .pp_mat[i, ]

    if (sum(is.na(y)) >= .pp_na_limit)
      return(list(status = "dropped_high_na", row = .pp_na_row,
                  r1_series = .pp_na_series, var_series = .pp_na_series, skew_series = .pp_na_series))

    if (edge_na_run(y) > .pp_EDGE_NA_LIMIT || edge_na_run(rev(y)) > .pp_EDGE_NA_LIMIT)
      return(list(status = "dropped_edge_na", row = .pp_na_row,
                  r1_series = .pp_na_series, var_series = .pp_na_series, skew_series = .pp_na_series))

    y <- as.numeric(zoo::na.approx(y, na.rm = FALSE, rule = 2))

    ts_y   <- ts(y, start = .pp_series_start, frequency = 12)
    decomp <- tryCatch(
      stl(ts_y, s.window = "periodic", robust = TRUE),
      error = function(e) NULL
    )
    if (is.null(decomp))
      return(list(status = "dropped_stl_fail", row = .pp_na_row,
                  r1_series = .pp_na_series, var_series = .pp_na_series, skew_series = .pp_na_series))

    seasonal <- as.numeric(decomp$time.series[, "seasonal"])
    y_deseas <- y - seasonal

    # One seasonal offset per calendar month (s.window = "periodic")
    seas_by_month <- seasonal[match(1:12, .pp_col_months)]
    names(seas_by_month) <- sprintf("seas_%02d", 1:12)

    y_gs  <- y_deseas[.pp_gs_idx]
    vpd_gs <- .pp_vpd_mat[i, ]
    ppt_gs <- .pp_ppt_mat[i, ]

    # Fit on 2000-2011, predict 2000-2014
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

    r1_vec  <- rep(NA_real_, n)
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

    # unname() keeps column names like "tau_r1" instead of "tau_r1.tau"
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
      # Full rolling series, saved for the timeline chart in 12
      r1_series   = unname(round(r1_vec[valid],   4)),
      var_series  = unname(round(var_vec[valid],  4)),
      skew_series = unname(round(skew_vec[valid], 4))
    )
  }

  # process_pixel lives in .GlobalEnv so future doesn't also ship this whole frame
  needed <- list(.pp_mat = mat, .pp_vpd_mat = vpd_mat, .pp_ppt_mat = ppt_mat,
                 .pp_gs_idx = gs_idx, .pp_window_idx = window_idx,
                 .pp_window_dates = window_dates, .pp_na_row = na_row,
                 .pp_na_series = na_series,
                 .pp_WIN = WIN, .pp_EDGE_NA_LIMIT = EDGE_NA_LIMIT,
                 .pp_na_limit = na_limit, .pp_col_months = col_months,
                 .pp_baseline_gs = baseline_gs, .pp_series_start = series_start)
  list2env(needed, envir = globalenv())
  environment(process_pixel) <- globalenv()
  on.exit(rm(list = c(names(needed), ".pp_p"), envir = globalenv()),
          add = TRUE)

  cat(sprintf("\nRunning STL + climate correction + CSD for %d pixels across %d workers (%s)...\n",
              n_pixels, N_WORKERS, var_label))
  t_start <- proc.time()["elapsed"]

  with_progress({
    assign(".pp_p", progressor(along = seq_len(n_pixels)), envir = globalenv())
    # No RNG is used, so silence future's false "UNRELIABLE VALUE" warning
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

  # One row per pixel, one column per window end date
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

# Single-threaded BLAS per worker to avoid oversubscribing cores (set before plan())
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

future::plan(future::sequential)
