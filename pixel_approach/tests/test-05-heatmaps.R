library(testthat)
library(tidyverse)
library(sf)
library(patchwork)

# Load only plot_heatmap()/save_combined() and sn_boundary from
# 05_heatmaps.R, skipping the real d_ndvi/d_ndwi CSV loads (large, and
# unnecessary -- these tests use small synthetic data instead) and
# stopping before the four real save_combined(...) calls at the bottom.
exprs    <- parse("../src/05_heatmaps.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_data_loads <- grepl("ndvi_pixel_tau\\.csv|ndwi_pixel_tau\\.csv", src_text)
stop_at  <- which(grepl("^save_combined\\(", src_text))[1]
stopifnot(!is.na(stop_at))

old_wd <- getwd()
setwd("../..")
for (i in seq_len(stop_at - 1)) {
  if (!skip_data_loads[i]) eval(exprs[[i]], envir = .GlobalEnv)
}
setwd(old_wd)

make_synth_data <- function(n = 20) {
  data.frame(
    lon      = seq(-120, -119, length.out = n),
    lat      = seq(37, 38, length.out = n),
    tau_r1   = seq(-1, 1, length.out = n),
    tau_var  = seq(-0.5, 0.5, length.out = n),
    Actual_Drop_15_16 = seq(0, 0.5, length.out = n)
  )
}

test_that("plot_heatmap returns a ggplot object", {
  d <- make_synth_data()
  p <- plot_heatmap(d, "tau_r1", "TEST", "#2166AC", "#D73027", limits = c(-1, 1))
  expect_s3_class(p, "ggplot")
})

test_that("plot_heatmap uses the fixed limits when provided", {
  d <- make_synth_data()
  p <- plot_heatmap(d, "tau_r1", "TEST", "#2166AC", "#D73027", limits = c(-1, 1))
  scale_limits <- p$scales$scales[[1]]$limits
  expect_equal(scale_limits, c(-1, 1))
})

test_that("plot_heatmap auto-scales (NULL limits) when none are given", {
  d <- make_synth_data()
  p <- plot_heatmap(d, "Actual_Drop_15_16", "TEST", "#2166AC", "#D73027")
  scale_limits <- p$scales$scales[[1]]$limits
  expect_null(scale_limits)
})

test_that("plot_heatmap's legend name defaults to the column name", {
  d <- make_synth_data()
  p <- plot_heatmap(d, "tau_var", "TEST", "#2166AC", "#D73027", limits = c(-1, 1))
  expect_equal(p$scales$scales[[1]]$name, "tau_var")
})

test_that("plot_heatmap's legend name can be overridden", {
  d <- make_synth_data()
  p <- plot_heatmap(d, "Actual_Drop_15_16", "TEST", "#2166AC", "#D73027",
                     legend_name = "Custom Name")
  expect_equal(p$scales$scales[[1]]$name, "Custom Name")
})

test_that("plot_heatmap's title matches panel_title", {
  d <- make_synth_data()
  p <- plot_heatmap(d, "tau_r1", "MY PANEL TITLE", "#2166AC", "#D73027", limits = c(-1, 1))
  expect_equal(p$labels$title, "MY PANEL TITLE")
})

test_that("plot_heatmap correctly maps color to the specified column, not another one", {
  d <- make_synth_data()
  p_r1  <- plot_heatmap(d, "tau_r1",  "TEST", "#2166AC", "#D73027", limits = c(-1, 1))
  p_var <- plot_heatmap(d, "tau_var", "TEST", "#2166AC", "#D73027", limits = c(-1, 1))

  built_r1  <- ggplot_build(p_r1)$data[[1]]$colour
  built_var <- ggplot_build(p_var)$data[[1]]$colour
  # Different source columns (tau_r1 vs tau_var have different value
  # ranges/ordering) must produce different rendered colors.
  expect_false(identical(built_r1, built_var))
})

test_that("save_combined writes a file and works with both collect_guides settings", {
  d <- make_synth_data()
  p1 <- plot_heatmap(d, "tau_r1", "A", "#2166AC", "#D73027", limits = c(-1, 1))
  p2 <- plot_heatmap(d, "tau_r1", "B", "#2166AC", "#D73027", limits = c(-1, 1))

  out_collected <- tempfile(fileext = ".pdf")
  save_combined(p1, p2, "TEST TITLE", out_collected, collect_guides = TRUE)
  expect_true(file.exists(out_collected))

  out_kept <- tempfile(fileext = ".pdf")
  save_combined(p1, p2, "TEST TITLE", out_kept, collect_guides = FALSE)
  expect_true(file.exists(out_kept))
})
