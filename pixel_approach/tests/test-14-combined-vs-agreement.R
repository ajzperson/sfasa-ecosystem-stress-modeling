library(testthat)
library(tidyverse)
library(patchwork)
library(ggtext)

# Tests for 14_combined_vs_agreement.R, run on synthetic inputs (charts skipped).
exprs    <- parse("../src/14_combined_vs_agreement.R")
src_text <- vapply(exprs, function(e) paste(deparse(e), collapse = "\n"), character(1))
skip_lines <- grepl(
  paste("^combined <-", "^ndvi <-", "^ndwi <-", "^beetle <-", "^fire <-",
        "^write_csv\\(", "^plot_df <-", "^group_n <-", "^group_x_labels <-",
        "^p <- ggplot\\(", "^build_pie <-", "^n_total <-", "^pie_alone <-",
        "^pie_combined <-", "^build_note <-", "^note_left <-",
        "^note_right <-", "^p_full <-", "^ggsave\\(",
        sep = "|"),
  src_text
)

# Pixels 1-10 Both, 11-20 NDVI+NDWI_4 only, 21-30 NDVI_3 ∩ NDWI_3 only, 31-40 Neither
combined <- tibble(pixel_id = 1:40,
                    combined_cluster = c(rep(4, 20), rep(2, 10), rep(1, 10)))
ndvi <- tibble(pixel_id = 1:40,
               ndvi_cluster = c(rep(3, 10), rep(1, 10), rep(3, 10), rep(1, 10)))
ndwi <- tibble(pixel_id = 1:40,
               ndwi_cluster = c(rep(3, 10), rep(1, 10), rep(3, 10), rep(1, 10)))

fire_burned <- c(rep(TRUE, 7), rep(FALSE, 3),    # Both: 7/10
                 rep(TRUE, 4), rep(FALSE, 6),    # NDVI+NDWI_4 only: 4/10
                 rep(TRUE, 3), rep(FALSE, 7),    # NDVI_3 ∩ NDWI_3 only: 3/10
                 rep(TRUE, 1), rep(FALSE, 9))    # Neither: 1/10
beetle_mortality <- c(rep(TRUE, 5), rep(FALSE, 5),
                       rep(TRUE, 3), rep(FALSE, 7),
                       rep(TRUE, 2), rep(FALSE, 8),
                       rep(TRUE, 1), rep(FALSE, 9))

beetle <- tibble(pixel_id = 1:40,
                  tpa_2015 = ifelse(beetle_mortality, 1, 0), tpa_2016 = 0)
fire <- tibble(pixel_id = 1:40,
               burned_2015 = fire_burned, burned_2016 = FALSE)

for (i in seq_along(exprs)) {
  if (!skip_lines[i]) eval(exprs[[i]])
}

test_that("4-way overlap group classification matches the constructed fixture", {
  expect_setequal(d$pixel_id[d$group == "Both"], 1:10)
  expect_setequal(d$pixel_id[d$group == "\"NDVI+NDWI_4\" only"], 11:20)
  expect_setequal(d$pixel_id[d$group == "\"NDVI_3 ∩ NDWI_3\" only"], 21:30)
  expect_setequal(d$pixel_id[d$group == "Neither"], 31:40)
})

test_that("Jaccard index matches hand-computed overlap / union", {
  expect_equal(n_a, 20)
  expect_equal(n_b, 20)
  expect_equal(n_both, 10)
  expect_equal(n_union, 30)
  expect_equal(jaccard, 10 / 30, tolerance = 1e-6)
})

test_that("summary_tbl computes the correct fire/beetle rate per group", {
  s <- summary_tbl
  expect_equal(s$pct_fire[s$group == "Both"], 70)
  expect_equal(s$pct_fire[s$group == "Neither"], 10)
  expect_equal(s$pct_beetle[s$group == "Both"], 50)
  expect_equal(s$pct_beetle[s$group == "Neither"], 10)
})

test_that("relative risk (Both vs Neither) matches hand-computed ratios", {
  expect_equal(rr_fire, 7, tolerance = 1e-6)
  expect_equal(rr_beetle, 5, tolerance = 1e-6)
})

test_that("wrap_def formats a bold markdown term/definition pair", {
  expect_equal(wrap_def("NDVI_3", "NDVI Most Stressed Cluster"),
               "**NDVI_3 = NDVI Most Stressed Cluster**")
})
