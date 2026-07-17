library(testthat)
library(terra)

# 02d_pixel_tree_mortality.R is a top-level script, not a set of functions,
# so these tests replicate its key logic on small synthetic rasters rather
# than sourcing the script against the real multi-GB LCMAP files. The
# buffer/batching extraction mechanics are identical to 02c_pixel_fire.R
# (already covered in test-02c-fire.R) -- this file focuses on what's
# distinct here: the LCACHG tree-cover-loss code classification.

TREE_LOSS_CODES <- c(41, 42, 43, 45, 46, 47, 48)  # same as the real script

test_that("LCACHG reclassification flags only prior-class-4 (tree) -> other transitions", {
  # 1-8: single-digit, no-change codes (e.g. 4 = "was and still is tree" --
  # NOT a loss). 41-48: two-digit (prior)(current) transition codes; only
  # {41,42,43,45,46,47,48} mean "was tree, now something else" (44 doesn't
  # appear as a real value since staying tree is encoded as the single
  # digit "4", not "44"). Also include a transition INTO tree (14 = prior
  # class 1, now tree) and a same-family-but-different-prior-class
  # transition (52) to confirm those are correctly excluded too.
  codes <- c(1, 4, 8, 41, 42, 43, 45, 46, 47, 48, 14, 52)
  expected_loss <- c(0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0)

  r <- rast(nrows = 1, ncols = length(codes), xmin = 0, xmax = length(codes), ymin = 0, ymax = 1,
            crs = "EPSG:32610")
  values(r) <- codes
  r_bin <- app(r, function(v) as.numeric(v %in% TREE_LOSS_CODES))

  expect_equal(values(r_bin)[, 1], expected_loss)
})

test_that("NA from no-overlap extraction is treated as FALSE (no loss), matching 02d's edge-case handling", {
  pixel_meta <- data.frame(
    pixel_id = 1:3,
    lost_tree_cover_2015 = c(TRUE, NA, FALSE),
    lost_tree_cover_2016 = c(NA, FALSE, TRUE)
  )
  year_cols <- c("lost_tree_cover_2015", "lost_tree_cover_2016")
  for (col in year_cols) pixel_meta[[col]][is.na(pixel_meta[[col]])] <- FALSE

  expect_equal(pixel_meta$lost_tree_cover_2015, c(TRUE, FALSE, FALSE))
  expect_equal(pixel_meta$lost_tree_cover_2016, c(FALSE, FALSE, TRUE))
  expect_false(anyNA(pixel_meta[year_cols]))
})
