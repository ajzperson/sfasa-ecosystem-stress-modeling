# Runs every test-*.R file directly in pixel_approach/tests/ -- does NOT
# recurse into tests/deprecated/, so tests for deprecated scripts (04, 02d)
# must be run individually, e.g. testthat::test_file("pixel_approach/tests/deprecated/test-04-actual-drop.R").
testthat::test_dir("pixel_approach/tests")
