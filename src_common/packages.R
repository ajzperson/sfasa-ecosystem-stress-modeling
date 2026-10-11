# Installs the R packages this project uses.

install.packages(c(
  # Data manipulation & visualization
  "tidyverse",     # ggplot2, dplyr, tidyr, etc.
  "lubridate",     # date/time handling
  "ggplot2",       # advanced plotting

  # Geospatial & climate data
  "terra",         # raster data (NDVI, NDWI, VPD, PRISM, MTBS, SWE)
  "sf",            # spatial features
  "prism",         # PRISM climate data download
  "remotes",       # needed to install packages from GitHub

  # NA interpolation
  "zoo",           # na.approx()

  # Parallelization
  "future",        # parallel backend (plan())
  "future.apply",  # parallel lapply (future_lapply)
  "progressr",     # progress bars across future workers

  # Clustering
  "cluster",       # silhouette(), alongside kmeans()

  # Charts
  "patchwork",     # combine multiple ggplots into one figure
  "cowplot",       # draw_plot()/draw_line() composite figures
  "ggtext",        # markdown-formatted plot text (element_markdown, geom_textbox)

  # Testing
  "testthat"       # unit tests (pixel_approach/tests/)
))

# climateR is not on CRAN — install from GitHub:
remotes::install_github("mikejohnson51/climateR")

# Not currently used by any script
install.packages(c(
  "ncdf4",         # NetCDF file handling
  "tseries",       # time-series statistics
  "forecast",      # ARIMA, forecasting
  "car",           # VIF for multicollinearity checks
  "torch",         # deep neural networks
  "caret",         # ML utilities & preprocessing
  "multcompView",  # Tukey HSD letters
  "ggsignif",      # significance brackets
  "tigris",        # US Census county boundaries
  "elevatr"        # elevation extraction (USGS 3DEP)
))
remotes::install_github("earlywarningtoolbox/earlywarnings-R")
