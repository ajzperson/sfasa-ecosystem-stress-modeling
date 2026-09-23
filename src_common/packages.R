# Installs R package dependencies for this project's climate downloads,
# geospatial processing, CSD time-series analysis, parallelization,
# clustering, charts, and testing.

install.packages(c(
  # Data manipulation & visualization
  "tidyverse",     # ggplot2, dplyr, tidyr, etc.
  "lubridate",     # date/time handling
  "ggplot2",       # advanced plotting

  # Geospatial & climate data
  "terra",         # raster data (NDVI, NDWI, VPD, PRISM, MTBS, SWE)
  "sf",            # spatial features
  "tigris",        # US Census county boundaries
  "prism",         # PRISM climate data download
  "elevatr",       # elevation extraction (USGS 3DEP)
  # climateR — install from GitHub (not on CRAN), see below.
  # AOI is not needed — sf objects work directly.

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

# earlywarnings is not on CRAN — install from GitHub:
remotes::install_github("earlywarningtoolbox/earlywarnings-R")

# --- Not currently used by any script in this repo — kept in case of future work ---
install.packages(c(
  "ncdf4",         # NetCDF file handling
  "tseries",       # time-series statistics
  "forecast",      # ARIMA, forecasting
  "car",           # VIF (variance inflation factor) for multicollinearity checks
  "torch",         # deep neural networks
  "caret",         # ML utilities & preprocessing
  "multcompView",  # Tukey HSD letters -- loaded, but never actually called
  "ggsignif"       # significance brackets -- dead code in 06_clustering.R
))
