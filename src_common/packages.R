install.packages(c(
  # Data manipulation & visualization
  "tidyverse",     # ggplot2, dplyr, tidyr, etc.
  "lubridate",     # date/time handling
  "ggplot2",       # advanced plotting

  # Geospatial & climate data
  "terra",         # raster data (PRISM + NDVI processing)
  "sf",            # spatial features
  "tigris",        # US Census county boundaries
  "ncdf4",         # NetCDF file handling (NDVI)
  "prism",         # PRISM climate data download
  # climateR + AOI — install from GitHub (not on CRAN):
  # remotes::install_github("mikejohnson51/climateR")
  # remotes::install_github("mikejohnson51/AOI")  # optional, sf objects work directly

  "remotes",       # needed to install packages from GitHub

  # Time-series analysis
  "tseries",       # time-series statistics
  "forecast",      # ARIMA, forecasting
  "zoo",           # rolling windows, time-series

  # Regression diagnostics
  "car",           # VIF (variance inflation factor) for multicollinearity checks
  "patchwork",     # combine multiple ggplots into one figure
  "multcompView",  # compact letter display (a/b/c) for Tukey HSD post-hoc
  "ggsignif",      # significance brackets on ggplot boxplots

  # Machine learning
  "torch",         # deep neural networks
  "caret"          # ML utilities & preprocessing
))

# earlywarnings is not on CRAN — install from GitHub:
remotes::install_github("earlywarningtoolbox/earlywarnings-R")
