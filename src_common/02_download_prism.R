library(prism)
library(tidyverse)

# Downloads PRISM monthly climate data (precipitation, mean/min/max
# temperature), 2000-2026, CONUS-wide (PRISM has no AOI-clip option).
#
# Output: data_common/raw/prism/prismz/

cat("Downloading PRISM climate data...\n")

prism_set_dl_dir("data_common/raw/prism/prismz")

cat("Downloading precipitation data (2000-2026)...\n")
get_prism_monthlys(type = "ppt", years = 2000:2026, keepZip = FALSE)

cat("Downloading mean temperature data (2000-2026)...\n")
get_prism_monthlys(type = "tmean", years = 2000:2026, keepZip = FALSE)

cat("Downloading minimum temperature data (2000-2026)...\n")
get_prism_monthlys(type = "tmin", years = 2000:2026, keepZip = FALSE)

cat("Downloading maximum temperature data (2000-2026)...\n")
get_prism_monthlys(type = "tmax", years = 2000:2026, keepZip = FALSE)

cat("\n✓ PRISM data downloaded to data_common/raw/prism/prismz/\n")
