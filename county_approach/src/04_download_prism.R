library(prism)
library(tidyverse)

cat("Downloading PRISM climate data...\n")

# Set download directory
prism_set_dl_dir("data_common/raw/prism/prismz")

# Download monthly precipitation (2000-2026)
cat("Downloading precipitation data (2000-2026)...\n")
get_prism_monthlys(type = "ppt", years = 2000:2026, keepZip = FALSE)

# Download monthly mean temperature (2000-2026)
cat("Downloading mean temperature data (2000-2026)...\n")
get_prism_monthlys(type = "tmean", years = 2000:2026, keepZip = FALSE)

# Download monthly minimum temperature (2000-2026)
cat("Downloading minimum temperature data (2000-2026)...\n")
get_prism_monthlys(type = "tmin", years = 2000:2026, keepZip = FALSE)

# Download monthly maximum temperature (2000-2026)
cat("Downloading maximum temperature data (2000-2026)...\n")
get_prism_monthlys(type = "tmax", years = 2000:2026, keepZip = FALSE)

cat("\n✓ PRISM data downloaded to data_common/raw/prism/prismz/\n")
cat("  Next: Run 05_process_prism.R to extract Sierra Nevada values\n")
