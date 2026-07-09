# SFASA Project
**Using Time-Series Statistics ~~and Deep Neural Networks~~ to Predict Ecosystem Collapse**

Case Study: The 2012–2016 Sierra Nevada Canopy Collapse

---

## Deliverables
- Written research paper
- Working R code / notebooks
- Visualizations

---

## Data Source

1. **EPA Boundary** — Sierra Nevada ecoregion polygon (spatial framework for all extractions)
2. **Kaggle** — California Wildfire Incidents (validation data)
3. **VPD** — Vapor Pressure Deficit (TerraClimate)
4. **PRISM** — Temperature + Precipitation
5. **NDVI** — Normalized Difference Vegetation Index (NASA MODIS)
6. **NDWI** — Normalized Difference Water Index / vegetation water content (NASA MODIS) *(not currently used)*
7. **SWE** — Snow Water Equivalent (TerraClimate) *(downloaded, not currently used)*
8. **USFS IDS** — Forest Damage / Tree Mortality (Insect and Disease Survey)

All raw downloads are shared across approaches and live under `data_common/raw/`; the extracted boundary/cluster polygons live under `docs_common/`.

### Raw Data Inventory

| Data | File | Download | Spatial scope |
|------|------|----------|---------------|
| EPA Boundary | `data_common/raw/ecoregion/us_eco_l3.zip` (~28 MB) | [EPA Level III and IV Ecoregions of the Continental United States](https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states)<br>Download "US Level III Ecoregions shapefile without state boundaries" (28 MB) | Full CONUS |
| Kaggle | `data_common/raw/kaggle/California_Fire_Incidents_2013-2019_raw.csv` | [Kaggle](https://www.kaggle.com/datasets/ananthu017/california-wildfire-incidents-20132020) | California statewide |
| VPD | `data_common/raw/vpd/TerraClimate_vpd_2000-2025.nc` (~36 MB) | `climateR::getTerraClim(varname = "vpd", startDate = "2000-01-01", endDate = "2025-12-31")` | 2000-01-01 to 2025-12-01 · 312 months · ~4 km |
| PRISM | `data_common/raw/prism/prismz.zip` (Google Drive, ~1,004.8 MB) | `get_prism_monthlys(type = "ppt/tmean/tmin/tmax", years = 2000:2026)` | 2000-01-01 to 2026-05-01 · 317 months · Full CONUS · ~4 km |
| NDVI | `data_common/raw/ndvi/MOD13A1.061_500m_aid0001.nc` (Google Drive, ~906 MB) | [NASA AppEEARS](https://appeears.earthdatacloud.nasa.gov/) (NASA Earthdata account required)<br>Area Sample request<br>Product: Terra MODIS Vegetation Indices (NDVI & EVI), MOD13A1.061, 500m, 16-day (2000-02-18 to Present)<br>Layer: `_500m_16_days_NDVI`<br>Upload polygon: `docs_common/sierra_nevada_boundary_AppEEARS.geojson`<br>Date range: 2000-01-01 to 2026-05-31<br>Output: NetCDF-4<br>Projection: MODIS Sinusoidal | 2000-02-18 to 2026-05-25 · 605 16-day time steps · clipped to `sierra_nevada_boundary_AppEEARS.geojson` (AppEEARS upload polygon) · 500m |
| NDWI | `data_common/raw/ndwi/MOD09A1.061_sur_refl_b02_b06.nc` (pending download) | [NASA AppEEARS](https://appeears.earthdatacloud.nasa.gov/) (NASA Earthdata account required)<br>Area Sample request<br>Product: Terra MODIS Surface Reflectance Bands 1-7, MOD09A1.061, 500m, 8-day (2000-02-24 to Present)<br>Layers: `sur_refl_b02` (NIR, 841–876nm) + `sur_refl_b06` (SWIR, 1628–1652nm)<br>Upload polygon: `docs_common/sierra_nevada_boundary.geojson`<br>Date range: 2000-01-01 to 2026-05-31<br>Output: NetCDF-4<br>Projection: MODIS Sinusoidal<br>NDWI = (B02 − B06) / (B02 + B06) computed locally after download | 2000-02-24 to Present · 8-day composites · clipped to `sierra_nevada_boundary.geojson` · 500m |
| SWE *(unused)* | `data_common/raw/swe/TerraClimate_swe_2000-2025.nc` | `climateR::getTerraClim(varname = "swe", startDate = "2000-01-01", endDate = "2025-12-31")` | 2000-01-01 to 2025-12-01 · 312 months · ~4 km |
| Forest Damage | `data_common/raw/forest_damage/CONUS_Region5_AllYears.gdb.zip` (Google Drive, ~166 MB) | [USDA Forest Service Detection Surveys](https://www.fs.usda.gov/science-technology/data-tools-products/fhp-mapping-reporting/detection-surveys)<br>Scroll to "R5 IDS DATA DOWNLOADS" section<br>Download "R5 CA All Years IDS Data (ZIP, 159 MB)"<br>[Direct URL](https://www.fs.usda.gov/foresthealth/docs/IDS_Data_for_Download/CONUS_Region5_AllYears.gdb.zip)<br>No login required | 1997–2025 · annual · Region 5 California · 348,829 damage polygons |

---

### EPA Sierra Nevada Boundary (`sierra_nevada_boundary.geojson`)

- **Boundary:** EPA Level III Ecoregion — Sierra Nevada (`US_L3CODE = "5"`, `US_L3NAME = "Sierra Nevada"`). Area: ~53,102 km². CRS: OGC:1.3:CRS84 (WGS84). Delineated by ecological characteristics (vegetation, geology, climate) rather than administrative lines — widely cited in peer-reviewed ecology literature.
- Filtered to `US_L3CODE == "5"`, saved as `docs_common/sierra_nevada_boundary.geojson` — the shared spatial mask used by every downstream extraction, regardless of approach.

---

### VPD — Vapor Pressure Deficit (TerraClimate)

- **Why VPD over PDSI:** PDSI drainage divisions don't align cleanly with sub-regions of interest within the Sierra Nevada. VPD is also a more direct measure of atmospheric water demand that drove tree mortality.

#### Raw data format
The raw file is a NetCDF-4 (`.nc`) file — a scientific format for multidimensional arrays, not inspectable like a CSV. It contains a 3D array: **312 monthly time steps × latitude grid × longitude grid**. TerraClimate is a global dataset at ~4 km resolution; only the Sierra Nevada extent was downloaded by passing the boundary polygon as AOI to `getTerraClim()`. Each cell in the grid holds a VPD value in kPa for that location and month. Readable in R with `terra::rast()` or `ncdf4::nc_open()`.

### PRISM — Temperature + Precipitation

- **Raw file:** `data_common/raw/prism/prismz.zip` (~1,004.8 MB) — too large for GitHub; stored on Google Drive. See `data_common/raw/prism/PRISM_source_location.md` for the download link.
- `ppt`: precipitation (rainfall + snowfall), in millimeters
- `tmean`: mean temperature (average of daily high and low), in degrees Celsius
- `tmin`: monthly minimum temperature, in degrees Celsius
- `tmax`: monthly maximum temperature, in degrees Celsius

```r
install.packages("prism")
library(prism)
prism_set_dl_dir("data_common/raw/prism/prismz")
get_prism_monthlys(type = "ppt",   years = 2000:2026, keepZip = FALSE)
get_prism_monthlys(type = "tmean", years = 2000:2026, keepZip = FALSE)
get_prism_monthlys(type = "tmin",  years = 2000:2026, keepZip = FALSE)
get_prism_monthlys(type = "tmax",  years = 2000:2026, keepZip = FALSE)
```

#### Raw data format
PRISM downloads as a collection of `.bil` (Band Interleaved by Line) raster files — one file per variable per month, saved under `data_common/raw/prism/prismz/`. Unlike a single NetCDF file, each month is its own file. Each `.bil` file is a 2D raster (latitude × longitude) covering the contiguous United States at ~4 km resolution, holding one value per pixel for that variable and month. PRISM is US-only (unlike TerraClimate which is global), but higher quality for the US because it incorporates weather station data and terrain adjustments.

### NDVI — Normalized Difference Vegetation Index

- File (`MOD13A1.061_500m_aid0001.nc`, ~906 MB) is too large for GitHub — stored on Google Drive
- See `data_common/raw/ndvi/NDVI_source_location.md` for the download link
- Download from Google Drive and place at `data_common/raw/ndvi/MOD13A1.061_500m_aid0001.nc` before running any downstream extraction

#### Raw data format
The raw file is a NetCDF-4 (`.nc`) file with the same 3D structure as VPD: **time steps × latitude × longitude**. Unlike VPD (monthly) and PRISM (monthly), NDVI is a 16-day composite — one layer every 16 days, giving 605 time steps from 2000-02-18 to 2026-05-25. Spatial resolution is ~500m per pixel, much finer than VPD (~4km) or PRISM (~4km). AppEEARS clips the file to the submitted polygon.

### SWE — Snow Water Equivalent (TerraClimate) *(not currently used)*
- **Variable:** `swe` (Snow Water Equivalent, mm)
- **Raw file:** `data_common/raw/swe/TerraClimate_swe_2000-2025.nc`

### USFS IDS — Forest Damage / Tree Mortality
- **File:** `data_common/raw/forest_damage/CONUS_Region5_AllYears.gdb.zip` (~166 MB) — too large for GitHub; stored on Google Drive (see `data_common/raw/forest_damage/FOREST_DAMAGE_source_location.md`)
- Raw layer (`DAMAGE_AREAS_FLAT`) covers all forest damage causal agents (insects, disease, abiotic), not just bark beetle — downstream processing filters to bark-beetle-specific damage causal agent codes.
- **Why acquired:** Bark beetle infestation is a candidate mechanism behind loss of resilience — used as ground-truth validation against modeled collapse signal.

---

## Tech Stack
- **Language:** R
- **Key packages:** see `src_common/packages.R`
- **Notebooks:** R Markdown (`.Rmd`) in RStudio
