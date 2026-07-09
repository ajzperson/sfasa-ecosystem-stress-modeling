# SFASA Project
**Using Time-Series Statistics ~~and Deep Neural Networks~~ to Predict Ecosystem Collapse**

Case Study: The 2012–2016 Sierra Nevada Canopy Collapse

---

## Deliverables
- Poster (including visualizations)
- Working R code / notebooks
- Written research paper

---

## Data Source

1. **EPA Boundary** — EPA Level III Ecoregions shapefile (full continental US).<br>
   EPA = Environmental Protection Agency (the U.S. federal agency).<br>
   Used to generate the Sierra Nevada ecoregion polygon (spatial framework/mask that every raw dataset gets clipped/extracted to).
2. **Kaggle** — California Wildfire Incidents (validation data)
3. **VPD** — Vapor Pressure Deficit (TerraClimate)
4. **PRISM** — Temperature + Precipitation
5. **NDVI** — Normalized Difference Vegetation Index (NASA MODIS)
6. **NDWI** — Normalized Difference Water Index (NASA MODIS)
7. **SWE** — Snow Water Equivalent (TerraClimate) *(downloaded, not currently used)*
8. **USFS IDS** — Forest Damage / Tree Mortality (Insect and Disease Survey)<br>
   USFS IDS = US Forest Service, Insect and Disease Survey (also called the Aerial Detection Survey/ADS).
9. **Fire Perimeters** — California historical fire perimeter polygons (CAL FIRE FRAP).

All raw downloads are shared across approaches and live under `data_common/raw/`.

### Raw Data Inventory

| Data | File | Download | Spatial scope |
|------|------|----------|---------------|
| EPA Boundary | `data_common/raw/ecoregion/us_eco_l3.zip` | [EPA Level III and IV Ecoregions of the Continental United States](https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states)<br>Download "US Level III Ecoregions shapefile without state boundaries (28 MB)" | Full CONUS |
| Kaggle | `data_common/raw/kaggle/California_Fire_Incidents_2013-2019_raw.csv` | [Kaggle](https://www.kaggle.com/datasets/ananthu017/california-wildfire-incidents-20132020) | California statewide |
| VPD | `data_common/raw/vpd/TerraClimate_vpd_2000-2025.nc` | `vpd_raw <- getTerraClim(`<br>`  AOI       = sn_aoi,`<br>`  varname   = "vpd",`<br>`  startDate = "2000-01-01",`<br>`  endDate   = "2025-12-31"`<br>`)` | 2000-01-01 to 2025-12-01 · 312 months · clipped to Sierra Nevada ecoregion boundary · ~4 km |
| PRISM | `data_common/raw/prism/prismz.zip` (Google Drive, ~2,103 MB)<br>See [PRISM_source_location.md](../data_common/raw/prism/PRISM_source_location.md) | `get_prism_monthlys(type = "ppt", years = 2000:2026, keepZip = FALSE)`<br>`get_prism_monthlys(type = "tmean", years = 2000:2026, keepZip = FALSE)`<br>`get_prism_monthlys(type = "tmin", years = 2000:2026, keepZip = FALSE)`<br>`get_prism_monthlys(type = "tmax", years = 2000:2026, keepZip = FALSE)` | 2000-01-01 to 2026-05-01 · 317 months · Full CONUS · ~4 km |
| NDVI | `data_common/raw/ndvi/MOD13A1.061_500m_aid0001.nc` (Google Drive, ~1,019 MB)<br>See [NDVI_source_location.md](../data_common/raw/ndvi/NDVI_source_location.md) | [NASA AppEEARS](https://appeears.earthdatacloud.nasa.gov/) (NASA Earthdata account required)<br>Area Sample request<br>Product: Terra MODIS Vegetation Indices (NDVI & EVI), MOD13A1.061, 500m, 16-day (2000-02-18 to Present)<br>Layer: `_500m_16_days_NDVI`<br>Upload polygon: `docs_common/sierra_nevada_boundary_AppEEARS.geojson`<br>Date range: 2000-01-01 to 2026-05-31<br>Output: NetCDF-4<br>Projection: MODIS Sinusoidal | 2000-02-18 to 2026-05-25 · 605 16-day time steps · clipped to `sierra_nevada_boundary_AppEEARS.geojson` (AppEEARS upload polygon) · 500m |
| NDWI | `data_common/raw/ndwi/MOD09A1.061_500m_aid0001.nc` (Google Drive, ~1,011 MB)<br>See [NDWI_source_location.md](../data_common/raw/ndwi/NDWI_source_location.md) | [NASA AppEEARS](https://appeears.earthdatacloud.nasa.gov/) (NASA Earthdata account required)<br>Area Sample request<br>Product: Terra MODIS Surface Reflectance Bands 1-7, MOD09A1.061, 500m, 8-day (2000-02-24 to Present)<br>Layers: `sur_refl_b02` (841–876nm) + `sur_refl_b05` (1230–1250nm)<br>Upload polygon: `docs_common/sierra_nevada_boundary.geojson`<br>Date range: 2000-01-01 to 2026-05-31<br>Output: NetCDF-4<br>Projection: MODIS Sinusoidal<br>NDWI = (B02 − B05) / (B02 + B05) computed locally after download, per [Gao 1996](https://doi.org/10.1016/S0034-4257(96)00067-3) — Band 5 (1240nm) matches Gao's original SWIR reference wavelength; Terra MODIS Band 6 has malfunctioned since 2006 | 2000-02-18 to 2026-05-25 · 8-day composites · clipped to `sierra_nevada_boundary.geojson` · 500m |
| SWE *(unused)* | `data_common/raw/swe/TerraClimate_swe_2000-2025.nc` | `swe_raw <- getTerraClim(`<br>`  AOI       = sn_aoi,`<br>`  varname   = "swe",`<br>`  startDate = "2000-01-01",`<br>`  endDate   = "2025-12-31"`<br>`)` | 2000-01-01 to 2025-12-01 · 312 months · clipped to Sierra Nevada ecoregion boundary · ~4 km |
| Forest Damage | `data_common/raw/forest_damage/CONUS_Region5_AllYears.gdb.zip` (Google Drive, ~166 MB)<br>See [FOREST_DAMAGE_source_location.md](../data_common/raw/forest_damage/FOREST_DAMAGE_source_location.md) | [USDA Forest Service Detection Surveys](https://www.fs.usda.gov/science-technology/data-tools-products/fhp-mapping-reporting/detection-surveys)<br>Scroll to "R5 IDS DATA DOWNLOADS" section<br>Download "R5 CA All Years IDS Data (ZIP, 159 MB)"<br>[Direct URL](https://www.fs.usda.gov/foresthealth/docs/IDS_Data_for_Download/CONUS_Region5_AllYears.gdb.zip)<br>No login required | 1997–2025 · annual · Region 5 California · 348,829 damage polygons |
| Fire Perimeters | `data_common/raw/fire_perimeters/California_Historic_Fire_Perimeters_-4934337877712688794.gpkg` (Google Drive, ~147 MB)<br>See [FIRE_PERIMETERS_source_location.md](../data_common/raw/fire_perimeters/FIRE_PERIMETERS_source_location.md) | [CAL FIRE Hub — California Historical Fire Perimeters](https://hub-calfire-forestry.hub.arcgis.com/datasets/CALFIRE-Forestry::california-historical-fire-perimeters/about?layer=0)<br>Download → GeoPackage | California statewide, historical fire perimeter polygons |

---

### EPA Boundary

- `src_common/00_process_ecoregion.R` → `docs_common/sierra_nevada_boundary.geojson`

**References**:
- [Ecoregions of California (Griffith et al. 2016, USGS)](https://pubs.usgs.gov/publication/ofr20161021) — [Sheet 2 of 2](https://pubs.usgs.gov/of/2016/1021/ofr20161021_sheet2.pdf)
- [Ecoregions of the Conterminous United States (Omernik 1987)](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-8306.1987.tb00149.x)
- [Level III/IV refinement (Omernik & Griffith 2014)](https://link.springer.com/article/10.1007/s00267-014-0364-1)

---

### VPD — Vapor Pressure Deficit (TerraClimate)

- **Why VPD over PDSI:** VPD is a more direct measure of atmospheric water demand that drove tree mortality.

#### Raw data format
The raw file is a NetCDF-4 (`.nc`) file — a scientific format for multidimensional arrays, not inspectable like a CSV. It contains a 3D array: **312 monthly time steps × latitude grid × longitude grid**. TerraClimate is a global dataset at ~4 km resolution. Each cell in the grid holds a VPD value in kPa for that location and month. Readable in R with `terra::rast()` or `ncdf4::nc_open()`.

**References**:
- [TerraClimate — Climatology Lab](https://www.climatologylab.org/terraclimate.html)
- [Abatzoglou et al. 2018, *Scientific Data*](https://www.nature.com/articles/sdata2017191)

### PRISM — Temperature + Precipitation

- `ppt`: precipitation (rainfall + snowfall), in millimeters
- `tmean`: mean temperature (average of daily high and low), in degrees Celsius
- `tmin`: monthly minimum temperature, in degrees Celsius
- `tmax`: monthly maximum temperature, in degrees Celsius

```r
library(prism)
prism_set_dl_dir("data_common/raw/prism/prismz")
get_prism_monthlys(type = "ppt", years = 2000:2026, keepZip = FALSE)
get_prism_monthlys(type = "tmean", years = 2000:2026, keepZip = FALSE)
get_prism_monthlys(type = "tmin", years = 2000:2026, keepZip = FALSE)
get_prism_monthlys(type = "tmax", years = 2000:2026, keepZip = FALSE)
```

#### Raw data format
PRISM downloads as a collection of `.bil` (Band Interleaved by Line) raster files — one file per variable per month, saved under `data_common/raw/prism/prismz/`. Unlike a single NetCDF file, each month is its own file. Each `.bil` file is a 2D raster (latitude × longitude) covering the contiguous United States at ~4 km resolution, holding one value per pixel for that variable and month. PRISM is US-only (unlike TerraClimate which is global), but higher quality for the US because it incorporates weather station data and terrain adjustments.

**References**:
- [PRISM Climate Group — Oregon State University](https://prism.oregonstate.edu/)
- [Parameter-elevation Regressions on Independent Slopes Model (PRISM) — Adaptation Clearinghouse](https://www.adaptationclearinghouse.org/resources/parameter-elevation-regressions-on-independent-slopes-model-prism.html)

### NDVI — Normalized Difference Vegetation Index

#### Raw data format
The raw file is a NetCDF-4 (`.nc`) file with the same 3D structure as VPD: **time steps × latitude × longitude**. Unlike VPD (monthly) and PRISM (monthly), NDVI is a 16-day composite — one layer every 16 days, giving 605 time steps from 2000-02-18 to 2026-05-25. Spatial resolution is ~500m per pixel, much finer than VPD (~4km) or PRISM (~4km). AppEEARS clips the file to the submitted polygon.

### SWE — Snow Water Equivalent (TerraClimate) *(not currently used)*
- **Variable:** `swe` (Snow Water Equivalent, mm)

### USFS IDS — Forest Damage / Tree Mortality
- Raw layer (`DAMAGE_AREAS_FLAT`) covers all forest damage causal agents (insects, disease, abiotic), not just bark beetle — downstream processing filters to bark-beetle-specific damage causal agent codes.
- **Why acquired:** Bark beetle infestation is a candidate mechanism behind loss of resilience — used as ground-truth validation against modeled collapse signal.

### Fire Perimeters — CAL FIRE FRAP Historical Fire Perimeters
- **Why acquired:** distinguishes fire-caused NDVI decline from drought/climate-driven decline.

---

## Tech Stack
- **Language:** R
- **Key packages:** see `src_common/packages.R`
- **Notebooks:** R Markdown (`.Rmd`) in RStudio
