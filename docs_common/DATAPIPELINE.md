# Acquiring Data

## Setup and How to Run

> **Working directory:** run every script from the repo root. In RStudio, open `sfasa-ecosystem-stress-modeling.Rproj`.

| Script | Input | Output |
| --- | --- | --- |
| `packages.R` | — | installs required R packages — **run this first** |
| **(manual)** | — | Go to each `data_common/raw/<source>/` folder, look for a `*_source_location.md` file, follow the Google Drive link inside it, and download into that folder — skip any source marked **not used**. |
| `00_process_ecoregion.R` (optional) | `data_common/raw/ecoregion/us_eco_l3.zip` | `data_common/processed/sierra_nevada_boundary.geojson`<br>• already tracked in git; only run this if you want to regenerate it yourself |
| `01a_download_vpd.R` (optional) | `data_common/processed/sierra_nevada_boundary.geojson` | `data_common/raw/vpd/TerraClimate_vpd_2000-2025.nc`<br>• already tracked in git; only run this if you want to redownload it yourself |
| `01b_download_swe.R` (optional) | `data_common/processed/sierra_nevada_boundary.geojson` | `data_common/raw/swe/TerraClimate_swe_2000-2025.nc`<br>• already tracked in git; only run this if you want to redownload it yourself |
| `02_download_prism.R` (optional) | — | `data_common/raw/prism/prismz/`<br>• only run this if you don't want to download PRISM manually via `PRISM_source_location.md` in the **(manual)** step above |

## Deep Dive of Data Sources

1. **EPA Boundary** — EPA Level III Ecoregions shapefile (full continental US)<br>
   EPA = Environmental Protection Agency (the US federal agency)<br>
   Used to generate the Sierra Nevada ecoregion polygon (spatial framework/mask that every raw dataset gets clipped/extracted to).
2. **Kaggle** — California Wildfire Incidents (downloaded, **not used**)
3. **VPD** — Vapor Pressure Deficit (TerraClimate)
4. **PRISM** — Temperature + Precipitation
5. **NDVI** — Normalized Difference Vegetation Index (NASA MODIS)
6. **NDWI** — Normalized Difference Water Index (NASA MODIS)
7. **SWE** — Snow Water Equivalent (TerraClimate) (**not used**)
8. **USFS IDS** — Forest Damage / Tree Mortality (Insect and Disease Survey)<br>
   USFS IDS = US Forest Service, Insect and Disease Survey (also called the Aerial Detection Survey/ADS)
9. **Fire Perimeters** — California historical fire perimeter polygons (CAL FIRE FRAP) (superseded by MTBS, **not used**)
10. **MTBS** — Monitoring Trends in Burn Severity (pixel-level burn severity)

### Raw Data Inventory

| Data | File | Download | Spatial scope |
|------|------|----------|---------------|
| EPA Boundary | `data_common/raw/ecoregion/us_eco_l3.zip` | [EPA Level III and IV Ecoregions of the Continental United States](https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states)<br>Download "US Level III Ecoregions shapefile without state boundaries" | Full CONUS |
| Kaggle (**not used**) | `data_common/raw/kaggle/California_Fire_Incidents_2013-2019_raw.csv` | [Kaggle](https://www.kaggle.com/datasets/ananthu017/california-wildfire-incidents-20132020) | California statewide |
| VPD | `data_common/raw/vpd/TerraClimate_vpd_2000-2025.nc` | `vpd_raw <- getTerraClim(`<br>`  AOI       = sn_aoi,`<br>`  varname   = "vpd",`<br>`  startDate = "2000-01-01",`<br>`  endDate   = "2025-12-31"`<br>`)` | 2000-01-01 to 2025-12-01 · 312 months · cropped to the Sierra Nevada ecoregion's bounding box · ~4km |
| PRISM | `data_common/raw/prism/prismz.zip`<br>See [PRISM_source_location.md](../data_common/raw/prism/PRISM_source_location.md) | `get_prism_monthlys(type = "ppt", years = 2000:2026, keepZip = FALSE)`<br>`get_prism_monthlys(type = "tmean", years = 2000:2026, keepZip = FALSE)`<br>`get_prism_monthlys(type = "tmin", years = 2000:2026, keepZip = FALSE)`<br>`get_prism_monthlys(type = "tmax", years = 2000:2026, keepZip = FALSE)` | 2000-01-01 to 2026-05-01 · 317 months · Full CONUS · ~4km |
| NDVI | `data_common/raw/ndvi/MOD13A1.061_500m_aid0001.nc`<br>See [NDVI_source_location.md](../data_common/raw/ndvi/NDVI_source_location.md) | [NASA AppEEARS](https://appeears.earthdatacloud.nasa.gov/) (NASA Earthdata account required)<br>Area Sample request<br>Product: Terra MODIS Vegetation Indices (NDVI & EVI), MOD13A1.061, 500m, 16-day (2000-02-18 to Present)<br>Layer: `_500m_16_days_NDVI`<br>Upload polygon: `data_common/processed/sierra_nevada_boundary.geojson`<br>Date range: 2000-01-01 to 2026-05-31<br>Output: NetCDF-4<br>Projection: MODIS Sinusoidal | 2000-02-18 to 2026-05-25 · 605 16-day time steps · clipped to `sierra_nevada_boundary.geojson` · 500m |
| NDWI | `data_common/raw/ndwi/MOD09A1.061_500m_aid0001.nc`<br>See [NDWI_source_location.md](../data_common/raw/ndwi/NDWI_source_location.md) | [NASA AppEEARS](https://appeears.earthdatacloud.nasa.gov/) (NASA Earthdata account required)<br>Area Sample request<br>Product: Terra MODIS Surface Reflectance Bands 1-7, MOD09A1.061, 500m, 8-day (2000-02-18 to Present)<br>Layers: `sur_refl_b02` (841–876nm) + `sur_refl_b05` (1230–1250nm)<br>Upload polygon: `data_common/processed/sierra_nevada_boundary.geojson`<br>Date range: 2000-01-01 to 2026-05-31<br>Output: NetCDF-4<br>Projection: MODIS Sinusoidal<br>NDWI = (B02 − B05) / (B02 + B05) computed locally after download, per [Gao 1996](https://doi.org/10.1016/S0034-4257(96)00067-3) — Gao's formula specifies 0.86μm and 1.24μm exactly, which Band 2 (858nm center) and Band 5 (1240nm center) match; Band 6 (1628–1652nm) is a different SWIR region used by other, separately-named indices | 2000-02-18 to 2026-05-25 · 1,206 8-day time steps · clipped to `sierra_nevada_boundary.geojson` · 500m |
| SWE (**not used**) | `data_common/raw/swe/TerraClimate_swe_2000-2025.nc` | `swe_raw <- getTerraClim(`<br>`  AOI       = sn_aoi,`<br>`  varname   = "swe",`<br>`  startDate = "2000-01-01",`<br>`  endDate   = "2025-12-31"`<br>`)` | 2000-01-01 to 2025-12-01 · 312 months · cropped to the Sierra Nevada ecoregion's bounding box · ~4km |
| USFS IDS | `data_common/raw/forest_damage/CONUS_Region5_AllYears.gdb.zip`<br>See [FOREST_DAMAGE_source_location.md](../data_common/raw/forest_damage/FOREST_DAMAGE_source_location.md) | [USDA Forest Service Detection Surveys](https://www.fs.usda.gov/science-technology/data-tools-products/fhp-mapping-reporting/detection-surveys)<br>Scroll to "R5 IDS DATA DOWNLOADS" section<br>Download "R5 CA All Years IDS Data (ZIP)"<br>[Direct URL](https://www.fs.usda.gov/foresthealth/docs/IDS_Data_for_Download/CONUS_Region5_AllYears.gdb.zip)<br>No login required | 1997–2025 · annual · 348,829 damage polygons (99.2% Region 5 California; 0.8% spillover from neighboring Regions 4/6) |
| Fire Perimeters (**not used**) | `data_common/raw/fire_perimeters/California_Historic_Fire_Perimeters_-4934337877712688794.gpkg`<br>See [FIRE_PERIMETERS_source_location.md](../data_common/raw/fire_perimeters/FIRE_PERIMETERS_source_location.md) | [CAL FIRE Hub — California Historical Fire Perimeters](https://hub-calfire-forestry.hub.arcgis.com/datasets/CALFIRE-Forestry::california-historical-fire-perimeters/about?layer=0)<br>Download → GeoPackage | California statewide, historical fire perimeter polygons |
| MTBS | `data_common/raw/mtbs/bsp_ugjds3epryqgif0ej5sb.zip` | [MTBS Burn Severity Direct Download](https://burnseverity.cr.usgs.gov/direct-download)<br>Click "MTBS" → select California → select years 2000–2019 (single multi-year request supported)<br>Request is processed asynchronously (~1 hour); download link sent via email | California, 2000–2019 · annual thematic burn severity mosaics · GeoTIFF · 30m |

---

### EPA Boundary

- **Variable:** none — polygon geometry only (ecoregion boundary), used as a spatial mask.
- `src_common/00_process_ecoregion.R` → `data_common/processed/sierra_nevada_boundary.geojson`

- **References:**
  - [Ecoregions of California (Griffith et al. 2016, USGS)](https://pubs.usgs.gov/publication/ofr20161021) — [Sheet 2 of 2](https://pubs.usgs.gov/of/2016/1021/ofr20161021_sheet2.pdf)
  - [Ecoregions of the Conterminous United States (Omernik 1987)](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-8306.1987.tb00149.x)
  - [Level III/IV refinement (Omernik & Griffith 2014)](https://link.springer.com/article/10.1007/s00267-014-0364-1)

### Kaggle — California Wildfire Incidents (**not used**)
- **Variable:** none — not loaded by any script.

### VPD — Vapor Pressure Deficit (TerraClimate)

- **Variable:** `vpd` (Vapor Pressure Deficit, kPa)
- `src_common/01a_download_vpd.R` → `data_common/raw/vpd/TerraClimate_vpd_2000-2025.nc`

- **References:**
  - [TerraClimate — Climatology Lab](https://www.climatologylab.org/terraclimate.html)
  - [Abatzoglou et al. 2018, *Scientific Data*](https://www.nature.com/articles/sdata2017191)

### PRISM — Temperature + Precipitation

- `ppt`: precipitation (rainfall + snowfall), in millimeters
- `tmean`: mean temperature (average of daily high and low), in degrees Celsius
- `tmin`: minimum temperature (daily minimum averaged over the month), in degrees Celsius
- `tmax`: maximum temperature (daily maximum averaged over the month), in degrees Celsius
- PRISM is US-only (unlike TerraClimate which is global), but higher quality for the US because it incorporates weather station data and terrain adjustments.
- `src_common/02_download_prism.R` → `data_common/raw/prism/prismz/`

- **References:**
  - [PRISM Climate Group — Oregon State University](https://prism.oregonstate.edu/)
  - [Parameter-elevation Regressions on Independent Slopes Model (PRISM) — Adaptation Clearinghouse](https://www.adaptationclearinghouse.org/resources/parameter-elevation-regressions-on-independent-slopes-model-prism.html)

### NDVI — Normalized Difference Vegetation Index (NASA MODIS)
- **Variable:** `NDVI` (unitless index, range -1 to 1)

### NDWI — Normalized Difference Water Index (NASA MODIS)
- **Variable:** `NDWI` (unitless index, range -1 to 1; NDWI = (B02 − B05) / (B02 + B05), see formula in table above)

### SWE — Snow Water Equivalent (TerraClimate) (**not used**)
- **Variable:** `swe` (Snow Water Equivalent, mm)
- `src_common/01b_download_swe.R` → `data_common/raw/swe/TerraClimate_swe_2000-2025.nc`

### USFS IDS — Forest Damage / Tree Mortality
- **Variable:** `LEGACY_TPA` (trees per acre killed)
- Raw layer (`DAMAGE_AREAS_FLAT_AllYears_CONUS_Rgn5`) covers all forest damage causal agents (insects, disease, abiotic), not just bark beetle — downstream processing filters to bark-beetle-specific damage causal agent codes.
- Bark beetle infestation is a candidate mechanism behind loss of resilience — used as ground-truth validation against the modeled early-warning signal.

### Fire Perimeters — CAL FIRE FRAP Historical Fire Perimeters (superseded by MTBS, **not used**)
- **Variable:** none — polygon geometry only (fire boundary), no per-pixel numeric field used.

### MTBS — Monitoring Trends in Burn Severity
- **Variable:** categorical severity class code per pixel, per MTBS's own raster metadata (this pipeline treats 3/4 as "burned"):

| Code | Meaning |
|------|---------|
| 0 | Background/No Data |
| 1 | Unburned/Underburned to Low Burn Severity |
| 2 | Low burn severity |
| 3 | Moderate burn severity |
| 4 | High Burn Severity |
| 5 | Increased Greenness/Increased Vegetation Response |
| 6 | Non-Processing Area Mask |
- MTBS gives an actual severity class per 30m pixel.
- Organized as one GeoTIFF mosaic per US state per year; 500m study pixels are aggregated from the underlying 30m MTBS pixels (~17×17 per study pixel).

- **References:**
  - [MTBS home](https://www.mtbs.gov/)
  - [MTBS Product Descriptions](https://www.mtbs.gov/product-descriptions)
