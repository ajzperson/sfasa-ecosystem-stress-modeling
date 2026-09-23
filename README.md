# SFASA 2026 Project
**Using Time-Series Statistics to Predict Ecosystem Collapse**

Case Study: The 2012–2016 Sierra Nevada Canopy Collapse

Tests whether Critical Slowing Down (CSD) statistics — trends in autocorrelation, variance, and skewness that theory predicts emerge as a system approaches a critical transition — computed from MODIS NDVI/NDWI time series (2000–2014) can serve as an early-warning signal for the 2012–2016 Sierra Nevada drought collapse. Pixels are clustered into stress tiers by their CSD trend, then checked against real 2015–2016 bark-beetle mortality (USFS ADS) and fire severity (MTBS): across every clustering variant (NDVI-only, NDWI-only, and NDVI+NDWI-combined), more-stressed tiers show a monotonically higher real-world damage rate.

Full narrative, method, and results: `pixel_approach/docs/sfasa2026_sierra_nevada_slides.pptx`

## Repo Layout

- `pixel_approach/` — primary analysis, at the 500m MODIS pixel level (251,971 pixels across the Sierra Nevada ecoregion). Covered by this README; full step-by-step pipeline design in `pixel_approach/docs/DESIGN.md`.
- `county_approach/` — a separate, county/cluster-level analysis. See `county_approach/docs/QUICKSTART.md`.
- `data_common/` — shared raw downloaded data, used by both approaches.
- `docs_common/` — shared docs; `docs_common/DATAPIPELINE.md` lists every raw data source and where to download it.

## Setup

1. `src_common/packages.R` — installs required R packages.
2. Download the raw data the pixel approach needs, per `docs_common/DATAPIPELINE.md`: EPA ecoregion boundary, NDVI, NDWI, PRISM, USFS IDS (forest damage / beetle mortality), MTBS.
3. `src_common/00_process_ecoregion.R` → `docs_common/sierra_nevada_boundary.geojson`
4. `src_common/01a_download_vpd.R` → `data_common/raw/vpd/TerraClimate_vpd_2000-2025.nc`

## Pixel Approach — Run Order

Full mechanism, parameters, and complete output list for every script below: `pixel_approach/docs/DESIGN.md`.

**1. Extract pixel-level inputs** (each pixel = one 500m MODIS grid cell; the shared grid is defined by the first script below)
- `01_extract_pixel_ndvi.R` — NDVI monthly matrix + `pixel_metadata.csv` (the shared pixel grid)
- `01_extract_pixel_ndwi.R` — NDWI monthly matrix, verified against the same grid
- `02a_pixel_vpd.R`, `02b_pixel_ppt.R` — climate-correction covariates
- `02c_pixel_fire.R` — MTBS burn exposure (fire confound + later validation)
- `02e_pixel_beetle.R` — USFS ADS beetle mortality (later validation)

**2. Compute CSD trend statistics**
- `03_pixel_stl_csd.R` — STL decomposition + climate correction + rolling CSD stats + Kendall's tau, per pixel
- `05_heatmaps.R` — spatial maps of tau_r1 / tau_var / tau_skew

**3. Cluster pixels into stress tiers**
- `06_clustering.R` — K-means per index (NDVI, NDWI separately)
- `06b_clustering_combined.R` — K-means on combined NDVI+NDWI features

**4. Validate clusters against real 2015–2016 damage**
- `08_beetle_validation.R`, `08b_beetle_validation_combined.R` — beetle mortality by cluster
- `09_fire_validation.R`, `09b_fire_validation_combined.R` — fire exposure by cluster
- `11_dashboard.R`, `11b_dashboard_combined.R` — one-page cluster map + validation summary
- `12_csd_rolling_timeline.R` — CSD trajectory over time, by cluster

**5. Compare clustering methods**
- `13_cluster_agreement_validation.R` — NDVI vs. NDWI cluster agreement
- `14_combined_vs_agreement.R` — combined clustering vs. NDVI/NDWI agreement
- `13b_agreement_venn.R`, `15_overlap_venn.R` — Venn-diagram charts; run after `14`, both read its output

**Supplementary**: `00_study_area_chart.R` — study-area map (run any time after step 1)

## Tests

`pixel_approach/tests/` — testthat unit tests for the pixel pipeline. Run all: `testthat::test_dir("pixel_approach/tests")`.
