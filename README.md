# SFASA 2026 Project
**Using Time-Series Statistics Predict Ecosystem Collapse**

Case Study: The 2012–2016 Sierra Nevada Canopy Collapse

## How to Run

## Prep

Run 1: `src_common/packages.R` — installs all required R packages.

Run 2: `src_common/00_process_ecoregion.R`
- Output: `docs_common/sierra_nevada_boundary.geojson`

Download each raw data source by following its `*_source_location.md` file (the raw data itself is too large for git):
- `data_common/raw/prism/PRISM_source_location.md` (or run `src_common/02_download_prism.R`)
- `data_common/raw/ndvi/NDVI_source_location.md`
- `data_common/raw/ndwi/NDWI_source_location.md`
- `data_common/raw/forest_damage/FOREST_DAMAGE_source_location.md`
- `data_common/raw/lcmap/LCMAP_source_location.md`
- `data_common/raw/fire_perimeters/FIRE_PERIMETERS_source_location.md`

Optional runs:
- `src_common/01a_download_vpd.R`
- `src_common/01b_download_swe.R`

## Pixel Approach

Run 1: `pixel_approach/src/01_extract_pixel_ndvi.R`
- Output:
  - `ndvi_matrix.rds`
  - `pixel_metadata.csv`

Run 2: `pixel_approach/src/01_extract_pixel_ndwi.R`
- Output: `ndwi_matrix.rds`

Run 3: `pixel_approach/src/02a_pixel_vpd.R`
- Output: `vpd_matrix.rds`

Run 4: `pixel_approach/src/02b_pixel_ppt.R`
- Output: `ppt_matrix.rds`

Run 5: `pixel_approach/src/02c_pixel_fire.R`
- Output:
  - `pixel_fire.csv`
  - `mtbs_baseline_severity.tif`
  - `mtbs_severity_2015.tif`
  - `mtbs_severity_2016.tif`
  - `mtbs_severity_2017.tif`
  - `mtbs_severity_2018.tif`
  - `mtbs_severity_2019.tif`

Run 6: `pixel_approach/src/02d_pixel_tree_mortality.R`
- Output: `pixel_tree_mortality.csv`

Run 7: `pixel_approach/src/02e_pixel_beetle.R`
- Output: `pixel_beetle.csv`

Run 8: `pixel_approach/src/03_pixel_stl_csd.R`
- Output:
  - `ndvi_pixel_tau.csv`
  - `ndwi_pixel_tau.csv`
  - `ndvi_rolling_r1.rds`
  - `ndvi_rolling_var.rds`
  - `ndvi_rolling_skew.rds`
  - `ndwi_rolling_r1.rds`
  - `ndwi_rolling_var.rds`
  - `ndwi_rolling_skew.rds`

Run 9: `pixel_approach/src/05_heatmaps.R`
- Output:
  - `heatmap_tau_r1.pdf`
  - `heatmap_tau_var.pdf`
  - `heatmap_tau_skew.pdf`

Run 10: `pixel_approach/src/06_clustering.R`
- Output:
  - `ndvi_pixel_clusters.csv`
  - `ndwi_pixel_clusters.csv`
  - `cluster_profile_ndvi.csv`
  - `cluster_profile_ndwi.csv`
  - `cluster_k_sweep.pdf`
  - `cluster_map.pdf`

Run 11: `pixel_approach/src/08_beetle_validation.R`
- Output:
  - `beetle_tpa_by_cluster.pdf`
  - `pixel_beetle_tpa_ndvi.csv`
  - `pixel_beetle_tpa_ndwi.csv`
  - `panel_beetle_ndvi.rds`
  - `panel_beetle_ndwi.rds`

Run 12: `pixel_approach/src/09_fire_validation.R`
- Output:
  - `fire_by_cluster.pdf`
  - `pixel_fire_ndvi.csv`
  - `pixel_fire_ndwi.csv`
  - `panel_fire_ndvi.rds`
  - `panel_fire_ndwi.rds`

Run 13: `pixel_approach/src/10_lcmap_validation.R`
- Output:
  - `lcmap_by_cluster.pdf`
  - `pixel_lcmap_ndvi.csv`
  - `pixel_lcmap_ndwi.csv`
  - `panel_lcmap_ndvi.rds`
  - `panel_lcmap_ndwi.rds`

Run 14: `pixel_approach/src/11_dashboard.R`
- Output:
  - `dashboard_ndvi.pdf`
  - `dashboard_ndwi.pdf`

Run 15: `pixel_approach/src/12_csd_rolling_timeline.R`
- Output: `cluster_csd_rolling_timeline.pdf`

Run 16: `pixel_approach/src/13_cluster_agreement_validation.R`
- Output:
  - `cluster_agreement_validation.pdf`
  - `cluster_agreement_map.pdf`