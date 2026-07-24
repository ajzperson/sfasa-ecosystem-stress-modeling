# Pixel Pipeline — Design Flow Chart

**Pixel** = one 500m MODIS grid cell (NDVI's and NDWI's native resolution) within the Sierra Nevada boundary — 251,971 of them total. Identified by `pixel_id` (row position in every matrix/CSV below) and a fixed lon/lat centroid, both set once in `pixel_metadata.csv` and reused unchanged by every downstream script.

Every dataset below is resolved onto this same grid before use — everything in this pipeline is at the pixel level: VPD and PRISM (~4km native resolution) are bilinear-interpolated at each pixel centroid; MTBS and LCMAP (30m GeoTIFFs) are aggregated up to this grid via a 250m-radius buffer around each centroid; Forest Damage (vector polygons, no fixed resolution) is joined by intersecting each pixel centroid with the polygon layer. `pixel_id` is what lets every one of these datasets join back together row-for-row.

## Generate monthly NDVI for the Sierra Nevada region (Feb 2000–Dec 2016) from "NDVI" in ["Raw Data Inventory"](../../docs_common/DATAPIPELINE.md#raw-data-inventory)
- QC-mask (keep quality bits 00/01/10, drop 11), crop/mask to the SN polygon, aggregate 16-day composites to monthly means.
- `pixel_approach/src/01_extract_pixel_ndvi.R` →
  - `pixel_approach/data/processed/ndvi_matrix.rds`
  - `pixel_approach/data/processed/pixel_metadata.csv`

## Generate monthly NDWI for the Sierra Nevada region (Feb 2000–Dec 2016) from "NDWI" in ["Raw Data Inventory"](../../docs_common/DATAPIPELINE.md#raw-data-inventory)
- Compute NDWI = (b02 − b05) / (b02 + b05), QC-mask (bit 00 only), crop/mask to the SN polygon, aggregate 8-day composites to monthly means, verify against NDVI's shared pixel grid.
- `pixel_approach/src/01_extract_pixel_ndwi.R` →
  - `pixel_approach/data/processed/ndwi_matrix.rds`

## Generate monthly VPD (Feb 2000–Dec 2016) at each pixel centroid from "VPD" in ["Raw Data Inventory"](../../docs_common/DATAPIPELINE.md#raw-data-inventory)
- Crop TerraClimate (~4km) to the pixel centroids' bounding box (padded 0.5°), bilinear-interpolate a value at each centroid.
- `pixel_approach/src/02a_pixel_vpd.R` →
  - `pixel_approach/data/processed/vpd_matrix.rds`

## Generate monthly precipitation (Feb 2000–Dec 2016) at each pixel centroid from "PRISM" in ["Raw Data Inventory"](../../docs_common/DATAPIPELINE.md#raw-data-inventory)
- Project pixel centroids to PRISM's CRS (NAD83), crop PRISM (~4km) to their bounding box (padded 0.5°), bilinear-interpolate a value at each centroid.
- `pixel_approach/src/02b_pixel_ppt.R` →
  - `pixel_approach/data/processed/ppt_matrix.rds`

## Generate per-pixel burn exposure from "MTBS" in ["Raw Data Inventory"](../../docs_common/DATAPIPELINE.md#raw-data-inventory)
- Reclassify severity to binary (moderate/high = burned), aggregate 30m cells into each 500m pixel via a 250m-radius buffer max; baseline (2000-2014) OR'd into one flag, validation years (2015-2019) kept per-year.
- `pixel_approach/src/02c_pixel_fire.R` →
  - `pixel_approach/data/processed/pixel_fire.csv`
  - `pixel_approach/data/processed/mtbs_baseline_severity.tif` (cache, avoids reprocessing MTBS baseline on rerun)
  - `pixel_approach/data/processed/mtbs_severity_2015.tif` (cache, one per validation year)
  - `pixel_approach/data/processed/mtbs_severity_2016.tif`
  - `pixel_approach/data/processed/mtbs_severity_2017.tif`
  - `pixel_approach/data/processed/mtbs_severity_2018.tif`
  - `pixel_approach/data/processed/mtbs_severity_2019.tif`

## Generate per-pixel beetle mortality from "Forest Damage" in ["Raw Data Inventory"](../../docs_common/DATAPIPELINE.md#raw-data-inventory)
- Filter USFS ADS damage polygons to Region 5, survey years 2015-2016, `DAMAGE_TYPE == "Mortality"`, and bark-beetle DCA codes only (17 codes, see table below); clip to the Sierra Nevada boundary; spatial-join each pixel centroid to overlapping polygons (`st_intersects`), taking max `tpa` (trees per acre killed, from `LEGACY_TPA`) per pixel per year.

| DCA_CODE | DCA_COMMON_NAME |
|---|---|
| 11000 | bark beetles |
| 11002 | western pine beetle |
| 11004 | Jeffrey pine beetle |
| 11006 | mountain pine beetle |
| 11007 | Douglas-fir beetle |
| 11019 | pinyon ips |
| 11027 | California fivespined ips |
| 11029 | pine engraver |
| 11030 | ips engraver beetles |
| 11035 | cedar & cypress bark beetles |
| 11036 | western cedar bark beetle |
| 11039 | twig beetles |
| 11041 | fir root bark beetle |
| 11048 | true fir bark beetles |
| 11049 | Douglas-fir engraver |
| 11050 | fir engraver |
| 11900 | unknown bark beetle |

- `pixel_approach/src/02e_pixel_beetle.R` →
  - `pixel_approach/data/processed/pixel_beetle.csv`

## Compute per-pixel CSD (critical slowing down) statistics from the NDVI/NDWI matrices + VPD + PPT

Runs STL decomposition + climate correction + rolling CSD statistics per pixel.

**Four periods, four different roles:**

| Period | Role | Used for |
|---|---|---|
| 2000–2011 | Healthy baseline | Fits the climate model (`lm`) |
| 2000–2014 | Full STL/CSD analysis window | STL fit, growing-season filter, climate-model prediction, rolling stats, Kendall's tau |
| 2012–2014 | Drought onset (inside 2000–2014) | Included in STL fit and predicted against, but excluded from fitting the climate model |
| 2015–2016 | Collapse/validation period | Excluded entirely from `03` (no STL, no climate fit) |

**Mechanism**
- index = NDVI or NDWI (shown below as `ndvi` for brevity; same VPD/PPT data used for both)
- step-by-step, for each pixel:
  1. 2000–2014 (all months): `ndvi` → run STL → `seasonal`
  2. 2000–2014 (all months): `ndvi − seasonal` → `ndvi_deseas`
  3. 2000–2014 (Mar–Nov): `ndvi_deseas` → filter to growing season → `ndvi_gs`
  4. 2000–2011 (Mar–Nov): `ndvi_gs, vpd, ppt` → fit `lm` → `b0, b1, b2`
  5. 2000–2014 (Mar–Nov): `vpd, ppt` → `b0 + b1·vpd + b2·ppt` → `predicted`
  6. 2000–2014 (Mar–Nov): `ndvi_gs − predicted` → `clim_resid`
  7. 2000–2014 (Mar–Nov): `clim_resid` → 54-month rolling window → `r1, var, skew`
     - `r1` skips the Nov→Mar cross-year gap so it isn't miscounted as a normal 1-month lag (Dec–Feb are already excluded by the growing-season filter, so Nov and Mar sit array-adjacent but are 4 calendar months apart).
     - `var`/`skew` don't need this — they're order-independent (computed from the unordered set of values in the window, not adjacent pairs).
  8. 2000–2014 (Mar–Nov): `r1, var, skew` over time → Kendall's tau → `tau_r1, tau_var, tau_skew`

**Parameters:**
- `WIN = 54` months rolling window
- `NA_THRESHOLD = 0.30` — pixels with ≥ 30% NA months are skipped
- `EDGE_NA_LIMIT = 6` — pixels with a leading/trailing contiguous NA run longer than
  6 months are skipped (data guardrail, see below)
- `N_WORKERS = 5` — parallel worker processes
- `FUTURE_MAX_SIZE_MB = 2048` — raises `future`'s default 500 MiB safety cap so the
  parallel workers are allowed to receive the data they need.

**Output:** `pixel_approach/src/03_pixel_stl_csd.R` →
- `ndvi_pixel_tau.csv`
- `ndwi_pixel_tau.csv`
- `ndvi_rolling_r1.rds`
- `ndvi_rolling_var.rds`
- `ndvi_rolling_skew.rds`
- `ndwi_rolling_r1.rds`
- `ndwi_rolling_var.rds`
- `ndwi_rolling_skew.rds`

## Visualize spatial CSD trends
- Map tau_r1, tau_var, tau_skew across the Sierra Nevada, NDVI vs. NDWI side by side, shared color scale.
- `pixel_approach/src/05_heatmaps.R` →
  - `heatmap_tau_r1.pdf`
  - `heatmap_tau_var.pdf`
  - `heatmap_tau_skew.pdf`

## Cluster pixels into stress tiers

Turns each pixel's CSD statistics into a stress tier (K-means), run independently per index.

**Mechanism, step by step** (index = NDVI or NDWI, run independently for each):
1. Load `tau_r1, tau_skew` for every pixel.
2. Exclude pixels burned (MTBS moderate+) at any point 2000-2014, using the `burned_baseline_00_14` flag from `pixel_fire.csv` (`02c`'s output).
   - A fire causes a sudden, discrete NDVI/NDWI crash unrelated to gradual drought stress, which would contaminate the tau_r1/tau_skew trend statistics computed over the same 2000-2014 window.
3. Compute `base_level` per pixel: mean NDVI/NDWI value over 2000–2014.
4. Residualize: `resid_tau_r1 = resid(lm(tau_r1 ~ base_level))`, `resid_tau_skew = resid(lm(tau_skew ~ base_level))`.
   - tau_r1 and tau_skew both correlated with a pixel's own baseline vegetation density — denser forest mechanically has more fuel to burn and more tree cover eligible to register a loss, regardless of any actual CSD dynamics. Pearson correlation with `base_level`, across all pixels with a valid tau value:

     | Index | cor(tau_r1, base_level) | cor(tau_skew, base_level) |
     |---|---|---|
     | NDVI | 0.30 | -0.15 |
     | NDWI | 0.31 | -0.19 |

     Residualizing against `base_level` removes this confound (post-residualization correlation with `base_level` ~0, verified).
5. Flip sign: `neg_resid_tau_skew = -resid_tau_skew`, so higher = more stressed, consistent with `resid_tau_r1`'s convention.
6. Z-score both features (`resid_tau_r1`, `neg_resid_tau_skew`).
7. K-means sweep K=2-8 (elbow + silhouette, for the chart), then fit K-means at K=3 (the chosen value).
8. Relabel clusters 1..3 by stress rank: `composite = mean(resid_tau_r1) - mean(resid_tau_skew)` per raw cluster, ascending — 1 = least stressed, 3 = most stressed.

**Parameters:**
- `K = 3`
- `SEED = 42`
- `nstart = 25`, `iter.max = 100`
- `INDICES = c("ndvi", "ndwi")` — which index/indices to cluster; run independently, once per index

**Output:** `pixel_approach/src/06_clustering.R` →
- `ndvi_pixel_clusters.csv`
- `ndwi_pixel_clusters.csv`
- `cluster_profile_ndvi.csv`
- `cluster_profile_ndwi.csv`
- `cluster_k_sweep.pdf`
- `cluster_map.pdf`

## Validate clusters against beetle mortality (2015-2016)

Checks whether the stress-tier clusters (from `06`) actually line up with real bark-beetle kill, using `02e`'s per-pixel mortality data as ground truth.

**Mechanism, step by step** (index = NDVI or NDWI, run independently for each):
1. A pixel counts as having beetle mortality if it had any bark-beetle kill in 2015 or 2016 (the years the ADS data reliably covers).
2. For each cluster, compute the % of its pixels with mortality — the actual validation question: do more-stressed clusters really show more real-world mortality?
3. Plot that percentage by cluster, with the average kill severity (TPA) shown for context on each bar.

**Output:** `pixel_approach/src/08_beetle_validation.R` →
- `beetle_tpa_by_cluster.pdf`
- `pixel_beetle_tpa_ndvi.csv`
- `pixel_beetle_tpa_ndwi.csv`
- `panel_beetle_ndvi.rds` (saved plot panel, read by `11_dashboard.R`)
- `panel_beetle_ndwi.rds` (saved plot panel, read by `11_dashboard.R`)

## Validate clusters against fire exposure (2015-2016)

Checks whether the stress-tier clusters (from `06`) line up with real MTBS burn severity, a second independent ground-truth signal alongside beetle mortality.

**Mechanism, step by step** (index = NDVI or NDWI, run independently for each):
1. A pixel counts as burned if MTBS recorded moderate+ severity in 2015 or 2016.
2. For each cluster, compute the % of its pixels burned — the validation question: do more-stressed clusters really show more real-world fire?
3. Plot that percentage by cluster.

**Output:** `pixel_approach/src/09_fire_validation.R` →
- `fire_by_cluster.pdf`
- `pixel_fire_ndvi.csv`
- `pixel_fire_ndwi.csv`
- `panel_fire_ndvi.rds` (saved plot panel, read by `11_dashboard.R`)
- `panel_fire_ndwi.rds` (saved plot panel, read by `11_dashboard.R`)

## Assemble a one-page validation summary per index
- Combine the cluster map with the fire/beetle bar panels into a single page.
- `pixel_approach/src/11_dashboard.R` →
  - `dashboard_ndvi.pdf`
  - `dashboard_ndwi.pdf`

## Cluster pixels into stress tiers, combined NDVI+NDWI

Exploratory alternative to the per-index clustering above: instead of clustering NDVI and NDWI independently (2 features each), combines both indices' features into one 4-feature K-means clustering. Does not modify `06_clustering.R` or its outputs — a separate approach, run and compared alongside it.

**Mechanism, step by step:**
1. Load `tau_r1, tau_skew` for both NDVI and NDWI; inner-join on `pixel_id` (only pixels with a valid tau for both indices).
2. Exclude pixels burned (MTBS moderate+) at any point 2000-2014, same as `06`.
3. Compute `base_level_ndvi`, `base_level_ndwi` per pixel: mean value over 2000–2014, each index against its own matrix.
4. Residualize each index's `tau_r1`/`tau_skew` against its own `base_level` (not the other index's) — same reasoning as `06`, applied twice.
5. Flip sign on both skew residuals so higher = more stressed.
6. Z-score all 4 features, K-means sweep K=2-8, fit K-means at the chosen K.
7. Relabel clusters by stress rank using PC1 of the same 4 z-scored features (not a raw unweighted sum) — NDVI's and NDWI's features are only weakly correlated (r1: r=0.20; skew: r=-0.15), so an unweighted sum isn't a well-supported way to combine them; PC1 lets the shared-variance structure set each feature's weight instead.

**Parameters:**
- `K = 4`, `SEED = 42`, `nstart = 25`, `iter.max = 100`

**Output:** `pixel_approach/src/06b_clustering_combined.R` →
- `combined_pixel_clusters.csv`
- `combined_cluster_profile.csv`
- `combined_cluster_k_sweep.pdf`
- `combined_cluster_map.pdf`

## Validate combined clusters against beetle mortality (2015-2016)

Same role as `08_beetle_validation.R`, but for the combined clustering — one cluster assignment per pixel, so one chart, not two panels.

**Output:** `pixel_approach/src/08b_beetle_validation_combined.R` →
- `combined_beetle_tpa_by_cluster.pdf`
- `combined_pixel_beetle_tpa.csv`
- `panel_beetle_combined.rds` (saved plot panel, read by `11b_dashboard_combined.R`)

## Validate combined clusters against fire exposure (2015-2016)

Same role as `09_fire_validation.R`, but for the combined clustering.

**Output:** `pixel_approach/src/09b_fire_validation_combined.R` →
- `combined_fire_by_cluster.pdf`
- `combined_pixel_fire.csv`
- `panel_fire_combined.rds` (saved plot panel, read by `11b_dashboard_combined.R`)

## Assemble a one-page validation summary, combined
- Same role as `11_dashboard.R`, but for the combined clustering — one cluster assignment per pixel, so one output file, not two.
- `pixel_approach/src/11b_dashboard_combined.R` →
  - `dashboard_combined.pdf`

## Trace CSD statistics over time, by cluster
- Reconstruct each pixel's climate-corrected residual from saved STL/climate coefficients, average the rolling level/r1/skew trajectories within each cluster (2005–2014) to see whether the trend predates the 2012 drought onset.
- `pixel_approach/src/12_csd_rolling_timeline.R` →
  - `cluster_csd_rolling_timeline.pdf`

## Check agreement between NDVI's and NDWI's cluster assignments
- Split pixels by whether NDVI, NDWI, or both flagged them "full CSD," compare real-world outcome rates and spatial pattern between the groups.
- `pixel_approach/src/13_cluster_agreement_validation.R` →
  - `cluster_agreement_validation.pdf`
  - `cluster_agreement_map.pdf`

---

## Script Outputs

| Script | Output |
|---|---|
| `pixel_approach/src/01_extract_pixel_ndvi.R` (calls `extract_pixel_matrix.R`) | `pixel_approach/data/processed/ndvi_matrix.rds`<br>`pixel_approach/data/processed/pixel_metadata.csv` |
| `pixel_approach/src/01_extract_pixel_ndwi.R` (calls `extract_pixel_matrix.R`) | `pixel_approach/data/processed/ndwi_matrix.rds` |
| `pixel_approach/src/02a_pixel_vpd.R` | `pixel_approach/data/processed/vpd_matrix.rds` |
| `pixel_approach/src/02b_pixel_ppt.R` | `pixel_approach/data/processed/ppt_matrix.rds` |
| `pixel_approach/src/02c_pixel_fire.R` | `pixel_approach/data/processed/pixel_fire.csv`<br>`pixel_approach/data/processed/mtbs_baseline_severity.tif`<br>`pixel_approach/data/processed/mtbs_severity_2015.tif`<br>`pixel_approach/data/processed/mtbs_severity_2016.tif`<br>`pixel_approach/data/processed/mtbs_severity_2017.tif`<br>`pixel_approach/data/processed/mtbs_severity_2018.tif`<br>`pixel_approach/data/processed/mtbs_severity_2019.tif` |
| `pixel_approach/src/02e_pixel_beetle.R` | `pixel_approach/data/processed/pixel_beetle.csv` |
| `pixel_approach/src/03_pixel_stl_csd.R` | `pixel_approach/data/processed/ndvi_pixel_tau.csv`<br>`pixel_approach/data/processed/ndwi_pixel_tau.csv`<br>`pixel_approach/data/processed/ndvi_rolling_r1.rds`<br>`pixel_approach/data/processed/ndvi_rolling_var.rds`<br>`pixel_approach/data/processed/ndvi_rolling_skew.rds`<br>`pixel_approach/data/processed/ndwi_rolling_r1.rds`<br>`pixel_approach/data/processed/ndwi_rolling_var.rds`<br>`pixel_approach/data/processed/ndwi_rolling_skew.rds` |
| `pixel_approach/src/05_heatmaps.R` | `pixel_approach/charts/heatmap_tau_r1.pdf`<br>`pixel_approach/charts/heatmap_tau_var.pdf`<br>`pixel_approach/charts/heatmap_tau_skew.pdf` |
| `pixel_approach/src/06_clustering.R` | `pixel_approach/data/processed/ndvi_pixel_clusters.csv`<br>`pixel_approach/data/processed/ndwi_pixel_clusters.csv`<br>`pixel_approach/data/processed/cluster_profile_ndvi.csv`<br>`pixel_approach/data/processed/cluster_profile_ndwi.csv`<br>`pixel_approach/charts/cluster_k_sweep.pdf`<br>`pixel_approach/charts/cluster_map.pdf` |
| `pixel_approach/src/08_beetle_validation.R` | `pixel_approach/charts/beetle_tpa_by_cluster.pdf`<br>`pixel_approach/data/processed/pixel_beetle_tpa_ndvi.csv`<br>`pixel_approach/data/processed/pixel_beetle_tpa_ndwi.csv`<br>`pixel_approach/data/processed/panel_beetle_ndvi.rds`<br>`pixel_approach/data/processed/panel_beetle_ndwi.rds` |
| `pixel_approach/src/09_fire_validation.R` | `pixel_approach/charts/fire_by_cluster.pdf`<br>`pixel_approach/data/processed/pixel_fire_ndvi.csv`<br>`pixel_approach/data/processed/pixel_fire_ndwi.csv`<br>`pixel_approach/data/processed/panel_fire_ndvi.rds`<br>`pixel_approach/data/processed/panel_fire_ndwi.rds` |
| `pixel_approach/src/10_lcmap_validation.R` | `pixel_approach/charts/lcmap_by_cluster.pdf`<br>`pixel_approach/data/processed/pixel_lcmap_ndvi.csv`<br>`pixel_approach/data/processed/pixel_lcmap_ndwi.csv`<br>`pixel_approach/data/processed/panel_lcmap_ndvi.rds`<br>`pixel_approach/data/processed/panel_lcmap_ndwi.rds` |
| `pixel_approach/src/11_dashboard.R` | `pixel_approach/charts/dashboard_ndvi.pdf`<br>`pixel_approach/charts/dashboard_ndwi.pdf` |
| `pixel_approach/src/12_csd_rolling_timeline.R` | `pixel_approach/charts/cluster_csd_rolling_timeline.pdf` |
| `pixel_approach/src/13_cluster_agreement_validation.R` | `pixel_approach/charts/cluster_agreement_validation.pdf`<br>`pixel_approach/charts/cluster_agreement_map.pdf` |

---

## Why SWE is used in the combined clustering but not the separate NDVI/NDWI clustering

`03b_pixel_stl_csd_swe.R` adds SWE (snow water equivalent) as a third climate-correction predictor alongside VPD and PPT, on top of the same 2000–2011 baseline regression `03_pixel_stl_csd.R` already uses. SWE meaningfully improves the climate-correction fit for both NDVI and NDWI, by a comparable amount for each (adjusted R² rose from 0.35 to 0.48 for NDVI and from 0.39 to 0.48 for NDWI), and the SWE coefficient is statistically significant for the large majority of pixels (83% NDVI, 75% NDWI) even with VPD and PPT already in the model — not an artifact of adding a third predictor, and not lopsided toward one index.

That improvement carries through cleanly to `06b_clustering_combined.R` (the 4-feature combined clustering): beetle mortality and MTBS fire exposure both still increase monotonically with stress rank, checked at both K=2 and K=4. It does not carry through cleanly to `06_clustering.R` (the 2-feature separate NDVI/NDWI clustering): with SWE-corrected input, NDWI's beetle validation and NDVI's fire validation both stopped increasing monotonically, with the middle cluster (of 3) dipping below the least-stressed cluster.

The reason traces back to `tau_skew`. Across this whole pipeline, `tau_r1` (rising autocorrelation) is the feature that actually tracks real ground-truth damage; `tau_skew` carries little independent signal on its own — this was established from the combined clustering's PCA loadings, where both skew features loaded near zero while both r1 features loaded strongly and roughly equally. In the 2-feature separate clustering, `resid_tau_r1` and `neg_resid_tau_skew` get equal weight by construction (both z-scored before `kmeans()`), so when SWE reshapes the residual dynamics enough to shift `tau_skew`, an equal-weighted clustering can end up separating two groups mainly along the uninformative skew axis rather than the informative r1 axis — producing clusters that aren't really ordered by real stress, even though the labeling logic still ranks them as if they were. The 4-feature combined clustering doesn't have this failure mode the same way: NDVI's and NDWI's r1 features reinforce each other (a real, if weak, cross-index correlation), so the informative axis has more relative pull in the combined feature space than either index's own noisier skew feature does alone.

So SWE is used in `06b` and not in `06`: it's a real, verified improvement to the climate correction itself, and the combined clustering is robust enough to benefit from it cleanly, but the separate clustering is not.
