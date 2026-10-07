# SFASA 2026 Project
**Using Time-Series Statistics to Predict Ecosystem Collapse**

Case Study: The 2012–2016 Sierra Nevada Canopy Collapse

California's 2012–2016 drought, the most severe multi-year drought in over 1,200 years, killed more than 100 million trees. Prolonged moisture deficits weakened tree defenses and triggered catastrophic bark beetle infestations, while accumulated dead biomass increased regional wildfire vulnerability. This study evaluates whether such stress can be detected from satellite and climate data before collapse becomes visible.

Critical slowing down (CSD) theory holds that as a stressed system approaches a critical transition, its restoring forces weaken and recovery from perturbation slows. This loss of resilience produces characteristic changes in time-series behavior: rising lag-1 autocorrelation, rising variance, and shifting skewness. Autocorrelation and skewness were used here; variance was excluded because its direction near a transition is not theoretically guaranteed. We hypothesized that trends in CSD statistics, calculated from satellite vegetation indices across a 2000–2014 analysis window, would predict higher rates of bark beetle mortality and fire damage during 2015–2016.

We analyzed monthly satellite NDVI (canopy greenness) and NDWI (canopy water content) from February 2000 to December 2014 for all 251,971 500-m pixels in the 53,102.7 km² Sierra Nevada ecoregion. Each pixel series was seasonally decomposed using STL (Seasonal-Trend decomposition using Loess) and adjusted for climate through per-pixel regression on vapor pressure deficit and precipitation, fit only on pre-drought growing seasons (2000–2011). Kendall's tau summarized the trend in each pixel's 54-month rolling autocorrelation and skewness over growing-season months. After excluding 30,563 pixels burned during the analysis window and residualizing against baseline vegetation density, K-means clustering assigned pixels to stress tiers, separately for NDVI and NDWI (K=3) and jointly on all four features (K=4). No data from 2015–2016 entered any step.

Damage rates increased monotonically across stress tiers in all three clusterings. Under the combined clustering, most-stressed pixels carried 1.79 times the fire risk and 1.35 times the beetle mortality risk of least-stressed pixels. The NDVI and NDWI clusterings produced comparable gradients on their own (fire RR 1.49 and 1.65; beetle mortality RR 1.29 and 1.35, most- versus least-stressed). Among the 219,637 pixels clustered by both indices, pixels flagged most-stressed by both NDVI and NDWI carried 1.37 times the fire risk and 1.14 times the beetle mortality risk of the weaker single-index flag (n = 31,146 vs. n = 59,447), exceeding either index's flag on its own.

These findings support CSD statistics derived from satellite time series as early-warning indicators of forest collapse, identifiable in advance of the damage. The analysis establishes a population-level association between stress tier and damage rate rather than pixel-level prediction; regression directly on CSD features, and deep learning applied to raw index time series, are natural next steps.

---

> **Working directory:** run every script from the repo root — all paths in the scripts are relative to it. In RStudio, open `sfasa-ecosystem-stress-modeling.Rproj` to set this automatically.

<div style="font-size: 1.5em; font-weight: 600;"><a href="docs_common/DATAPIPELINE.md">Acquiring Data</a></div>

<div style="font-size: 1.5em; font-weight: 600;"><a href="pixel_approach/docs/COMPUTATIONPIPELINE.md">Computing & Validating Early-Warning Signals Derived from Pixel Clustering</a></div>

<div style="font-size: 1.5em; font-weight: 600;"><a href="pixel_approach/docs/sfasa2026_sierra_nevada_slides.pptx">Full Narrative, Method & Results</a></div>

## Repo Layout

- `src_common/` — one-time setup/download scripts shared across the project: install R packages, process the ecoregion boundary, download VPD/SWE/PRISM climate data
- `data_common/` — raw input data
- `docs_common/` — shared docs: data-acquisition guide, background/method notes, the study-area boundary file
- `pixel_approach/` — the primary analysis, at the 500m MODIS pixel level (251,971 pixels across the Sierra Nevada ecoregion)
  - `src/` — analysis scripts
  - `data/` — processed data
  - `charts/` — chart outputs
  - `docs/` — documentation and slides
  - `tests/` — unit tests
  - `deprecated/` — old scripts no longer used
