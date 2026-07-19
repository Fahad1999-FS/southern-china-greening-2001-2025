# Vegetation Greening and Climate Covariability in Subtropical Southern China (2001–2025)

Data, code, and figures supporting the manuscript:

> **Spatiotemporal Dynamics of Vegetation Greening and Its Seasonally Structured Climate Covariability in Subtropical Southern China (2001–2025)**
> Mehmood K., et al. *(under review, Earth Science Informatics)*

---

## Overview

This repository documents a seasonally resolved reassessment of vegetation greening across subtropical southern China using a **single, internally consistent satellite record** (MODIS/061 MOD13A3 monthly NDVI) evaluated against ERA5-Land climate reanalysis. Restricting the analysis to one sensor avoids the inter-instrument discontinuities that can bias trends assembled from spliced AVHRR–MODIS composites.

Six complementary methods are applied so that each conclusion rests on convergent rather than single-method evidence:

| Analysis | Purpose |
|---|---|
| GLS with AR(1) errors | Regional trends, accounting for serial correlation |
| Mann–Kendall / Theil–Sen | Distribution-free pixel-wise trend detection |
| Pixel-wise Pearson correlation | Spatial NDVI–climate association |
| Climate-adjusted GLS-AR(1) | Conditional effects; does the trend survive climate adjustment? |
| Endpoint & split-period tests | Temporal stability |
| Morlet wavelet coherence | Time–frequency coupling (continuous / episodic) |

## Key results

- Regional NDVI increased significantly in all three aggregations: **0.032** (annual), **0.040** (dry), **0.025** (wet) NDVI decade⁻¹, all *p* < 0.001.
- Significant positive pixel-wise trends covered **99.63 %** (annual), **99.34 %** (dry), and **97.46 %** (wet) of the study area. No cell showed a significant decrease.
- Climate associations are **seasonally structured**: temperature dominates the annual (62.76 %) and wet-season (53.08 %) maps, while the dry season is dominated by negative NDVI–radiation correlations (61.09 %) with a smaller positive precipitation footprint (27.81 %).
- The temporal trend **persists after conditioning on all three climate covariates** (β = 0.85–0.97; all *p* < 10⁻⁶). Only dry-season solar radiation contributes independently (β = −0.165, *p* = 0.025).
- Wavelet coherence is **intermittent, not continuous** — 12.15 % (radiation), 9.47 % (precipitation), 3.12 % (temperature) of valid time–scale cells.

Together these indicate greening that is spatially pervasive, temporally stable, and not reducible to the interannual variability of the measured climate variables — a pattern consistent with a substantial land-management contribution.

## Important: dry-season period

Annual and wet-season analyses span **2001–2025**. Dry-season analyses use the **24 complete October–March seasons of 2002–2025**. The 2001 dry composite would contain only January–March and is therefore incomplete; it is excluded from every dry-season statistic, map, correlation, and figure. See [`docs/DRY_SEASON_CORRECTION.md`](docs/DRY_SEASON_CORRECTION.md).

ERA5-Land surface net solar radiation is converted from J m⁻² to **MJ m⁻²** (÷ 10⁶) after monthly aggregation; annual and seasonal values are summed energy totals.

## Repository layout

```
.
├── data/tables/          # All numeric results as CSV (42 files, ~250 KB)
│   ├── main/                     # annual + wet-season outputs (2001–2025)
│   ├── corrected_dry_2002_2025/  # corrected dry-season outputs
│   ├── wavelet/                  # coherence settings, summary, monthly inputs
│   └── robustness/               # endpoint, split-period, land-cover checks
├── figures/
│   ├── main/             # Fig. 1–8 as published
│   └── supplementary/    # Fig. S1
├── scripts/
│   ├── python/           # trend statistics, correlation areas, GLS, figures
│   └── R/                # dry-season correction pipeline
├── docs/                 # methods, data sources, correction note
└── manuscript/           # manuscript and supplementary (.docx)
```

Raster data (GeoTIFF) are **not** tracked here because of size; see [`docs/DATA.md`](docs/DATA.md) for how to regenerate or obtain them.

## Quick start

```bash
git clone https://github.com/<your-username>/southern-china-greening-2001-2025.git
cd southern-china-greening-2001-2025
pip install -r requirements.txt

cd scripts/python

# Climate-adjusted GLS-AR(1) + split-period robustness
python climate_adjusted_gls.py \
  --means ../../data/tables/corrected_dry_2002_2025/modis_era5land_regional_means_complete_dry_2002_2025.csv \
  --out-gls ../../data/tables/corrected_dry_2002_2025/gls_climate_adjusted.csv \
  --out-split ../../data/tables/corrected_dry_2002_2025/split_period.csv

# Regenerate figures
python make_figures.py \
  --endpoints ../../data/tables/corrected_dry_2002_2025/Table_MODIS_NDVI_Endpoint_Sensitivity_CompleteDry_2002_2025.csv \
  --gls ../../data/tables/corrected_dry_2002_2025/gls_climate_adjusted.csv \
  --outdir ../../figures

# Correlation-class areas (requires the GeoTIFF rasters)
python correlation_area_summary.py --root /path/to/project --csv correlation_areas.csv
```

`scripts/python/trend_stats.py` is a dependency-light module (numpy only) providing Mann–Kendall, Theil–Sen, and Prais–Winsten GLS-AR(1); run it directly for a self-test.

## Reproducibility note

`correlation_area_summary.py` computes latitude-weighted cell areas analytically and **reproduces the published annual values exactly** (precipitation 7.06 / 0.56 %, radiation 0.00 / 24.30 %), which validates the area weighting.

The climate-adjusted GLS coefficients here are fitted by **Prais–Winsten** estimation. Values obtained with R's `nlme::gls(correlation = corAR1())` use maximum likelihood and may differ in the second or third decimal; the conclusions are unaffected.

## Data sources

| Dataset | Product | Access |
|---|---|---|
| NDVI | MODIS/061 MOD13A3 (monthly, 1 km) | [NASA LP DAAC](https://doi.org/10.5067/MODIS/MOD13A3.006) |
| Climate | ECMWF ERA5-Land monthly aggregated | [Copernicus CDS](https://doi.org/10.5194/essd-13-4349-2021) |

Both are freely available. Analysis grid: 0.05° geographic, clipped to `SouthernChinaD.shp`.

## Citation

If you use this code or these data, please cite the manuscript (see [`CITATION.cff`](CITATION.cff)) and the underlying datasets.

## License

Code is released under the [MIT License](LICENSE). Derived data tables and figures are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Source MODIS and ERA5-Land data remain subject to their providers' terms.

## Contact


