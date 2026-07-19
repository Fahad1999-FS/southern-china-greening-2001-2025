# Data

## What is in this repository

**Included** — all numeric results as CSV (`data/tables/`, ~250 KB) and
publication figures as PNG (`figures/`). These are sufficient to reproduce every
table, statistic, and figure in the manuscript except the pixel-wise correlation
area summary, which needs the rasters.

**Not included** — GeoTIFF rasters and analysis-ready stacks (several hundred MB),
excluded via `.gitignore`.

## Source datasets

| Dataset | Product | DOI / access |
|---|---|---|
| NDVI | MODIS/061 MOD13A3, monthly, 1 km | https://doi.org/10.5067/MODIS/MOD13A3.006 |
| Climate | ECMWF ERA5-Land monthly aggregated | https://doi.org/10.5194/essd-13-4349-2021 |

Both are freely available. In this study they were accessed through Google Earth
Engine (`MODIS/061/MOD13A3`, `ECMWF/ERA5_LAND/MONTHLY_AGGR`), clipped to the
study boundary, and exported on a common 0.05° geographic grid.

## Study boundary

`SouthernChinaD.shp` defines the analysis mask, applied uniformly to all raster
clipping, area summaries, regional means and pixel-wise mapping. Extent
approximately 97.5–122° E, 18–32° N, covering Yunnan, Guizhou, Guangxi, Hunan,
Jiangxi, Fujian, Guangdong, Zhejiang, Hainan and Shanghai.

## Table directory guide

| Folder | Contents |
|---|---|
| `main/` | Annual and wet-season trends, regional means, MK/Sen area summaries, correlation summaries (2001–2025) |
| `corrected_dry_2002_2025/` | Corrected dry-season trends, regional means, MK/Sen areas, endpoint sensitivity |
| `wavelet/` | Coherence settings, Monte Carlo summary, monthly input series (raw, anomaly, detrended-standardised) |
| `robustness/` | Endpoint, split-period, climate-adjusted GLS, land-cover stratification, Terra–Aqua agreement, lagged associations |

## Regenerating rasters

Rasters can be regenerated from the source collections using the export settings
documented in `docs/METHODS.md` (0.05° grid, study-area mask, seasonal
compositing as defined). The correlation classification rasters
(`*_CleanClass_*.tif`) assign class 4 to non-significant cells and classes
1–2 / 6–7 to significant negative / positive correlations.

Contact the corresponding author for direct access to processed rasters.
