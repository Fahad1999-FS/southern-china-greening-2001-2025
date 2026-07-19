# R scripts

## `run_complete_dry_season_correction_2002_2025.R`

Regenerates all dry-season outputs using only the 24 complete October–March
seasons of 2002–2025, and confirms the ERA5-Land radiation conversion
(J m⁻² → MJ m⁻², ÷ 10⁶ after monthly aggregation).

Produces: corrected regional means and trend tables, pixel-wise MK/Sen maps and
area summaries, dry correlation rasters, endpoint sensitivity, and corrected
figures.

> **Place the script here.** It lives outside the tracked project directory in
> the original working tree; copy it into this folder before publishing the
> repository.

### Requirements

```r
install.packages(c("terra", "nlme", "trend", "dplyr", "readr", "ggplot2"))
```

### Run

```bash
Rscript run_complete_dry_season_correction_2002_2025.R
```

Adjust the input/output paths at the head of the script to match your local
project layout.

## Relationship to the Python scripts

The R pipeline is authoritative for the manuscript's regional trend statistics.
The Python scripts in `../python/` provide an independent reimplementation used
for the correlation-area summaries, the climate-adjusted models, the
split-period recomputation, and figure generation.
