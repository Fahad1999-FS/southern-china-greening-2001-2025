# Dry-season correction (2002–2025)

## Issue

The dry-season composite spans **October–December of the preceding year plus
January–March of the focal year**. Because the MODIS record used here begins in
January 2001, the 2001 dry composite could only contain January–March — three
months instead of six.

Including this half-length composite distorted dry-season statistics. The effect
was most visible in surface solar radiation, which is a **summed energy total**:
the 2001 value was roughly half a normal season, producing a spurious apparent
increase of **+460 MJ m⁻² decade⁻¹** (95 % CI −498 to 1419).

## Resolution

All dry-season statistics, maps, correlations, and figures were recomputed using
only the **24 complete October–March seasons of 2002–2025**. Annual and
wet-season analyses are unaffected and retain the full 2001–2025 window.

ERA5-Land radiation was also confirmed to be correctly converted from J m⁻² to
MJ m⁻² (÷ 10⁶) after monthly aggregation.

## Effect on results

| Quantity | Before (2001–2025) | After (2002–2025) |
|---|---|---|
| Dry NDVI GLS-AR(1) slope | 0.050 decade⁻¹ (n = 25) | **0.040 decade⁻¹ (n = 24)** |
| Dry NDVI Theil–Sen slope | 0.050 decade⁻¹ | **0.047 decade⁻¹** |
| Dry significant increasing area | 99.55 % | **99.34 %** |
| Dry radiation trend | +460.48 MJ m⁻² dec⁻¹ | **−31.96 MJ m⁻² dec⁻¹ (p = 0.352)** |
| Dry temperature trend | 0.304 °C dec⁻¹ (p = 0.028) | **0.231 °C dec⁻¹ (p = 0.124, n.s.)** |
| Dry NDVI–precipitation (+) | 53.52 % | **27.81 %** |
| Dry NDVI–temperature (+) | 44.57 % | **10.29 %** |
| Dry NDVI–radiation (−) | 9.80 % positive | **61.09 % negative** |
| Dry precipitation (climate-adjusted) | β = 0.236, p < 0.001 | **β = 0.096, p = 0.115 (n.s.)** |
| Dry solar radiation (climate-adjusted) | n.s. | **β = −0.165, p = 0.025** |

Two conclusions changed materially:

1. **Dry-season warming is no longer significant.** Significant temperature
   increase is confined to the annual and wet-season aggregations.
2. **Dry-season precipitation is no longer an independent contributor.** In the
   climate-adjusted model, solar radiation is now the only climate variable with
   an independent effect, and it is negative.

The central findings are unchanged: greening remains spatially pervasive,
temporally stable, steepest in the dry season, and not reducible to the
interannual variability of the measured climate variables.

## Reproducing

The R pipeline that regenerates the corrected dry-season outputs is
`scripts/R/run_complete_dry_season_correction_2002_2025.R`. Downstream Python
recomputation (correlation areas, climate-adjusted GLS, split-period, figures)
is in `scripts/python/`.
