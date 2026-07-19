# Methods

Condensed technical reference for the analysis. Equation numbers follow the manuscript.

---

## 1. Data

| Variable | Source | Native resolution | Units used |
|---|---|---|---|
| NDVI | MODIS/061 MOD13A3 (monthly) | 1 km | index |
| Precipitation (PPT) | ERA5-Land monthly aggregated | ~9 km | mm |
| 2 m temperature (Temp) | ERA5-Land monthly aggregated | ~9 km | °C |
| Surface net solar radiation (Sradi) | ERA5-Land monthly aggregated | ~9 km | MJ m⁻² |

All fields exported to a common **0.05° geographic grid** and clipped to `SouthernChinaD.shp`. The shared grid supports cell-by-cell comparison but does **not** increase the native information content of the climate fields — pixel-level correlations should be interpreted at the ERA5-Land scale.

Radiation is converted J m⁻² → MJ m⁻² by division by 10⁶ **after** monthly aggregation; annual and seasonal values are summed energy totals.

## 2. Seasonal compositing

| Aggregation | Months | Years |
|---|---|---|
| Annual | January–December | 2001–2025 |
| Dry | October–December (prior year) + January–March | **2002–2025** |
| Wet | April–September | 2001–2025 |

Composite value for cell and year *y*:

```
V̄(y) = (1/M) Σ_{k=1..M} v(y, k)                                    (1)
```

where *M* is the number of contributing months. The 2001 dry composite is incomplete (January–March only) and is excluded — see `DRY_SEASON_CORRECTION.md`.

## 3. Regional trend estimation

Interannual series exhibit temporal autocorrelation, which inflates apparent significance under OLS. Primary estimator is GLS with AR(1) errors:

```
V̄(t) = β₀ + β₁t + ε(t),     ε(t) = φ ε(t−1) + η(t)                  (2)
β̂ = (XᵀΣ⁻¹X)⁻¹ XᵀΣ⁻¹y                                              (3)
Σᵢⱼ = σ²/(1−φ²) · φ^|i−j|                                            (4)
```

Fitted by **Prais–Winsten** iteration (retains the first observation). Implemented in `scripts/python/trend_stats.py::gls_ar1`.

## 4. Non-parametric trend detection

Mann–Kendall statistic with tie correction, and the Theil–Sen slope:

```
S = Σ_{i<j} sgn(x_j − x_i)                                           (5)
Var(S) = [n(n−1)(2n+5) − Σ_p t_p(t_p−1)(2t_p+5)] / 18                (6)
Z = (S−1)/√Var(S) if S>0;  0 if S=0;  (S+1)/√Var(S) if S<0           (7)
τ = 2S / [n(n−1)]                                                    (8)
β_Sen = median{ (x_j − x_i)/(j − i) },  ∀ i<j                         (9)
```

A cell is classified as significantly increasing/decreasing where Mann–Kendall *p* < 0.05 and the Theil–Sen slope is positive/negative.

Implemented in `trend_stats.py::mann_kendall` and `::theil_sen`.

## 5. Spatial climate–NDVI association

Pearson correlation per valid cell, with a Student-*t* significance test:

```
r = Σ(xᵢ−x̄)(yᵢ−ȳ) / √[ Σ(xᵢ−x̄)² Σ(yᵢ−ȳ)² ]                        (10)
t = r √[ (n−2)/(1−r²) ]                                              (11)
```

Correlations are classified by sign, strength and significance, then summarised as **latitude-weighted** percentages of the study area. Cell area on a geographic grid:

```
A(φ) = R² · Δλ · [ sin(φ_top) − sin(φ_bottom) ],   R = 6371.0072 km
```

Implemented in `scripts/python/correlation_area_summary.py`. This reproduces the published annual values exactly, validating the weighting.

> These area fractions describe **spatial extent**, not field significance — they are not corrected for the multiplicity of pixel-level tests (Wilks 2016).

## 6. Climate-adjusted regression

Each variable standardised, then a multiple GLS-AR(1) model:

```
z = (x − μ)/σ                                                        (12)
V̄* = β₀ + β₁Year* + β₂PPT* + β₃Temp* + β₄Sradi* + ε(t)              (13)
```

Model support compared against a year-only model by AIC. Standardised coefficients quantify **conditional covariability**; they do not identify causal origin.

Implemented in `scripts/python/climate_adjusted_gls.py`.

## 7. Temporal robustness

- **Endpoint sensitivity** — recompute Theil–Sen / Mann–Kendall after dropping up to two years from the start, end, or both (windows of 21–25 years; 20–24 for the dry season).
- **Split-period** — compare trends before/after 2012.

## 8. Wavelet coherence

Monthly series (300 observations) deseasonalised, detrended, standardised:

```
x′(y,m) = x(y,m) − (1/N) Σ_{y′} x(y′,m)                             (14)
ψ₀(η) = π^(−1/4) e^{iω₀η} e^{−η²/2},   ω₀ = 6                       (15)
R²(s) = |S(s⁻¹ W^XY(s))|² / [ S(s⁻¹|W^X(s)|²) S(s⁻¹|W^Y(s)|²) ]      (16)
```

Periods 3–48 months, scale resolution *dj* = 0.05, significance against a red-noise (AR(1)) background from 1,000 Monte Carlo simulations at *p* < 0.05. The cone of influence is excluded. Computed with `pycwt`.

Because this operates on continuous monthly series rather than seasonal composites, it is **unaffected** by the dry-season compositing correction.

## 9. Reporting conventions

- Significance at *p* < 0.05 unless stated otherwise.
- Regional trends: GLS-AR(1) slope per decade with 95 % CI, alongside Theil–Sen slope and Mann–Kendall *p*.
- Pixel-wise results: areas and area percentages by class, not single-cell coefficients.
- Climate-adjusted models: standardised coefficients, AIC relative to a year-only baseline.
- Wavelet: percentage of valid out-of-cone time–scale cells above threshold, plus mean coherence over those cells.

## References for methods

Grinsted et al. (2004) · Hamed & Rao (1998) · Kendall (1975) · Mann (1945) · Sen (1968) · Theil (1950) · Torrence & Compo (1998) · Wilks (2016). Full details in the manuscript reference list.
