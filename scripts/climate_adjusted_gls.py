"""
climate_adjusted_gls.py
Climate-adjusted GLS-AR(1) models and split-period trend robustness for the
regional MODIS NDVI series.

Fits, for each seasonal aggregation:

    NDVI* = b0 + b1*Year* + b2*PPT* + b3*Temp* + b4*Sradi* + eps,   eps ~ AR(1)

where * denotes standardisation to zero mean and unit variance, and compares
model support against a year-only model via AIC. Also recomputes the
split-period Theil-Sen / Mann-Kendall trends.

Note
----
Coefficients are obtained by Prais-Winsten estimation of the AR(1) structure.
Values fitted with R's nlme::gls(correlation = corAR1()) use maximum
likelihood and may differ in the second or third decimal; conclusions are
unaffected.

Usage
-----
    python climate_adjusted_gls.py \
        --means ../../data/tables/corrected_dry_2002_2025/modis_era5land_regional_means_complete_dry_2002_2025.csv \
        --out-gls gls.csv --out-split split.csv

Author: Kaleem Mehmood
License: MIT
"""

from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from trend_stats import gls_ar1, mann_kendall, theil_sen, zscore

SEASONS = ["Annual", "Dry", "Wet"]
TERMS = ["(Intercept)", "Year", "Precipitation", "Temperature", "Solar radiation"]
SPLIT_BOUNDARY = 2012          # first period ends here, second starts at +1


def fit_season(d: pd.DataFrame) -> pd.DataFrame:
    """Climate-adjusted and year-only GLS-AR(1) for one seasonal series."""
    y = zscore(d["NDVI"].values)
    year_z = zscore(d["year"].values)

    X_clim = np.column_stack([
        np.ones(len(d)), year_z,
        zscore(d["PPT"].values), zscore(d["Temp"].values), zscore(d["Sradi"].values),
    ])
    X_year = np.column_stack([np.ones(len(d)), year_z])

    full = gls_ar1(y, X_clim)
    base = gls_ar1(y, X_year)

    return pd.DataFrame({
        "term": TERMS,
        "beta_std": full["beta"],
        "se": full["se"],
        "p": full["p"],
        "phi": full["phi"],
        "aic_year_only": base["aic"],
        "aic_year_climate": full["aic"],
        "delta_aic": full["aic"] - base["aic"],
    })


def split_period(d: pd.DataFrame) -> pd.DataFrame:
    """Theil-Sen / Mann-Kendall trends for each half of the record."""
    rows = []
    start = int(d["year"].min())
    for lo, hi in [(start, SPLIT_BOUNDARY), (SPLIT_BOUNDARY + 1, int(d["year"].max()))]:
        sub = d[(d["year"] >= lo) & (d["year"] <= hi)]
        if len(sub) < 3:
            continue
        mk = mann_kendall(sub["NDVI"].values)
        rows.append({
            "period": f"{lo}-{hi}", "n": len(sub),
            "sen_slope_per_decade": theil_sen(sub["NDVI"].values,
                                              sub["year"].values.astype(float)) * 10.0,
            "p_MK": mk["p"], "tau_MK": mk["tau"],
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--means", required=True, type=Path,
                    help="Regional means CSV (year, season, NDVI, PPT, Temp, Sradi)")
    ap.add_argument("--out-gls", type=Path, default=None)
    ap.add_argument("--out-split", type=Path, default=None)
    args = ap.parse_args()

    df = pd.read_csv(args.means)

    gls_all, split_all = [], []
    for season in SEASONS:
        d = df[df["season"] == season].sort_values("year").reset_index(drop=True)
        period = f"{int(d.year.min())}-{int(d.year.max())}"

        g = fit_season(d)
        g.insert(0, "period", period)
        g.insert(0, "season", season)
        gls_all.append(g)

        print(f"\n=== {season}  n={len(d)}  ({period})  phi={g['phi'].iloc[0]:.3f}  "
              f"dAIC={g['delta_aic'].iloc[0]:+.2f}")
        for _, r in g.iterrows():
            if r["term"] == "(Intercept)":
                continue
            star = " *" if r["p"] < 0.05 else ""
            print(f"    {r['term']:16s} beta={r['beta_std']:+.3f}  "
                  f"se={r['se']:.3f}  p={r['p']:.3g}{star}")

        s = split_period(d)
        s.insert(0, "season", season)
        split_all.append(s)

    gls_df = pd.concat(gls_all, ignore_index=True)
    split_df = pd.concat(split_all, ignore_index=True)

    print("\n=== Split-period robustness ===")
    print(split_df.to_string(index=False,
                             float_format=lambda v: f"{v:.4g}"))

    if args.out_gls:
        args.out_gls.parent.mkdir(parents=True, exist_ok=True)
        gls_df.to_csv(args.out_gls, index=False)
        print(f"\nWritten: {args.out_gls}")
    if args.out_split:
        args.out_split.parent.mkdir(parents=True, exist_ok=True)
        split_df.to_csv(args.out_split, index=False)
        print(f"Written: {args.out_split}")


if __name__ == "__main__":
    main()
