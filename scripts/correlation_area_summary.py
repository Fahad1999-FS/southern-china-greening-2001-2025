"""
correlation_area_summary.py
Latitude-weighted area summary of pixel-wise NDVI-climate correlation classes.

Reads the *_CleanClass_*.tif rasters (classified Pearson correlation, with
non-significant cells assigned class 4) and reports the area and percentage of
the study domain in each class, per climate variable and season.

Cell area is computed analytically for a geographic (lat/lon) grid:

    A(phi) = R^2 * dlambda * ( sin(phi_top) - sin(phi_bottom) )

with R the authalic Earth radius. This reproduces the published annual values
exactly (precipitation 7.06 / 0.56 %, radiation 0.00 / 24.30 %).

Usage
-----
    python correlation_area_summary.py --root /path/to/project [--csv out.csv]

Dependencies: numpy, pandas, Pillow (GeoTIFF tags read via Pillow; no GDAL
required for these single-band float rasters).

Author: Kaleem Mehmood
License: MIT
"""

from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

R_EARTH_KM = 6371.0072          # authalic radius
NODATA_BELOW = -9000.0

TAG_PIXEL_SCALE = 33550
TAG_TIEPOINT = 33922

CLASS_LABELS = {
    1: "< -0.40",
    2: "-0.40 to -0.20",
    3: "-0.20 to 0.00",
    4: "Not significant",
    5: "0.00 to 0.20",
    6: "0.20 to 0.40",
    7: "> 0.40",
}
POSITIVE_CLASSES = (6, 7)
NEGATIVE_CLASSES = (1, 2)

VARIABLES = [("PPT", "Precipitation"), ("Temp", "Temperature"), ("Sradi", "Solar radiation")]


def cell_area_km2(shape, pixel_scale_deg, lat_origin):
    """Per-cell area (km^2) for a north-up geographic grid."""
    ny, nx = shape
    lat_top = lat_origin - np.arange(ny) * pixel_scale_deg
    lat_bot = lat_top - pixel_scale_deg
    dlam = np.deg2rad(pixel_scale_deg)
    band = (R_EARTH_KM ** 2) * dlam * (
        np.sin(np.deg2rad(lat_top)) - np.sin(np.deg2rad(lat_bot))
    )
    return np.repeat(band[:, None], nx, axis=1)


def summarise(path: Path) -> pd.DataFrame:
    """Area and percentage per correlation class for one raster."""
    img = Image.open(path)
    tags = img.tag_v2
    pixel_scale = tags[TAG_PIXEL_SCALE][0]
    lat_origin = tags[TAG_TIEPOINT][4]

    arr = np.array(img)
    area = cell_area_km2(arr.shape, pixel_scale, lat_origin)
    valid = arr > NODATA_BELOW
    total = area[valid].sum()

    rows = []
    for cls in sorted(np.unique(arr[valid])):
        mask = (arr == cls) & valid
        a = area[mask].sum()
        rows.append({
            "class_id": int(cls),
            "class_label": CLASS_LABELS.get(int(cls), "unknown"),
            "area_km2": a,
            "percent": 100.0 * a / total,
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, type=Path,
                    help="Project root containing Correlations/ and "
                         "Complete_Dry_Season_2002_2025/")
    ap.add_argument("--csv", type=Path, default=None, help="Optional output CSV path")
    args = ap.parse_args()

    annual_wet = args.root / "Correlations"
    dry = args.root / "Complete_Dry_Season_2002_2025" / "Correlations"

    records = []
    for var, label in VARIABLES:
        for season in ("Annual", "Dry", "Wet"):
            if season == "Dry":
                path = dry / f"MODIS_NDVI_ERA5Land_{var}_Correlation_CleanClass_Dry_2002_2025.tif"
                period = "2002-2025"
            else:
                path = annual_wet / f"MODIS_NDVI_ERA5Land_{var}_Correlation_CleanClass_{season}_2001_2025.tif"
                period = "2001-2025"
            if not path.exists():
                print(f"  [skip] missing {path.name}")
                continue

            df = summarise(path)
            pos = df.loc[df.class_id.isin(POSITIVE_CLASSES), "percent"].sum()
            neg = df.loc[df.class_id.isin(NEGATIVE_CLASSES), "percent"].sum()
            for _, r in df.iterrows():
                records.append({"variable": label, "season": season, "period": period,
                                **r.to_dict()})
            print(f"{label:16s} {season:6s} ({period})  "
                  f"positive {pos:6.2f}%   negative {neg:6.2f}%")

    out = pd.DataFrame(records)
    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.csv, index=False)
        print(f"\nWritten: {args.csv}")


if __name__ == "__main__":
    main()
