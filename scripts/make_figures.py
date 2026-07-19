"""
make_figures.py
Reproduces the endpoint-sensitivity figure (main-text Fig. 7) and the
climate-adjusted GLS coefficient figure (Fig. S1) in the journal style used
throughout the manuscript: muted academic palette, bold black typography,
panel labels A/B/C, and no in-figure title (titles live in the captions).

Usage
-----
    python make_figures.py --endpoints <endpoint_sensitivity.csv> \
                           --gls <gls_coefficients.csv> \
                           --outdir ../../figures

Author: Kaleem Mehmood
License: MIT
"""

from __future__ import annotations
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams

# --- shared journal style -------------------------------------------------
rcParams.update({
    "font.family": "DejaVu Sans",
    "font.weight": "bold",
    "axes.labelweight": "bold",
    "axes.edgecolor": "black",
    "axes.linewidth": 1.3,
    "text.color": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
})

SEASON_COLOUR = {"Annual": "#2F5C8A", "Dry": "#B4652A", "Wet": "#2E7D5B"}
PANEL = {"Annual": "A", "Dry": "B", "Wet": "C"}
SEASONS = ["Annual", "Dry", "Wet"]
DPI = 400


def _style_axes(ax, season):
    ax.grid(axis="y", ls=":", lw=0.9, color="0.55", alpha=0.9, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", which="major", length=5, width=1.3)
    for sp in ax.spines.values():
        sp.set_linewidth(1.3)
        sp.set_color("black")
    ax.text(-0.02, 1.045, PANEL[season], transform=ax.transAxes,
            fontsize=18, fontweight="bold", color="black", va="bottom", ha="left")
    ax.text(0.055, 1.05, season, transform=ax.transAxes,
            fontsize=13.5, fontweight="bold", color="black", va="bottom", ha="left")


def figure_endpoint(csv: Path, out: Path):
    """Main-text Fig. 7: Theil-Sen slope for every endpoint window."""
    d = pd.read_csv(csv)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.4))

    for ax, season in zip(axes, SEASONS):
        sub = d[d["season"] == season]
        labels = [f"{int(a)}–{int(b)}"
                  for a, b in zip(sub["start_year"], sub["end_year"])]
        x = np.arange(len(sub))
        y = sub["slope_per_decade_Sen"].values

        ax.bar(x, y, color=SEASON_COLOUR[season], alpha=0.92, edgecolor="black",
               linewidth=1.1, width=0.72, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=55, ha="right",
                           fontsize=10.5, fontweight="bold")
        for t in ax.get_yticklabels():
            t.set_fontsize(11.5)
            t.set_fontweight("bold")
        ax.set_xlabel("Analysis window (start–end year)",
                      fontsize=12.5, fontweight="bold", labelpad=8)
        if season == "Annual":
            ax.set_ylabel("Theil–Sen slope (NDVI decade$^{-1}$)",
                          fontsize=13, fontweight="bold", labelpad=8)
        ax.set_ylim(0, y.max() * 1.22)
        for xi, yi in zip(x, y):
            ax.text(xi, yi + y.max() * 0.02, f"{yi:.3f}", ha="center", va="bottom",
                    fontsize=10, fontweight="bold", color="black",
                    rotation=90, zorder=4)
        _style_axes(ax, season)

    fig.tight_layout(w_pad=2.6)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Written: {out}")


def figure_gls(csv: Path, out: Path):
    """Fig. S1: standardised GLS-AR(1) coefficients with 95% CIs."""
    d = pd.read_csv(csv)
    terms = ["Year", "Precipitation", "Temperature", "Solar radiation"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8))

    for ax, season in zip(axes, SEASONS):
        sub = d[(d["season"] == season) & (d["term"].isin(terms))]
        sub = sub.set_index("term").loc[terms].reset_index()
        beta = sub["beta_std"].values
        se = sub["se"].values
        pvals = sub["p"].values
        y = np.arange(len(sub))[::-1]

        ax.errorbar(beta, y, xerr=1.96 * se, fmt="o", color=SEASON_COLOUR[season],
                    ecolor="black", elinewidth=1.6, capsize=5, capthick=1.6,
                    markersize=10, markeredgecolor="black", markeredgewidth=1.2,
                    zorder=3)
        ax.axvline(0, color="black", lw=1.4, ls="--", zorder=2)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["term"], fontsize=11.5, fontweight="bold")
        for t in ax.get_xticklabels():
            t.set_fontsize(11)
            t.set_fontweight("bold")
        ax.set_xlabel("Standardised coefficient (β)",
                      fontsize=12.5, fontweight="bold", labelpad=8)
        for yi, bi, pi in zip(y, beta, pvals):
            ax.text(bi, yi + 0.22, ("*" if pi < 0.05 else "") + f" {bi:.3f}",
                    ha="center", va="bottom", fontsize=10,
                    fontweight="bold", color="black")
        ax.grid(axis="x", ls=":", lw=0.9, color="0.55")
        ax.set_axisbelow(True)
        ax.tick_params(length=5, width=1.3)
        ax.set_ylim(-0.6, len(sub) - 0.2)
        for sp in ax.spines.values():
            sp.set_linewidth(1.3)
        ax.text(-0.02, 1.05, PANEL[season], transform=ax.transAxes,
                fontsize=18, fontweight="bold", va="bottom")
        ax.text(0.07, 1.06, season, transform=ax.transAxes,
                fontsize=13.5, fontweight="bold", va="bottom")

    fig.tight_layout(w_pad=2.4)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Written: {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--endpoints", type=Path, default=None,
                    help="Endpoint sensitivity CSV")
    ap.add_argument("--gls", type=Path, default=None,
                    help="GLS coefficient CSV from climate_adjusted_gls.py")
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    if args.endpoints:
        figure_endpoint(args.endpoints, args.outdir / "main" / "Fig07_endpoint_sensitivity.png")
    if args.gls:
        figure_gls(args.gls, args.outdir / "supplementary" / "FigS1_gls_coefficients.png")
    if not args.endpoints and not args.gls:
        ap.error("supply --endpoints and/or --gls")


if __name__ == "__main__":
    main()
