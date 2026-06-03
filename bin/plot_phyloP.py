#!/usr/bin/env python3
"""
plot_phylop.py — Visualise max_phylop / median_phylop distributions and
scatter-plot them against genomic / expression features.

Outputs
-------
  <outdir>/max_phylop_dist.png
  <outdir>/med_phylop_dist.png
  <outdir>/scatter_max_phylop.png   (one panel per feature, 4-wide layout)
  <outdir>/scatter_med_phylop.png

Usage
-----
python plot_phylop.py data.csv \\
    -c gene_id \\
    -g BRCA1 TP53 MYC \\
    -o ./plots
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ── per-feature colour palette (accent, light-fill) ──────────────────────────
PALETTE = {
    "utr3_length":     ("#2563EB", "#BFDBFE"),   # blue
    "utr5_length":     ("#16A34A", "#BBF7D0"),   # green
    "utr3_gc":         ("#DC2626", "#FEE2E2"),   # red
    "utr5_gc":         ("#7C3AED", "#EDE9FE"),   # purple
    "human_half_life": ("#2563EB", "#BFDBFE"),
    "half_life":       ("#2563EB", "#BFDBFE"),
    "mouse_half_life": ("#16A34A", "#BBF7D0"),
    "human_te":        ("#DC2626", "#FEE2E2"),
    "mean_te":         ("#DC2626", "#FEE2E2"),
    "mouse_te":        ("#7C3AED", "#EDE9FE"),
}

HIGHLIGHT_ACCENT = "#F59E0B"   # amber — genes of interest
HIGHLIGHT_LIGHT  = "#FEF3C7"

FEATURE_LABELS = {
    "utr3_length":     "3\u2032 UTR Length (nt)",
    "utr5_length":     "5\u2032 UTR Length (nt)",
    "utr3_gc":         "3\u2032 UTR GC Content",
    "utr5_gc":         "5\u2032 UTR GC Content",
    "human_half_life": "Half-Life\nHuman (z-score)",
    "half_life":       "Half-Life\nHuman (z-score)",
    "mouse_half_life": "Half-Life\nMouse (z-score)",
    "human_te":        "Translation Efficiency\nHuman (z-score)",
    "mean_te":         "Translation Efficiency\nHuman (z-score)",
    "mouse_te":        "Translation Efficiency\nMouse (z-score)",
}

PHYLOP_COLS = ["max_phylop", "median_phylop"]
PHYLOP_NICE = {"max_phylop": "Max PhyloP", "median_phylop": "Median PhyloP"}


# ── helpers ───────────────────────────────────────────────────────────────────

def load_csv(path):
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    return df

def present_features(df):
    return {k: v for k, v in FEATURE_LABELS.items() if k in df.columns}

def present_phylop(df):
    return [c for c in PHYLOP_COLS if c in df.columns]

def highlight_mask(df, gene_col, genes):
    if not gene_col or not genes or gene_col not in df.columns:
        return None
    return df[gene_col].astype(str).isin([str(g) for g in genes])

def _style_ax(ax):
    """Strip top/right spines, subtle tick params — matches reference style."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)

def _sig_stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


# ── distribution plots ────────────────────────────────────────────────────────

def plot_distribution(df, pcol, mask, gene_col, outdir):
    """Histogram + KDE for one phylop column → <pcol>_dist.png"""
    data = df[pcol].dropna()
    nice = PHYLOP_NICE[pcol]
    accent, light = "#2563EB", "#BFDBFE"

    fig, ax = plt.subplots(figsize=(7, 5))
    fig.patch.set_facecolor("white")

    # histogram
    ax.hist(data, bins=60, color=light, edgecolor=accent,
            linewidth=0.4, density=True, zorder=2, alpha=0.85)

    # KDE
    kde_x = np.linspace(data.min(), data.max(), 500)
    kde   = stats.gaussian_kde(data)
    ax.plot(kde_x, kde(kde_x), color=accent, linewidth=2, zorder=3)

    # median line
    med = data.median()
    ax.axvline(med, color="#DC2626", linewidth=1.8, linestyle="--",
               label=f"Median = {med:.3f}", zorder=4)

    # highlighted genes — vertical lines + rotated labels
    if mask is not None:
        hi_idx = df.index[mask & df[pcol].notna()]
        # draw lines first so we know the y-scale
        for idx in hi_idx:
            v = df.at[idx, pcol]
            ax.axvline(v, color=HIGHLIGHT_ACCENT, linewidth=1.3,
                       alpha=0.9, zorder=5)
        # re-draw after autoscale settles
        ymax = ax.get_ylim()[1]
        for idx in hi_idx:
            v   = df.at[idx, pcol]
            lbl = str(df.at[idx, gene_col]) if gene_col else ""
            ax.text(v, ymax * 0.97, lbl,
                    rotation=90, va="top", ha="right",
                    fontsize=7, color=HIGHLIGHT_ACCENT,
                    fontweight="semibold", zorder=6)

    _style_ax(ax)
    ax.set_xlabel(nice, fontsize=11)
    ax.set_ylabel("Density", fontsize=11)
    ax.set_title(f"Distribution of {nice}",
                 fontsize=13, fontweight="bold", pad=10)
    ax.legend(fontsize=9, frameon=True, framealpha=0.85)

    fig.tight_layout()
    out = os.path.join(outdir, f"{pcol}_dist.png")
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


# ── scatter grid ──────────────────────────────────────────────────────────────

def plot_scatter_grid(df, pcol, features, mask, gene_col, outdir):
    """
    One figure, 4-panels-wide layout (like reference), one panel per feature.
    Rows wrap automatically.  → scatter_<pcol>.png
    """
    feat_keys = list(features.keys())
    n     = len(feat_keys)
    if n == 0:
        return
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(5 * ncols, 5 * nrows),
                             squeeze=False)
    fig.patch.set_facecolor("white")
    nice = PHYLOP_NICE[pcol]

    for i, feat in enumerate(feat_keys):
        r, c  = divmod(i, ncols)
        ax    = axes[r][c]
        flabel = FEATURE_LABELS[feat]
        accent, light = PALETTE.get(feat, ("#2563EB", "#BFDBFE"))

        sub = df[[pcol, feat]].dropna()
        x   = sub[feat].values
        y   = sub[pcol].values

        # scatter
        ax.scatter(x, y, color=light, edgecolors=accent, linewidths=0.4,
                   alpha=0.35, s=18, rasterized=True, zorder=2)

        # regression trend line
        if len(sub) > 2:
            slope, intercept, r_val, p_val, _ = stats.linregress(x, y)
            x_line = np.linspace(x.min(), x.max(), 300)
            ax.plot(x_line, slope * x_line + intercept,
                    color=accent, linewidth=2, zorder=3)
            r2  = r_val ** 2
            sig = _sig_stars(p_val)
            ax.text(0.05, 0.95,
                    f"R² = {r2:.3f}\np = {p_val:.2e} {sig}",
                    transform=ax.transAxes, fontsize=8.5,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                              edgecolor=accent, alpha=0.85))

        # highlighted genes
        if mask is not None:
            hi = df.loc[mask, [feat, pcol]].dropna()
            ax.scatter(hi[feat], hi[pcol],
                       color=HIGHLIGHT_LIGHT, edgecolors=HIGHLIGHT_ACCENT,
                       linewidths=0.8, s=55, zorder=4)
            if gene_col:
                for idx in hi.index:
                    ax.annotate(
                        str(df.at[idx, gene_col]),
                        (df.at[idx, feat], df.at[idx, pcol]),
                        textcoords="offset points", xytext=(5, 4),
                        fontsize=7, color=HIGHLIGHT_ACCENT,
                        fontweight="semibold", zorder=5,
                    )

        # give the y-axis a bit of headroom so the dense top cluster
        # isn't clipped — add 15 % padding above the data max
        ymin, ymax_auto = ax.get_ylim()
        data_ymax = sub[pcol].max()
        padded_ymax = data_ymax + 0.15 * (data_ymax - sub[pcol].min())
        ax.set_ylim(ymin, max(ymax_auto, padded_ymax))

        _style_ax(ax)
        ax.set_xlabel(flabel, fontsize=10)
        ax.set_ylabel(nice, fontsize=10)
        ax.set_title(flabel.replace("\n", " "),
                     fontsize=11, fontweight="semibold")

    # hide unused panels
    for j in range(n, nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r][c].set_visible(False)

    fig.suptitle(f"{nice} — Feature Scatter Plots",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    out = os.path.join(outdir, f"scatter_{pcol}.png")
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot PhyloP distributions and feature scatter plots.")
    p.add_argument("csv", help="Input CSV file")
    p.add_argument("-c", "--gene_id_column", default=None,
                   help="Column name containing gene IDs")
    p.add_argument("-g", "--genes_of_interest", nargs="*", default=[],
                   metavar="GENE",
                   help="Gene IDs to highlight (space-separated)")
    p.add_argument("-o", "--outdir", default=".",
                   help="Output directory for PNGs (default: current dir)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = load_csv(args.csv)

    phylop_cols = present_phylop(df)
    if not phylop_cols:
        sys.exit("ERROR: CSV contains neither 'max_phylop' nor 'median_phylop'.")

    features = present_features(df)
    gene_col  = args.gene_id_column.lower() if args.gene_id_column else None
    if gene_col and gene_col not in df.columns:
        print(f"WARNING: column '{gene_col}' not found — highlighting disabled.",
              file=sys.stderr)
        gene_col = None

    mask = highlight_mask(df, gene_col, args.genes_of_interest)
    if args.genes_of_interest and mask is not None:
        print(f"Highlighting {mask.sum()} / {len(args.genes_of_interest)} "
              "requested genes.")

    for pcol in phylop_cols:
        plot_distribution(df, pcol, mask, gene_col, args.outdir)
        if features:
            plot_scatter_grid(df, pcol, features, mask, gene_col, args.outdir)
        else:
            print("No recognised feature columns found — skipping scatter plots.")


if __name__ == "__main__":
    main()
