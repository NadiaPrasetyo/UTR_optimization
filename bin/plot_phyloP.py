#!/usr/bin/env python3
"""
plot_phylop.py — Visualise UTR-specific PhyloP distributions and
scatter-plot them against genomic / expression features.

Outputs (one set per UTR × metric combination)
-------
  <outdir>/utr5_max_phylop_dist.png
  <outdir>/utr5_med_phylop_dist.png
  <outdir>/utr3_max_phylop_dist.png
  <outdir>/utr3_med_phylop_dist.png
  <outdir>/scatter_utr5_max_phylop.png
  <outdir>/scatter_utr5_med_phylop.png
  <outdir>/scatter_utr3_max_phylop.png
  <outdir>/scatter_utr3_med_phylop.png

Usage
-----
python plot_phylop.py data.csv \\
    -c display_name \\
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

# ── UTR PhyloP columns and display names ─────────────────────────────────────

# Each entry: (column_name, nice_label, utr_tag)
UTR_PHYLOP_COLS = [
    ("utr5_max_phylop", "5′ UTR Max PhyloP",    "utr5"),
    ("utr5_med_phylop", "5′ UTR Median PhyloP", "utr5"),
    ("utr3_max_phylop", "3′ UTR Max PhyloP",    "utr3"),
    ("utr3_med_phylop", "3′ UTR Median PhyloP", "utr3"),
]

# Colour scheme: utr5 = blue family, utr3 = green family
UTR_PALETTE = {
    "utr5": ("#2563EB", "#BFDBFE"),   # blue  accent / light
    "utr3": ("#16A34A", "#BBF7D0"),   # green accent / light
}

# ── per-feature colour palette ────────────────────────────────────────────────
FEATURE_PALETTE = {
    "utr3_length":  ("#16A34A", "#BBF7D0"),   # green
    "utr5_length":  ("#2563EB", "#BFDBFE"),   # blue
    "utr3_gc":      ("#DC2626", "#FEE2E2"),   # red
    "utr5_gc":      ("#7C3AED", "#EDE9FE"),   # purple
    "cds_length":   ("#0891B2", "#CFFAFE"),   # cyan
    "cds_gc":       ("#D97706", "#FEF3C7"),   # amber
    "half_life":    ("#2563EB", "#BFDBFE"),
    "mean_te":      ("#DC2626", "#FEE2E2"),
}

FEATURE_LABELS = {
    "utr3_length": "3\u2032 UTR Length (nt)",
    "utr5_length": "5\u2032 UTR Length (nt)",
    "utr3_gc":     "3\u2032 UTR GC Content",
    "utr5_gc":     "5\u2032 UTR GC Content",
    "cds_length":  "CDS Length (nt)",
    "cds_gc":      "CDS GC Content",
    "half_life":   "Half-Life (z-score)",
    "mean_te":     "Translation Efficiency (z-score)",
}

HIGHLIGHT_ACCENT = "#F59E0B"   # amber
HIGHLIGHT_LIGHT  = "#FEF3C7"


# ── helpers ───────────────────────────────────────────────────────────────────

def load_csv(path):
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    return df

def present_utr_phylop(df):
    """Return list of (col, nice, utr_tag) tuples that actually exist in df."""
    return [(col, nice, tag)
            for col, nice, tag in UTR_PHYLOP_COLS
            if col in df.columns]

def present_features(df):
    """Return {col: label} for feature columns present in df."""
    return {k: v for k, v in FEATURE_LABELS.items() if k in df.columns}

def highlight_mask(df, gene_col, genes):
    if not gene_col or not genes or gene_col not in df.columns:
        return None
    return df[gene_col].astype(str).isin([str(g) for g in genes])

def _style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)

def _sig_stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


# ── distribution plots ────────────────────────────────────────────────────────

def plot_distribution(df, pcol, nice, utr_tag, mask, gene_col, outdir):
    """Histogram + KDE for one UTR PhyloP column → <pcol>_dist.png"""
    data = df[pcol].dropna()
    accent, light = UTR_PALETTE[utr_tag]

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

    # highlighted genes
    if mask is not None:
        hi_idx = df.index[mask & df[pcol].notna()]
        for idx in hi_idx:
            v = df.at[idx, pcol]
            ax.axvline(v, color=HIGHLIGHT_ACCENT, linewidth=1.3,
                       alpha=0.9, zorder=5)
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

def plot_scatter_grid(df, pcol, nice, utr_tag, features, mask, gene_col, outdir):
    """
    One figure, up to 4 panels wide, one panel per feature.
    → scatter_<pcol>.png
    """
    feat_keys = list(features.keys())
    n = len(feat_keys)
    if n == 0:
        return

    # For UTR-specific scatter plots, skip the matching-UTR structural features
    # when they would be trivially correlated (optional — remove filter if
    # you want all features regardless).
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))

    # Use the UTR palette for the trend line; per-feature colour for dots
    utr_accent, _ = UTR_PALETTE[utr_tag]

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(5 * ncols, 5 * nrows),
                             squeeze=False)
    fig.patch.set_facecolor("white")

    for i, feat in enumerate(feat_keys):
        r, c    = divmod(i, ncols)
        ax      = axes[r][c]
        flabel  = FEATURE_LABELS[feat]
        f_accent, f_light = FEATURE_PALETTE.get(feat, ("#2563EB", "#BFDBFE"))

        sub = df[[pcol, feat]].dropna()
        x   = sub[pcol].values   # UTR phylop on x-axis
        y   = sub[feat].values   # feature on y-axis

        # scatter (feature colour for the dots)
        ax.scatter(x, y, color=f_light, edgecolors=f_accent, linewidths=0.4,
                   alpha=0.35, s=18, rasterized=True, zorder=2)

        # regression trend line (UTR colour)
        if len(sub) > 2:
            slope, intercept, r_val, p_val, _ = stats.linregress(x, y)
            x_line = np.linspace(x.min(), x.max(), 300)
            ax.plot(x_line, slope * x_line + intercept,
                    color=utr_accent, linewidth=2, zorder=3)
            r2  = r_val ** 2
            sig = _sig_stars(p_val)
            ax.text(0.05, 0.95,
                    f"R² = {r2:.3f}\np = {p_val:.2e} {sig}",
                    transform=ax.transAxes, fontsize=8.5,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                              edgecolor=utr_accent, alpha=0.85))

        # highlighted genes
        if mask is not None:
            hi = df.loc[mask, [pcol, feat]].dropna()
            ax.scatter(hi[pcol], hi[feat],
                       color=HIGHLIGHT_LIGHT, edgecolors=HIGHLIGHT_ACCENT,
                       linewidths=0.8, s=55, zorder=4)
            if gene_col:
                for idx in hi.index:
                    ax.annotate(
                        str(df.at[idx, gene_col]),
                        (df.at[idx, pcol], df.at[idx, feat]),
                        textcoords="offset points", xytext=(5, 4),
                        fontsize=7, color=HIGHLIGHT_ACCENT,
                        fontweight="semibold", zorder=5,
                    )

        # 15 % right-side padding
        xmin, xmax_auto = ax.get_xlim()
        data_xmax  = sub[pcol].max()
        data_xrange = data_xmax - sub[pcol].min()
        padded_xmax = data_xmax + 0.15 * data_xrange
        ax.set_xlim(xmin, max(xmax_auto, padded_xmax))

        _style_ax(ax)
        ax.set_xlabel(nice, fontsize=10)
        ax.set_ylabel(flabel, fontsize=10)
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
        description="Plot UTR-specific PhyloP distributions and feature scatter plots.")
    p.add_argument("csv", help="Input CSV file")
    p.add_argument("-c", "--gene_id_column", default=None,
                   help="Column name containing gene IDs (e.g. display_name)")
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

    utr_phylop = present_utr_phylop(df)
    if not utr_phylop:
        sys.exit(
            "ERROR: CSV contains none of the expected UTR PhyloP columns "
            "(utr5_max_phylop, utr5_med_phylop, utr3_max_phylop, utr3_med_phylop)."
        )

    features = present_features(df)
    gene_col  = args.gene_id_column.lower() if args.gene_id_column else None
    if gene_col and gene_col not in df.columns:
        print(f"WARNING: column '{gene_col}' not found — highlighting disabled.",
              file=sys.stderr)
        gene_col = None

    mask = highlight_mask(df, gene_col, args.genes_of_interest)
    if args.genes_of_interest and mask is not None:
        found = mask.sum()
        print(f"Highlighting {found} / {len(args.genes_of_interest)} requested genes.")

    print(f"Found UTR PhyloP columns: {[c for c,_,_ in utr_phylop]}")
    print(f"Found feature columns:    {list(features.keys())}")

    for pcol, nice, utr_tag in utr_phylop:
        plot_distribution(df, pcol, nice, utr_tag, mask, gene_col, args.outdir)
        if features:
            plot_scatter_grid(df, pcol, nice, utr_tag,
                              features, mask, gene_col, args.outdir)
        else:
            print("No recognised feature columns found — skipping scatter plots.")


if __name__ == "__main__":
    main()