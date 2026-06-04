#!/usr/bin/env python3
"""
UTR Alignment Distribution Visualizer
--------------------------------------
Plots distribution statistics from bovine/human/mouse UTR alignment CSV data.

Usage:
    python plot_utr_alignments.py <csv_file> [--output <output_dir>] [--dpi <dpi>]

Example:
    python plot_utr_alignments.py alignments.csv --output ./plots --dpi 150
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns

# ── Palette ───────────────────────────────────────────────────────────────────
PALETTE = {
    "bg":      "#ffffff",
    "panel":   "#f7f8fc",
    "border":  "#dde1ec",
    "grid":    "#eceef5",
    "accent1": "#2563eb",   # 3′ UTR – cobalt blue
    "accent2": "#e11d48",   # 5′ UTR – crimson
    "accent3": "#059669",   # shared  – emerald
    "text":    "#111827",
    "subtext": "#6b7280",
    "fill1":   "#dbeafe",   # light blue fill
    "fill2":   "#ffe4e6",   # light rose fill
}

plt.rcParams.update({
    "figure.facecolor":       PALETTE["bg"],
    "axes.facecolor":         PALETTE["panel"],
    "axes.edgecolor":         PALETTE["border"],
    "axes.labelcolor":        PALETTE["text"],
    "axes.spines.top":        False,
    "axes.spines.right":      False,
    "axes.linewidth":         0.9,
    "xtick.color":            PALETTE["subtext"],
    "ytick.color":            PALETTE["subtext"],
    "xtick.labelsize":        8,
    "ytick.labelsize":        8,
    "text.color":             PALETTE["text"],
    "grid.color":             PALETTE["grid"],
    "grid.linewidth":         0.7,
    "axes.grid":              True,
    "axes.grid.axis":         "y",
    "font.family":            "sans-serif",
    "font.sans-serif":        ["DejaVu Sans"],
    "figure.dpi":             100,
    "savefig.facecolor":      PALETTE["bg"],
})

NUMERIC_COLS = {
    "3utr_pid_mean":        ("3′ UTR % Identity",    PALETTE["accent1"], PALETTE["fill1"]),
    "3utr_coverage_mean":   ("3′ UTR Coverage",      PALETTE["accent1"], PALETTE["fill1"]),
    "3utr_pctdiff_mean":    ("3′ UTR % Difference",  PALETTE["accent1"], PALETTE["fill1"]),
    "5utr_pid_mean":        ("5′ UTR % Identity",    PALETTE["accent2"], PALETTE["fill2"]),
    "5utr_coverage_mean":   ("5′ UTR Coverage",      PALETTE["accent2"], PALETTE["fill2"]),
    "5utr_pctdiff_mean":    ("5′ UTR % Difference",  PALETTE["accent2"], PALETTE["fill2"]),
}

# ── Helpers ───────────────────────────────────────────────────────────────────

def load_data(path: str) -> pd.DataFrame:
    p = Path(path)
    if not p.exists():
        sys.exit(f"[ERROR] File not found: {path}")
    try:
        df = pd.read_csv(p, sep=None, engine="python")
    except Exception as exc:
        sys.exit(f"[ERROR] Could not parse CSV: {exc}")
    df.columns = df.columns.str.strip().str.lower()
    missing = [c for c in NUMERIC_COLS if c not in df.columns]
    if missing:
        print(f"[WARN] Missing columns (will be skipped): {missing}")
    return df


def summary_stats(df: pd.DataFrame, col: str) -> str:
    s = df[col].dropna()
    return (f"n={len(s):,}   median={s.median():.2f}   "
            f"mean={s.mean():.2f}   σ={s.std():.2f}")


def style_ax(ax):
    """Apply clean, consistent axis chrome."""
    ax.set_facecolor(PALETTE["panel"])
    ax.tick_params(length=3, width=0.7)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(PALETTE["border"])
        ax.spines[spine].set_linewidth(0.9)


# ── Plot 1 : individual distributions ────────────────────────────────────────

def plot_distributions(df: pd.DataFrame, output_dir: Path, dpi: int):
    cols_present = [c for c in NUMERIC_COLS if c in df.columns]
    n = len(cols_present)
    if n == 0:
        print("[WARN] No numeric columns found – skipping distribution plot.")
        return

    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(6.5 * ncols, 4.2 * nrows),
                             facecolor=PALETTE["bg"])
    axes_flat = axes.flatten() if n > 1 else [axes]

    fig.suptitle("UTR Alignment — Individual Distributions",
                 fontsize=15, color=PALETTE["text"],
                 fontweight="bold", y=1.02)

    for ax, col in zip(axes_flat, cols_present):
        label, color, fill = NUMERIC_COLS[col]
        data = df[col].dropna()
        style_ax(ax)

        # filled histogram
        ax.hist(data, bins=40, color=fill, edgecolor=color,
                linewidth=0.4, density=True, zorder=2)

        # KDE line
        try:
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(data, bw_method="scott")
            xs = np.linspace(data.min(), data.max(), 400)
            ax.plot(xs, kde(xs), color=color, linewidth=2.2, zorder=3)
        except ImportError:
            pass

        # median dashed line
        med = data.median()
        ax.axvline(med, color=color, linewidth=1.2,
                   linestyle="--", alpha=0.8, zorder=4)
        ax.text(med + (data.max() - data.min()) * 0.01,
                ax.get_ylim()[1] * 0.93,
                f"median\n{med:.2f}",
                color=color, fontsize=7, va="top", fontweight="bold")

        ax.set_title(label, color=PALETTE["text"],
                     fontsize=10, fontweight="semibold", pad=7)
        ax.set_xlabel("Value", fontsize=8, color=PALETTE["subtext"])
        ax.set_ylabel("Density", fontsize=8, color=PALETTE["subtext"])

        ax.text(0.99, 0.99, summary_stats(df, col),
                transform=ax.transAxes, ha="right", va="top",
                fontsize=6.2, color=PALETTE["subtext"],
                bbox=dict(boxstyle="round,pad=0.3",
                          fc=PALETTE["bg"], ec=PALETTE["border"],
                          alpha=0.85, linewidth=0.6))

    for ax in axes_flat[n:]:
        ax.set_visible(False)

    fig.tight_layout(pad=1.5)
    out = output_dir / "01_distributions.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    plt.close(fig)


# ── Plot 2 : 3′ vs 5′ paired comparison ──────────────────────────────────────

def plot_paired_comparison(df: pd.DataFrame, output_dir: Path, dpi: int):
    pairs = [
        ("% Identity",   "3utr_pid_mean",      "5utr_pid_mean"),
        ("Coverage",     "3utr_coverage_mean",  "5utr_coverage_mean"),
        ("% Difference", "3utr_pctdiff_mean",   "5utr_pctdiff_mean"),
    ]
    pairs = [(lbl, c3, c5) for lbl, c3, c5 in pairs
             if c3 in df.columns and c5 in df.columns]
    if not pairs:
        print("[WARN] Not enough columns for paired comparison – skipping.")
        return

    fig, axes = plt.subplots(1, len(pairs),
                             figsize=(5.5 * len(pairs), 5.8),
                             facecolor=PALETTE["bg"])
    if len(pairs) == 1:
        axes = [axes]

    fig.suptitle("3′ UTR vs 5′ UTR — Paired Comparisons",
                 fontsize=15, color=PALETTE["text"],
                 fontweight="bold", y=1.02)

    region_palette = {"3′ UTR": PALETTE["accent1"],
                      "5′ UTR": PALETTE["accent2"]}

    for ax, (label, c3, c5) in zip(axes, pairs):
        style_ax(ax)
        ax.grid(True, axis="y", color=PALETTE["grid"],
                linewidth=0.7, zorder=0)

        tidy = pd.DataFrame({
            "value":  pd.concat([df[c3].dropna(), df[c5].dropna()],
                                ignore_index=True),
            "region": (["3′ UTR"] * df[c3].dropna().shape[0] +
                       ["5′ UTR"] * df[c5].dropna().shape[0]),
        })

        sns.boxplot(
            data=tidy, x="region", y="value",
            palette={"3′ UTR": PALETTE["fill1"],
                     "5′ UTR": PALETTE["fill2"]},
            width=0.48, linewidth=1.1,
            medianprops=dict(linewidth=2),
            flierprops=dict(marker="o", markersize=2,
                            markerfacecolor=PALETTE["subtext"],
                            alpha=0.4, linewidth=0),
            ax=ax, zorder=3,
        )
        # recolour box edges by region
        for patch, key in zip(ax.patches,
                              ["3′ UTR", "5′ UTR"] * 10):
            patch.set_edgecolor(region_palette.get(key, PALETTE["border"]))
            patch.set_linewidth(1.3)

        sns.stripplot(
            data=tidy.sample(min(len(tidy), 500), random_state=42),
            x="region", y="value",
            palette=region_palette,
            alpha=0.22, size=3, jitter=0.18, ax=ax, zorder=2,
        )

        ax.set_title(label, fontsize=11, color=PALETTE["text"],
                     fontweight="semibold", pad=8)
        ax.set_xlabel("")
        ax.set_ylabel("Value", fontsize=9, color=PALETTE["subtext"])

        # annotate n
        for i, (col, key) in enumerate([(c3, "3′ UTR"), (c5, "5′ UTR")]):
            n = df[col].dropna().shape[0]
            ax.text(i, ax.get_ylim()[0] - (ax.get_ylim()[1] -
                    ax.get_ylim()[0]) * 0.06,
                    f"n={n:,}", ha="center", va="top",
                    fontsize=7.5, color=PALETTE["subtext"])

    fig.tight_layout(pad=1.5)
    out = output_dir / "02_paired_comparison.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    plt.close(fig)


# ── Plot 3 : correlation heatmap ──────────────────────────────────────────────

def plot_correlation(df: pd.DataFrame, output_dir: Path, dpi: int):
    cols_present = [c for c in NUMERIC_COLS if c in df.columns]
    if len(cols_present) < 2:
        return

    corr = df[cols_present].corr()
    labels = [NUMERIC_COLS[c][0] for c in cols_present]

    fig, ax = plt.subplots(figsize=(9, 7), facecolor=PALETTE["bg"])
    ax.set_facecolor(PALETTE["bg"])
    fig.suptitle("Pearson Correlation — Alignment Metrics",
                 fontsize=14, color=PALETTE["text"], fontweight="bold")

    cmap = sns.diverging_palette(225, 15, s=75, l=55, as_cmap=True)
    mask = np.zeros_like(corr, dtype=bool)
    mask[np.triu_indices_from(mask)] = True

    sns.heatmap(
        corr, mask=mask, cmap=cmap, vmin=-1, vmax=1,
        annot=True, fmt=".2f",
        annot_kws={"size": 9.5, "color": PALETTE["text"]},
        linewidths=1.2, linecolor=PALETTE["bg"],
        xticklabels=labels, yticklabels=labels,
        ax=ax,
        cbar_kws={"shrink": 0.72, "label": "Pearson r"},
        square=True,
    )

    ax.tick_params(axis="x", rotation=30, colors=PALETTE["text"], labelsize=8.5)
    ax.tick_params(axis="y", rotation=0,  colors=PALETTE["text"], labelsize=8.5)
    ax.spines[:].set_visible(False)

    fig.tight_layout(pad=1.5)
    out = output_dir / "03_correlation_heatmap.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    plt.close(fig)


# ── Plot 4 : identity vs coverage scatter ────────────────────────────────────

def plot_identity_vs_coverage(df: pd.DataFrame, output_dir: Path, dpi: int):
    pairs = [
        ("3′ UTR", "3utr_pid_mean", "3utr_coverage_mean",
         PALETTE["accent1"], PALETTE["fill1"]),
        ("5′ UTR", "5utr_pid_mean", "5utr_coverage_mean",
         PALETTE["accent2"], PALETTE["fill2"]),
    ]
    pairs = [(lbl, x, y, c, f) for lbl, x, y, c, f in pairs
             if x in df.columns and y in df.columns]
    if not pairs:
        return

    fig, axes = plt.subplots(1, len(pairs),
                             figsize=(6 * len(pairs), 5.5),
                             facecolor=PALETTE["bg"])
    if len(pairs) == 1:
        axes = [axes]

    fig.suptitle("% Identity vs Coverage",
                 fontsize=14, color=PALETTE["text"], fontweight="bold")

    for ax, (label, xcol, ycol, color, fill) in zip(axes, pairs):
        style_ax(ax)
        ax.grid(True, color=PALETTE["grid"], linewidth=0.7, zorder=0)

        sub = df[[xcol, ycol]].dropna()
        ax.scatter(sub[xcol], sub[ycol],
                   c=fill, edgecolors=color,
                   linewidths=0.4, alpha=0.65, s=18, zorder=3)

        if len(sub) > 2:
            m, b = np.polyfit(sub[xcol], sub[ycol], 1)
            xs = np.linspace(sub[xcol].min(), sub[xcol].max(), 200)
            ax.plot(xs, m * xs + b, color=color, linewidth=1.8,
                    linestyle="--", alpha=0.85,
                    label=f"slope = {m:.3f}", zorder=4)
            ax.legend(fontsize=8, framealpha=0.7,
                      edgecolor=PALETTE["border"],
                      facecolor=PALETTE["bg"])

        ax.set_title(label, fontsize=11, color=PALETTE["text"],
                     fontweight="semibold")
        ax.set_xlabel("% Identity", fontsize=9, color=PALETTE["subtext"])
        ax.set_ylabel("Coverage",   fontsize=9, color=PALETTE["subtext"])

        r = sub.corr().iloc[0, 1]
        ax.text(0.04, 0.97, f"r = {r:.3f}",
                transform=ax.transAxes, fontsize=9.5,
                color=color, va="top", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3",
                          fc=fill, ec=color,
                          alpha=0.85, linewidth=0.8))

    fig.tight_layout(pad=1.5)
    out = output_dir / "04_identity_vs_coverage.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    plt.close(fig)


# ── Plot 5 : violin overview ──────────────────────────────────────────────────

def plot_violin_all(df: pd.DataFrame, output_dir: Path, dpi: int):
    cols_present = [c for c in NUMERIC_COLS if c in df.columns]
    if not cols_present:
        return

    tidy_rows = []
    for col in cols_present:
        label, color, fill = NUMERIC_COLS[col]
        for val in df[col].dropna():
            tidy_rows.append({"metric": label, "value": val,
                               "color": color, "fill": fill})
    tidy = pd.DataFrame(tidy_rows)

    fig, ax = plt.subplots(figsize=(max(10, len(cols_present) * 2.2), 6),
                           facecolor=PALETTE["bg"])
    style_ax(ax)
    ax.grid(True, axis="y", color=PALETTE["grid"], linewidth=0.7, zorder=0)

    fill_palette = {NUMERIC_COLS[c][0]: NUMERIC_COLS[c][2]
                    for c in cols_present}
    line_palette = {NUMERIC_COLS[c][0]: NUMERIC_COLS[c][1]
                    for c in cols_present}

    vp = sns.violinplot(data=tidy, x="metric", y="value",
                        palette=fill_palette, inner="box",
                        linewidth=1.2, cut=0, ax=ax, zorder=3)

    # recolour violin outlines
    for i, col in enumerate(vp.collections):
        metric = cols_present[i // 2] if i // 2 < len(cols_present) else None
        if metric:
            col.set_edgecolor(NUMERIC_COLS[metric][1])

    ax.set_title("All Metrics — Distribution Overview",
                 fontsize=14, color=PALETTE["text"],
                 fontweight="bold", pad=10)
    ax.set_xlabel("")
    ax.set_ylabel("Value", fontsize=9, color=PALETTE["subtext"])
    ax.tick_params(axis="x", rotation=18, colors=PALETTE["text"])
    ax.tick_params(axis="y", colors=PALETTE["subtext"])

    fig.tight_layout(pad=1.5)
    out = output_dir / "05_violin_overview.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"  Saved → {out}")
    plt.close(fig)


# ── CLI ────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualise UTR alignment distribution data from a CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("csv_file",
                        help="Path to the alignment CSV.")
    parser.add_argument("--output", "-o", default="./utr_plots",
                        help="Output directory (default: ./utr_plots).")
    parser.add_argument("--dpi", type=int, default=150,
                        help="Figure resolution (default: 150).")
    parser.add_argument("--no-show", action="store_true",
                        help="Do not open an interactive window.")
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n[INFO] Loading  → {args.csv_file}")
    df = load_data(args.csv_file)
    print(f"[INFO] Rows: {len(df):,}  |  Columns: {list(df.columns)}\n")

    numeric_present = [c for c in NUMERIC_COLS if c in df.columns]
    if numeric_present:
        print("[INFO] Numeric summary:")
        print(df[numeric_present].describe().round(3).to_string())
        print()

    print("[INFO] Generating plots …")
    plot_distributions(df, output_dir, args.dpi)
    plot_paired_comparison(df, output_dir, args.dpi)
    plot_correlation(df, output_dir, args.dpi)
    plot_identity_vs_coverage(df, output_dir, args.dpi)
    plot_violin_all(df, output_dir, args.dpi)

    print(f"\n[DONE] All plots saved to → {output_dir.resolve()}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()