import argparse
import csv
import pandas as pd
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from scipy.stats import ks_2samp


def run_ks_tests(df, output_dir):
    """
    Run KS tests for each feature against hl_human, hl_mouse, te_human, te_mouse.
    Splits data into high/low groups by median of the target variable and compares
    the feature distributions between groups.
    """
    features = ['utr3_length', 'utr5_length', 'utr3_gc', 'utr5_gc']
    targets  = ['hl_human', 'hl_mouse', 'te_human', 'te_mouse']

    feature_labels = {
        'utr3_length': "3' UTR Length",
        'utr5_length': "5' UTR Length",
        'utr3_gc':     "3' UTR GC Content",
        'utr5_gc':     "5' UTR GC Content",
    }
    target_labels = {
        'hl_human': 'Half-Life (Human)',
        'hl_mouse': 'Half-Life (Mouse)',
        'te_human': 'Translation Efficiency (Human)',
        'te_mouse': 'Translation Efficiency (Mouse)',
    }

    results = []

    print("\n" + "="*70)
    print("KS TEST RESULTS")
    print("="*70)
    print(f"{'Feature':<22} {'Target':<30} {'KS Stat':>8} {'p-value':>12} {'Sig':>5}")
    print("-"*70)

    for target in targets:
        col = df[target].dropna()
        median_val = col.median()

        for feature in features:
            valid = df[[feature, target]].dropna()
            high_group = valid.loc[valid[target] >= median_val, feature]
            low_group  = valid.loc[valid[target] <  median_val, feature]

            ks_stat, p_val = ks_2samp(high_group, low_group)
            sig = ('***' if p_val < 0.001 else
                   '**'  if p_val < 0.01  else
                   '*'   if p_val < 0.05  else
                   'ns')

            results.append({
                'Feature':      feature_labels[feature],
                'Target':       target_labels[target],
                'KS Statistic': round(ks_stat, 4),
                'p-value':      round(p_val, 6),
                'Significance': sig,
                'n_high':       len(high_group),
                'n_low':        len(low_group),
            })

            print(f"{feature_labels[feature]:<22} {target_labels[target]:<30} "
                  f"{ks_stat:>8.4f} {p_val:>12.2e} {sig:>5}")

    print("="*70)
    print("Significance codes: *** p<0.001  ** p<0.01  * p<0.05  ns not significant")
    print("Groups split by median of target variable (high >= median, low < median)\n")

    ks_df = pd.DataFrame(results)
    ks_path = os.path.join(output_dir, 'ks_test_results.csv')
    ks_df.to_csv(ks_path, index=False)
    print(f"KS test results saved to: {ks_path}")
    return ks_df


def plot_scatter_with_trendline(df, output_dir):
    """
    One figure per feature, with four panels — one per target variable
    (hl_human, hl_mouse, te_human, te_mouse).
    Each panel shows a scatter plot + linear regression trend line,
    annotated with R² and p-value.
    Output files: scatter_utr5_length.png, scatter_utr3_length.png, etc.
    """
    features = ['utr5_length', 'utr3_length', 'utr5_gc', 'utr3_gc']
    targets  = ['hl_human', 'hl_mouse', 'te_human', 'te_mouse']

    feature_labels = {
        'utr3_length': "3' UTR Length (nt)",
        'utr5_length': "5' UTR Length (nt)",
        'utr3_gc':     "3' UTR GC Content",
        'utr5_gc':     "5' UTR GC Content",
    }
    target_labels = {
        'hl_human': 'Half-Life\nHuman (z-score)',
        'hl_mouse': 'Half-Life\nMouse (z-score)',
        'te_human': 'Translation Efficiency\nHuman (z-score)',
        'te_mouse': 'Translation Efficiency\nMouse (z-score)',
    }

    # One colour per target so panels within each figure are visually distinct
    palette = {
        'hl_human': ('#2563EB', '#BFDBFE'),   # blue
        'hl_mouse': ('#16A34A', '#BBF7D0'),   # green
        'te_human': ('#DC2626', '#FEE2E2'),   # red
        'te_mouse': ('#7C3AED', '#EDE9FE'),   # purple
    }

    # One figure per feature, four panels across (one per target)
    for feature in features:
        fig, axes = plt.subplots(1, 4, figsize=(20, 5))
        fig.suptitle(feature_labels[feature], fontsize=14, fontweight='bold', y=1.02)

        for ax, target in zip(axes, targets):
            accent, light = palette[target]

            sub = df[[feature, target]].dropna()
            x = sub[feature].values
            y = sub[target].values

            # Scatter
            ax.scatter(x, y, alpha=0.35, s=18, color=light, edgecolors=accent,
                       linewidths=0.4, rasterized=True)

            # Linear regression trend line
            slope, intercept, r_val, p_val, _ = stats.linregress(x, y)
            x_line = np.linspace(x.min(), x.max(), 300)
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, color=accent, linewidth=2)

            # Annotation: R² and p-value
            r2  = r_val ** 2
            sig = ('***' if p_val < 0.001 else
                   '**'  if p_val < 0.01  else
                   '*'   if p_val < 0.05  else
                   'ns')
            ax.text(0.05, 0.95,
                    f"R² = {r2:.3f}\np = {p_val:.2e} {sig}",
                    transform=ax.transAxes, fontsize=8.5,
                    verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                              edgecolor=accent, alpha=0.85))

            ax.set_xlabel(feature_labels[feature], fontsize=10)
            ax.set_ylabel('z-score', fontsize=10)
            ax.set_title(target_labels[target], fontsize=11, fontweight='semibold')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.tick_params(labelsize=8)

        plt.tight_layout()
        out_path = os.path.join(output_dir, f'scatter_{feature}.png')
        fig.savefig(out_path, dpi=180, bbox_inches='tight')
        plt.close(fig)
        print(f"Scatter plot saved: {out_path}")



def plot_gc_content_outliers(df, output_dir, outlier_zscore=2.0):
    """
    Scatter plot of 5' UTR GC content (x) vs 3' UTR GC content (y).
    Points that are outliers in EITHER axis (|z-score| > outlier_zscore)
    are highlighted and labelled with their UniProt accession.
    Uses adjustText to avoid label overlaps where possible.
    Output: scatter_utr_gc_content.png
    """
    sub = df[['uniprot_accession', 'utr5_gc', 'utr3_gc']].dropna()
    x = sub['utr5_gc'].values
    y = sub['utr3_gc'].values
    labels = sub['uniprot_accession'].values

    # Identify positive outliers: points > outlier_zscore SDs above the mean
    # in EITHER 5' or 3' UTR GC content (high in at least one axis)
    zx = stats.zscore(x)
    zy = stats.zscore(y)
    is_outlier = (zx > outlier_zscore) | (zy > outlier_zscore)

    fig, ax = plt.subplots(figsize=(9, 7))

    # Background points
    ax.scatter(x[~is_outlier], y[~is_outlier],
               alpha=0.4, s=20, color='#60A5FA', edgecolors='#1D4ED8',
               linewidths=0.3, rasterized=True, label='Normal')

    # Outlier points
    ax.scatter(x[is_outlier], y[is_outlier],
               alpha=0.85, s=55, color='#FCA5A5', edgecolors='#DC2626',
               linewidths=0.8, zorder=5, label=f'Outlier (z > {outlier_zscore} in either)')

    # Labels placed directly beside each point, no arrows
    texts = []
    for xi, yi, lbl in zip(x[is_outlier], y[is_outlier], labels[is_outlier]):
        texts.append(ax.text(
            xi, yi, f'  {lbl}',
            fontsize=7, color='#7F1D1D', va='center',
        ))

    try:
        from adjustText import adjust_text
        adjust_text(texts, x=x[is_outlier], y=y[is_outlier], ax=ax)
    except ImportError:
        pass  # adjustText optional; labels still appear without it

    ax.set_xlabel("5' UTR GC Content (%)", fontsize=12)
    ax.set_ylabel("3' UTR GC Content (%)", fontsize=12)
    ax.set_title("5' UTR vs 3' UTR GC Content", fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.85)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    n_outliers = is_outlier.sum()
    ax.text(0.02, 0.98, f"n = {len(sub)}  |  outliers = {n_outliers}",
            transform=ax.transAxes, fontsize=8, va='top', color='grey')

    plt.tight_layout()
    out_path = os.path.join(output_dir, 'scatter_utr_gc_content.png')
    fig.savefig(out_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f"GC content outlier plot saved: {out_path}")

    # Also save the outlier list as CSV
    outlier_df = sub[is_outlier].copy()
    outlier_df['zscore_utr5_gc'] = zx[is_outlier].round(3)
    outlier_df['zscore_utr3_gc'] = zy[is_outlier].round(3)
    outlier_csv = os.path.join(output_dir, 'gc_outliers.csv')
    outlier_df.to_csv(outlier_csv, index=False)
    print(f"GC outlier list saved: {outlier_csv}")

def main(input_file, output_file, sort):
    # open and read the input file
    # uniprot_accession,cds,utr5,utr3,hl_human,hl_mouse,te_human,te_mouse,
    # cds_length,cds_gc,utr5_length,utr5_gc,utr3_length,utr3_gc
    df = pd.read_csv(input_file)

    # remove any data with no utr5 or utr3
    df = df[(df['utr5'] != '') & (df['utr3'] != '')]

    # filter the data
    df = df[df['utr5_length'] <= 100]
    df = df[df['utr5_length'] >= 10]
    df = df[df['utr3_length'] <= 600]
    df = df[df['utr3_length'] >= 10]

    # also filter any lines where the HL and TE are all negative
    df = df[(df['hl_human'] >= 0) | (df['hl_mouse'] >= 0) | (df['te_human'] >= 0) | (df['te_mouse'] >= 0)]

    # sort the data according to the [sort]
    sort_cases = {
        'length':                ['cds_length', 'utr5_length', 'utr3_length'],
        'half-life':             ['hl_human', 'hl_mouse'],
        'transcript-efficiency': ['te_human', 'te_mouse'],
    }
    if "length" in sort:
        for sort_case in sort_cases['length']:
            df = df.sort_values(by=sort_case, ascending=False)
    if "half-life" in sort:
        for sort_case in sort_cases['half-life']:
            df = df.sort_values(by=sort_case, ascending=False)
    if "transcript-efficiency" in sort:
        for sort_case in sort_cases['transcript-efficiency']:
            df = df.sort_values(by=sort_case, ascending=False)

    # write the filtered/sorted data to the output file
    df.to_csv(output_file, index=False)
    print(f"Filtered data written to: {output_file}  ({len(df)} rows)")

    # ── Analysis output directory (same folder as the output CSV) ──────────
    output_dir = os.path.dirname(output_file) or '.'

    # ── KS tests ───────────────────────────────────────────────────────────
    run_ks_tests(df, output_dir)

    # ── Scatter plots ──────────────────────────────────────────────────────
    plot_scatter_with_trendline(df, output_dir)

    # ── GC content outlier plot ────────────────────────────────────────────
    plot_gc_content_outliers(df, output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extract length and GC content from a CSV summary file, '
                    'run KS tests, and plot feature scatter plots.')
    parser.add_argument('-i', '--input_file', required=True, type=str,
                        help='Input CSV file')
    parser.add_argument('-o', '--output_file', default='data/seq_length_gc.csv',
                        type=str, help='Output CSV file (default: data/seq_length_gc.csv)')
    parser.add_argument('--sort',
                        choices=['length', 'half-life', 'transcript-efficiency'],
                        nargs='+', default='length', type=str,
                        help='Column(s) to sort by (default: length)')
    args = parser.parse_args()

    # ensure output directory exists
    out_dir = os.path.dirname(args.output_file)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    main(args.input_file, args.output_file, args.sort)