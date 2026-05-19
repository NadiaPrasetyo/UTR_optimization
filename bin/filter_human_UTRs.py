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
    Run KS tests for each UTR feature against half_life and mean_te.
    Splits data into high/low groups by median of the target variable and compares
    the feature distributions between groups.
    """
    features = ['utr5_length', 'utr3_length', 'cds_length',
                'utr5_gc',    'utr3_gc',    'cds_gc']
    targets  = ['half_life', 'mean_te']

    feature_labels = {
        'utr5_length': "5' UTR Length",
        'utr3_length': "3' UTR Length",
        'cds_length':  "CDS Length",
        'utr5_gc':     "5' UTR GC Content",
        'utr3_gc':     "3' UTR GC Content",
        'cds_gc':      "CDS GC Content",
    }
    target_labels = {
        'half_life': 'Half-Life',
        'mean_te':   'Translation Efficiency',
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


def plot_scatter_with_trendline(df, output_dir, feature_group="hl_stats"):
    """
    Scatter plots with linear regression trendlines.

    In both modes:
      X-axis = sequence feature (length or GC content)  ← predictor
      Y-axis = biological target (half_life or mean_te) ← response

    hl_stats  — one figure per sequence feature (6 figures), two panels each
                (half_life | mean_te).
    utr_stats — one figure per biological target (2 figures), six panels each
                (utr5_length | utr3_length | cds_length |
                 utr5_gc     | utr3_gc     | cds_gc).

    Output files: scatter_<feature>.png  or  scatter_<target>.png
    """

    # ── shared metadata ────────────────────────────────────────────────
    bio_targets  = ['half_life', 'mean_te']

    bio_labels = {
        'half_life': 'Half-Life (z-score)',
        'mean_te':   'Translation Efficiency (z-score)',
    }

    # One colour per biological target
    palette = {
        'half_life': ('#2563EB', '#BFDBFE'),   # blue
        'mean_te':   ('#DC2626', '#FEE2E2'),   # red
    }

    def _draw_panel(ax, x_vals, y_vals, xlabel, ylabel, title, accent, light):
        """Draw scatter + regression on a single axes."""
        ax.scatter(x_vals, y_vals, alpha=0.35, s=18,
                   color=light, edgecolors=accent,
                   linewidths=0.4, rasterized=True)

        slope, intercept, r_val, p_val, _ = stats.linregress(x_vals, y_vals)
        x_line = np.linspace(x_vals.min(), x_vals.max(), 300)
        ax.plot(x_line, slope * x_line + intercept, color=accent, linewidth=2)

        r2  = r_val ** 2
        sig = ('***' if p_val < 0.001 else
               '**'  if p_val < 0.01  else
               '*'   if p_val < 0.05  else 'ns')
        ax.text(0.05, 0.95,
                f"R² = {r2:.3f}\np = {p_val:.2e} {sig}",
                transform=ax.transAxes, fontsize=8.5, verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                          edgecolor=accent, alpha=0.85))

        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11, fontweight='semibold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=8)

    utr_features = ['utr5_length', 'utr5_gc', 'utr3_length', 'utr3_gc']
    utr_panel_titles = {
        'utr5_length': "5' UTR Length (nt)",
        'utr5_gc':     "5' UTR GC Content (%)",
        'utr3_length': "3' UTR Length (nt)",
        'utr3_gc':     "3' UTR GC Content (%)",
    }

    # ── hl_stats: one figure per bio target, 4 panels (utr5_length |
    #             utr5_gc | utr3_length | utr3_gc)
    if feature_group == "hl_stats":
        for bio_tgt in bio_targets:
            accent, light = palette[bio_tgt]
            fig, axes = plt.subplots(1, 4, figsize=(22, 5))
            fig.suptitle(bio_labels[bio_tgt],
                         fontsize=14, fontweight='bold', y=1.02)

            for ax, utr_feat in zip(axes, utr_features):
                sub = df[[utr_feat, bio_tgt]].dropna()
                _draw_panel(ax,
                            x_vals=sub[utr_feat].values,
                            y_vals=sub[bio_tgt].values,
                            xlabel=utr_panel_titles[utr_feat],
                            ylabel=bio_labels[bio_tgt],
                            title=utr_panel_titles[utr_feat],
                            accent=accent, light=light)

            plt.tight_layout()
            out_path = os.path.join(output_dir, f'scatter_{bio_tgt}.png')
            fig.savefig(out_path, dpi=180, bbox_inches='tight')
            plt.close(fig)
            print(f"Scatter plot saved: {out_path}")

    # ── utr_stats: one figure per UTR feature, 2 panels (half_life | mean_te)
    elif feature_group == "utr_stats":
        for utr_feat in utr_features:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            fig.suptitle(utr_panel_titles[utr_feat],
                         fontsize=14, fontweight='bold', y=1.02)

            for ax, bio_tgt in zip(axes, bio_targets):
                accent, light = palette[bio_tgt]
                sub = df[[utr_feat, bio_tgt]].dropna()
                _draw_panel(ax,
                            x_vals=sub[utr_feat].values,
                            y_vals=sub[bio_tgt].values,
                            xlabel=utr_panel_titles[utr_feat],
                            ylabel=bio_labels[bio_tgt],
                            title=bio_labels[bio_tgt],
                            accent=accent, light=light)

            plt.tight_layout()
            out_path = os.path.join(output_dir, f'scatter_{utr_feat}.png')
            fig.savefig(out_path, dpi=180, bbox_inches='tight')
            plt.close(fig)
            print(f"Scatter plot saved: {out_path}")


def plot_gc_content_outliers(df, output_dir, outlier_zscore=2.0):
    """
    Produces two plots and one CSV from UTR GC-content data.

    4-D scatter (3D + colour):
      X     = 5' UTR GC content (%)
      Y     = 3' UTR GC content (%)
      Z     = 5' UTR length (nt)
      Colour= 3' UTR length (nt)  via plasma colormap

    2-D scatter:
      X     = 5' UTR GC content (%)
      Y     = 3' UTR GC content (%)

    Outliers are defined by GC-content z-scores only.
    Baseline references (Human HBB-201, Bovine HBB-204) are plotted
    at their actual UTR values using the same axes.

    Required columns in df:
        utr5_gc, utr3_gc, utr5_length, utr3_length, accession

    Outputs:
        scatter_utr_gc_content_4d.png
        scatter_utr_gc_content_2d.png
        gc_outliers.csv
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 — registers projection
    from scipy import stats
    from matplotlib.lines import Line2D

    # ------------------------------------------------------------------ #
    #  Helper                                                              #
    # ------------------------------------------------------------------ #
    def _gc(seq):
        seq = seq.replace('\n', '').replace(' ', '').upper()
        if not seq:
            return np.nan
        return (seq.count('G') + seq.count('C')) / len(seq) * 100

    # ------------------------------------------------------------------ #
    #  Baseline reference sequences                                        #
    # ------------------------------------------------------------------ #
    baselines = {
        'Human HBB-201': {
            'utr5': 'ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACC',
            'utr3': (
                'GCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACT'
                'AAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTA'
                'TTTTCATTGCAA'
            ),
        },
        'Bovine HBB-204': {
            'utr5': 'CTTACACTTGCTTCTGACACAACCGTGTTCACTAGCAACTACACAAACAGACACC',
            'utr3': (
                'GCTCCCCTTCCTGATTTTCAGGAAAGGTCTTTTCATCCTCAGAGCCCAAAAACTGAATATG'
                'GAAAAATTATGAAGCGTTTTGTGCATCTTGCCTCTGCCTAATAAAGACATTTATTTTCATT'
                'GCACTGGTGTATTT'
            ),
        },
    }

    baseline_styles = {
        'Human HBB-201': dict(
            facecolor='#16A34A', edgecolor='#14532D',
            marker='*', s=400, zorder=10,
            label='Human HBB-201 (baseline)',
        ),
        'Bovine HBB-204': dict(
            facecolor='#F59E0B', edgecolor='#78350F',
            marker='*', s=400, zorder=10,
            label='Bovine HBB-204 (baseline)',
        ),
    }

    for name, seqs in baselines.items():
        raw5 = seqs['utr5'].replace('\n', '').replace(' ', '')
        raw3 = seqs['utr3'].replace('\n', '').replace(' ', '')
        baselines[name]['gc5']  = _gc(raw5)
        baselines[name]['gc3']  = _gc(raw3)
        baselines[name]['len5'] = len(raw5)
        baselines[name]['len3'] = len(raw3)

    # ------------------------------------------------------------------ #
    #  Main dataset                                                        #
    # ------------------------------------------------------------------ #
    needed = ['utr5_gc', 'utr3_gc', 'utr5_length', 'utr3_length', 'accession']
    sub    = df.dropna(subset=needed).reset_index(drop=True)

    x      = sub['utr5_gc'].values
    y      = sub['utr3_gc'].values
    z      = sub['utr5_length'].values
    len3   = sub['utr3_length'].values
    labels = sub['accession'].values

    # Outlier mask — GC-content z-scores only
    zx         = stats.zscore(x)
    zy         = stats.zscore(y)
    is_outlier = (np.abs(zx) > outlier_zscore) | (np.abs(zy) > outlier_zscore)

    # ------------------------------------------------------------------ #
    #  Shared colour normalisation (dataset + baselines)                  #
    # ------------------------------------------------------------------ #
    all_len3 = np.concatenate([len3, [b['len3'] for b in baselines.values()]])
    norm     = mcolors.Normalize(vmin=all_len3.min(), vmax=all_len3.max())
    cmap     = cm.plasma

    # ================================================================== #
    #  Plot 1 — 4-D scatter (3-D axes + colour)                          #
    # ================================================================== #
    fig = plt.figure(figsize=(11, 8))
    ax  = fig.add_subplot(111, projection='3d')

    ax.scatter(
        x[~is_outlier], y[~is_outlier], z[~is_outlier],
        c=len3[~is_outlier], cmap=cmap, norm=norm,
        s=22, alpha=0.50, edgecolors='none',
        depthshade=True, rasterized=True,
    )
    ax.scatter(
        x[is_outlier], y[is_outlier], z[is_outlier],
        c=len3[is_outlier], cmap=cmap, norm=norm,
        s=80, alpha=0.92, edgecolors='#DC2626', linewidths=0.8,
        depthshade=True, zorder=5,
    )

    x_off = (x.max() - x.min()) * 0.02
    y_off = (y.max() - y.min()) * 0.01
    z_off = (z.max() - z.min()) * 0.01
    for xi, yi, zi, lbl in zip(x[is_outlier], y[is_outlier],
                                z[is_outlier], labels[is_outlier]):
        ax.text(xi + x_off, yi + y_off, zi + z_off, lbl,
                fontsize=6.5, color='#7F1D1D', va='bottom', ha='left')

    for name, seqs in baselines.items():
        style = baseline_styles[name]
        ax.scatter(
            seqs['gc5'], seqs['gc3'], seqs['len5'],
            c=[cmap(norm(seqs['len3']))],
            marker='*', s=400,
            edgecolors=style['edgecolor'], linewidths=0.9,
            depthshade=False, zorder=10,
        )
        ax.text(seqs['gc5'], seqs['gc3'], seqs['len5'],
                f"      {name}",
                fontsize=7, color=style['edgecolor'], fontweight='bold')

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.022, pad=0.12, shrink=0.55)
    cbar.set_label("3' UTR length (nt)", fontsize=9)

    from matplotlib.lines import Line2D
    cat_handles = [
        Line2D([0], [0], linestyle='none', marker='o', markersize=7,
               markerfacecolor='#A78BFA', markeredgewidth=0, label='Normal'),
        Line2D([0], [0], linestyle='none', marker='o', markersize=9,
               markerfacecolor='#F87171', markeredgecolor='#DC2626',
               markeredgewidth=0.6,
               label=f'Outlier (|z| > {outlier_zscore}, GC only)'),
        Line2D([0], [0], linestyle='none', marker='*', markersize=13,
               markerfacecolor='#16A34A', markeredgecolor='#14532D',
               label='Human HBB-201 (baseline)'),
        Line2D([0], [0], linestyle='none', marker='*', markersize=13,
               markerfacecolor='#F59E0B', markeredgecolor='#78350F',
               label='Bovine HBB-204 (baseline)'),
    ]
    ax.legend(handles=cat_handles, fontsize=8, loc='upper left', framealpha=0.85)

    ax.set_xlabel("5' UTR GC content (%)", fontsize=10, labelpad=8)
    ax.set_ylabel("3' UTR GC content (%)", fontsize=10, labelpad=8)
    ax.set_zlabel("5' UTR length (nt)",    fontsize=10, labelpad=8)
    ax.set_title(
        "UTR GC content & length  |  colour = 3' UTR length (nt)",
        fontsize=11, fontweight='bold', pad=14,
    )
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False
        pane.set_edgecolor('#CCCCCC')
    ax.grid(True, linestyle='--', linewidth=0.4, alpha=0.5)
    fig.text(0.01, 0.98, f"n = {len(sub)}  |  outliers = {is_outlier.sum()}",
             fontsize=8, va='top', color='grey')

    plt.tight_layout()
    out_4d = os.path.join(output_dir, 'scatter_utr_gc_content_4d.png')
    fig.savefig(out_4d, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f"4-D GC content plot saved: {out_4d}")

    # ================================================================== #
    #  Plot 2 — 2-D scatter                                               #
    # ================================================================== #
    fig, ax = plt.subplots(figsize=(9, 7))

    ax.scatter(x[~is_outlier], y[~is_outlier],
               alpha=0.4, s=20, color='#60A5FA', edgecolors='#1D4ED8',
               linewidths=0.3, rasterized=True, label='Normal')
    ax.scatter(x[is_outlier], y[is_outlier],
               alpha=0.85, s=55, color='#FCA5A5', edgecolors='#DC2626',
               linewidths=0.8, zorder=5,
               label=f'Outlier (|z| > {outlier_zscore} in either)')

    texts = []
    for xi, yi, lbl in zip(x[is_outlier], y[is_outlier], labels[is_outlier]):
        texts.append(ax.text(xi, yi, f'  {lbl}',
                             fontsize=7, color='#7F1D1D', va='center'))

    for name, seqs in baselines.items():
        style = baseline_styles[name]
        ax.scatter(
            seqs['gc5'], seqs['gc3'],
            marker=style['marker'], s=style['s'], zorder=style['zorder'],
            color=style['facecolor'], edgecolors=style['edgecolor'],
            linewidths=0.9, label=style['label'],
        )
        texts.append(ax.text(
            seqs['gc5'], seqs['gc3'],
            f"  {name}\n  5′={seqs['gc5']:.1f}%  3′={seqs['gc3']:.1f}%",
            fontsize=7.5, color=style['edgecolor'],
            va='center', fontweight='bold',
        ))

    try:
        from adjustText import adjust_text
        adjust_text(
            texts,
            x=np.concatenate([x[is_outlier], [b['gc5'] for b in baselines.values()]]),
            y=np.concatenate([y[is_outlier], [b['gc3'] for b in baselines.values()]]),
            ax=ax,
            arrowprops=dict(arrowstyle='-', color='grey', lw=0.5),
        )
    except ImportError:
        pass

    ax.set_xlabel("5' UTR GC Content (%)", fontsize=12)
    ax.set_ylabel("3' UTR GC Content (%)", fontsize=12)
    ax.set_title("5' UTR vs 3' UTR GC Content", fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, framealpha=0.85)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.text(0.02, 0.98, f"n = {len(sub)}  |  outliers = {is_outlier.sum()}",
            transform=ax.transAxes, fontsize=8, va='top', color='grey')

    plt.tight_layout()
    out_2d = os.path.join(output_dir, 'scatter_utr_gc_content_2d.png')
    fig.savefig(out_2d, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f"2-D GC content plot saved: {out_2d}")

    # ================================================================== #
    #  Outlier CSV                                                         #
    # ================================================================== #
    outlier_df = sub[is_outlier].copy()
    outlier_df.insert(outlier_df.columns.get_loc('utr5_gc') + 1,
                      'zscore_utr5_gc', zx[is_outlier].round(3))
    outlier_df.insert(outlier_df.columns.get_loc('utr3_gc') + 1,
                      'zscore_utr3_gc', zy[is_outlier].round(3))
    outlier_csv = os.path.join(output_dir, 'gc_outliers.csv')
    outlier_df.to_csv(outlier_csv, index=False)
    print(f"GC outlier list saved: {outlier_csv}")

    print("\nBaseline reference:")
    for name, seqs in baselines.items():
        print(
            f"  {name}:  5' GC={seqs['gc5']:.1f}%  len={seqs['len5']}nt"
            f"   |   3' GC={seqs['gc3']:.1f}%  len={seqs['len3']}nt"
        )


def main(input_file, output_file, sort, utr5_min, utr5_max, utr3_min, utr3_max, feature_group):
    """
    Load, filter, sort the input CSV, then run KS tests and scatter plots.

    Expected input columns:
        accession, ensembl_gene_id, ensembl_transcript_id,
        utr3, utr5, cds,
        utr5_length, utr5_gc, utr3_length, utr3_gc,
        cds_length, cds_gc,
        mean_te, half_life
    """
    df = pd.read_csv(input_file)

    # Drop rows with empty UTR sequences
    df = df[(df['utr5'].notna() & (df['utr5'] != '')) &
            (df['utr3'].notna() & (df['utr3'] != ''))]

    # Apply UTR length filters (from CLI args)
    df = df[df['utr5_length'].between(utr5_min, utr5_max)]
    df = df[df['utr3_length'].between(utr3_min, utr3_max)]

    # Keep rows where at least one of hl/te is non-negative
    df = df[(df['half_life'] >= 0) | (df['mean_te'] >= 0)]

    # Sort
    sort_cases = {
        'length':                ['utr5_length', 'utr3_length'],
        'half-life':             ['half_life'],
        'transcript-efficiency': ['mean_te'],
    }
    for sort_key in (sort if isinstance(sort, list) else [sort]):
        for col in sort_cases.get(sort_key, []):
            df = df.sort_values(by=col, ascending=False)

    df.to_csv(output_file, index=False)
    print(f"Filtered data written to: {output_file}  ({len(df)} rows)")

    output_dir = os.path.dirname(output_file) or '.'

    run_ks_tests(df, output_dir)
    plot_scatter_with_trendline(df, output_dir, feature_group)
    plot_gc_content_outliers(df, output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Filter a UTR feature CSV, run KS tests, and plot scatter figures.')
    parser.add_argument('-i', '--input_file', required=True, type=str,
                        help='Input CSV file')
    parser.add_argument('-o', '--output_file', default='data/seq_length_gc.csv',
                        type=str, help='Output CSV file (default: data/seq_length_gc.csv)')

    # ── 5' UTR length filter ───────────────────────────────────────────
    parser.add_argument('--utr5-min', default=10, type=int, metavar='NT',
                        help="Minimum 5' UTR length in nt (default: 10)")
    parser.add_argument('--utr5-max', default=100, type=int, metavar='NT',
                        help="Maximum 5' UTR length in nt (default: 100)")

    # ── 3' UTR length filter ───────────────────────────────────────────
    parser.add_argument('--utr3-min', default=10, type=int, metavar='NT',
                        help="Minimum 3' UTR length in nt (default: 10)")
    parser.add_argument('--utr3-max', default=600, type=int, metavar='NT',
                        help="Maximum 3' UTR length in nt (default: 600)")

    parser.add_argument('--sort',
                        choices=['length', 'half-life', 'transcript-efficiency'],
                        nargs='+', default=['length'], type=str,
                        help='Column(s) to sort by (default: length)')
    parser.add_argument('--group-plots',
                        choices=['hl_stats', 'utr_stats'],
                        default='hl_stats', type=str,
                        help='Feature group to plot (default: hl_stats)')

    args = parser.parse_args()

    # Validate filter ranges
    if args.utr5_min >= args.utr5_max:
        parser.error(f"--utr5-min ({args.utr5_min}) must be less than "
                     f"--utr5-max ({args.utr5_max})")
    if args.utr3_min >= args.utr3_max:
        parser.error(f"--utr3-min ({args.utr3_min}) must be less than "
                     f"--utr3-max ({args.utr3_max})")

    # Ensure output directory exists
    out_dir = os.path.dirname(args.output_file)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    main(
        args.input_file,
        args.output_file,
        args.sort,
        args.utr5_min,
        args.utr5_max,
        args.utr3_min,
        args.utr3_max,
        args.group_plots,
    )