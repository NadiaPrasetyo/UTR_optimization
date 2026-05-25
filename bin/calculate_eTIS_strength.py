import pandas as pd
import numpy as np
import sys
import os
import argparse
import matplotlib.pyplot as plt
from scipy import stats
PALETTE = {
    'half_life': ('#2563EB', '#BFDBFE'),   # blue  (accent, light)
    'mean_te':   ('#DC2626', '#FEE2E2'),   # red
}

def plot_etis_strength(df, output_file):
    def make_scatter(ax, x, y, xlabel, ylabel, title, accent, light):
        ax.scatter(
            x, y,
            marker='o', s=18, alpha=0.35,
            color=light, edgecolors=accent,
            linewidths=0.4, rasterized=True,
        )

        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() > 2:
            slope, intercept, r_val, p_val, _ = stats.linregress(x[mask], y[mask])
            xs = np.linspace(x[mask].min(), x[mask].max(), 300)
            ax.plot(xs, slope * xs + intercept, color=accent, linewidth=2)

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

    base = os.path.splitext(output_file)[0]

    for y_col, ylabel, suffix in [
        ('half_life', 'Half-Life (z-score)',               '_hl'),
        ('mean_te',   'Translation Efficiency (z-score)',  '_te'),
    ]:
        accent, light = PALETTE[y_col]
        fig, ax = plt.subplots(figsize=(6, 5))
        make_scatter(
            ax,
            x=df['eTIS_strength'].values,
            y=df[y_col].values,
            xlabel='eTIS Strength',
            ylabel=ylabel,
            title=f'eTIS Strength vs. {ylabel}',
            accent=accent,
            light=light,
        )
        fig.tight_layout()
        out_path = f"{base}{suffix}.png"
        fig.savefig(out_path, dpi=180, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved plot to {out_path}")
    

def main(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    #eTIS strength = 100 − (predicted leaky scanning / maximum predicted leaky scanning) × 100
    max_leaky_scanning = df['predictions_GFP'].max()
    df['eTIS_strength'] = 100 - (df['predictions_GFP'] / max_leaky_scanning) * 100
    df.to_csv(output_file, sep='\t', index=False)

    print(f"Done! Calculated the eTIS strength on {len(df)} rows based on maximum predicted leaky scanning: {max_leaky_scanning}")

    plot_etis_strength(df, output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate eTIS strength')
    parser.add_argument('-i', '--input', required=True, type=str, help='Input TSV file')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output TSV file')
    args = parser.parse_args()

    if not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    main(args.input, args.output)