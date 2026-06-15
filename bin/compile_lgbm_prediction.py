import argparse
import pandas as pd
import numpy as np
import os
import time
import json
import requests
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats

ENSEMBL_REST_BASE = "https://rest.ensembl.org"
ENSEMBL_LOOKUP_BATCH_URL = ENSEMBL_REST_BASE + "/lookup/id"
ENSEMBL_BATCH_SIZE = 100
ENSEMBL_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}
# Rate-limiting / retry settings
ENSEMBL_SLEEP_S = 0.1
ENSEMBL_RETRY_MAX = 3
ENSEMBL_RETRY_SLEEP_S = 5


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-dir', required=True, type=str,
                        help='Input directory that contains the prediction files')
    parser.add_argument('-tsv', '--tsv-file', required=True, type=str,
                        help='Input TSV file with data on transcript IDs and half-lives')
    parser.add_argument('-o', '--output', required=True, type=str,
                        help='Output CSV file')
    return parser.parse_args()  # fixed: was missing ()


def batch_lookup_ensembl(transcript_ids: list[str]) -> dict[str, str]:
    """
    Query Ensembl /lookup/id in batches.
    Returns a dict mapping transcript_id -> gene_id (Parent field).
    """
    transcript_to_gene = {}
    ids = list(set(transcript_ids))

    for batch_start in range(0, len(ids), ENSEMBL_BATCH_SIZE):
        batch = ids[batch_start: batch_start + ENSEMBL_BATCH_SIZE]
        payload = json.dumps({"ids": batch})

        for attempt in range(1, ENSEMBL_RETRY_MAX + 1):
            try:
                resp = requests.post(
                    ENSEMBL_LOOKUP_BATCH_URL,
                    headers=ENSEMBL_HEADERS,
                    data=payload,
                    timeout=30,
                )
                if resp.status_code == 200:
                    result = resp.json()
                    for tid, info in result.items():
                        if info and "Parent" in info:
                            transcript_to_gene[tid] = info["Parent"]
                    break
                elif resp.status_code == 429:
                    # Rate-limited — wait longer
                    retry_after = int(resp.headers.get("Retry-After", ENSEMBL_RETRY_SLEEP_S))
                    print(f"  Rate-limited. Waiting {retry_after}s ...")
                    time.sleep(retry_after)
                else:
                    print(f"  HTTP {resp.status_code} on attempt {attempt}: {resp.text[:200]}")
                    time.sleep(ENSEMBL_RETRY_SLEEP_S)
            except requests.RequestException as e:
                print(f"  Request error on attempt {attempt}: {e}")
                time.sleep(ENSEMBL_RETRY_SLEEP_S)

        time.sleep(ENSEMBL_SLEEP_S)

    return transcript_to_gene


def transcript_to_geneID(df: pd.DataFrame, transcript_col: str) -> pd.DataFrame:
    """
    Add a gene_id column to *df* by looking up each transcript_id via Ensembl.
    The new column is named 'gene_id'.
    """
    unique_transcripts = df[transcript_col].dropna().unique().tolist()
    print(f"  Looking up {len(unique_transcripts)} unique transcript IDs via Ensembl ...")
    mapping = batch_lookup_ensembl(unique_transcripts)
    print(f"  Resolved {len(mapping)} / {len(unique_transcripts)} transcripts to gene IDs.")
    df = df.copy()
    df["gene_id"] = df[transcript_col].map(mapping)
    return df


def _aggregate_predictions(pred_df: pd.DataFrame, gene_id_col: str = "gene_id",
                            hl_col: str = "predicted_halflife") -> pd.Series:
    """
    After mapping to gene IDs, aggregate (mean) duplicate gene entries.
    Returns a Series indexed by gene_id.
    """
    return pred_df.groupby(gene_id_col)[hl_col].mean()


def scatter_with_stats(ax, x, y, xlabel, ylabel, title, color):
    """Draw a scatter plot with regression line, R², and Pearson r."""
    mask = ~(np.isnan(x) | np.isnan(y))
    x_clean, y_clean = x[mask], y[mask]

    ax.scatter(x_clean, y_clean, alpha=0.45, s=18, color=color, edgecolors="none")

    if len(x_clean) > 2:
        slope, intercept, r_value, p_value, _ = stats.linregress(x_clean, y_clean)
        x_line = np.linspace(x_clean.min(), x_clean.max(), 200)
        ax.plot(x_line, slope * x_line + intercept, color="black", linewidth=1.2, linestyle="--")
        r2 = r_value ** 2
        p_str = f"{p_value:.2e}"
        ax.text(0.05, 0.93,
                f"r = {r_value:.3f}   R² = {r2:.3f}   p = {p_str}\nn = {mask.sum()}",
                transform=ax.transAxes, fontsize=8, verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))

    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.xaxis.set_major_locator(ticker.MaxNLocator(6))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(6))


def make_scatter_plots(merged_df: pd.DataFrame, output_path: str):
    """
    Create two scatter plots:
      - human_half_life vs human_pred_hl
      - mouse_half_life vs mouse_pred_hl
    Saved alongside the output CSV as a PNG.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    fig.suptitle("Predicted vs Measured mRNA Half-Life", fontsize=13, fontweight="bold", y=1.01)

    scatter_with_stats(
        axes[0],
        x=merged_df["human_half_life"].values.astype(float),
        y=merged_df["human_pred_hl"].values.astype(float),
        xlabel="Measured half-life (human)",
        ylabel="Predicted half-life (human)",
        title="Human",
        color="#3A86FF",
    )

    scatter_with_stats(
        axes[1],
        x=merged_df["mouse_half_life"].values.astype(float),
        y=merged_df["mouse_pred_hl"].values.astype(float),
        xlabel="Measured half-life (mouse)",
        ylabel="Predicted half-life (mouse)",
        title="Mouse",
        color="#FF6B6B",
    )

    plt.tight_layout()
    plot_path = os.path.splitext(output_path)[0] + "_scatter.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Scatter plots saved to: {plot_path}")


def main(input_dir, tsv_file, output):
    # ------------------------------------------------------------------ #
    # 1. Load input files
    # ------------------------------------------------------------------ #
    print("Loading input files ...")
    df = pd.read_csv(tsv_file, sep='\t')
    # columns: bovine_gene bovine_id human_gene_id human_transcript_id
    #          mouse_return_gene human_half_life mouse_half_life
    #          human_te mouse_te human_RBH mouse_RBH confidence_tier ensembl_gene_id

    human_pred = pd.read_csv(os.path.join(input_dir, 'human/metrics/predictions.tsv'), sep='\t')
    bovine_pred = pd.read_csv(os.path.join(input_dir, 'bovine/metrics/predictions.tsv'), sep='\t')
    mouse_pred  = pd.read_csv(os.path.join(input_dir, 'mouse/metrics/predictions.tsv'), sep='\t')
    # columns for each: transcript_id  predicted_halflife

    # Strip version suffixes (e.g. ENST00000428931.6 -> ENST00000428931)
    for pred_df in (human_pred, bovine_pred, mouse_pred):
        pred_df["transcript_id"] = pred_df["transcript_id"].str.split(".").str[0]

    # ------------------------------------------------------------------ #
    # 2. Map transcript IDs → gene IDs for bovine and mouse predictions
    #    (human predictions already use transcript IDs, which map directly
    #     to human_transcript_id in the TSV)
    # ------------------------------------------------------------------ #

    # --- Bovine ---
    print("\nMapping bovine transcripts → gene IDs ...")
    bovine_pred = transcript_to_geneID(bovine_pred, "transcript_id")
    bov_series = _aggregate_predictions(bovine_pred, gene_id_col="gene_id")

    # --- Mouse ---
    print("\nMapping mouse transcripts → gene IDs ...")
    mouse_pred = transcript_to_geneID(mouse_pred, "transcript_id")
    mouse_series = _aggregate_predictions(mouse_pred, gene_id_col="gene_id")

    # --- Human (transcript-level join, no API call needed) ---
    # human_pred has transcript_id; TSV has human_transcript_id
    human_pred_series = human_pred.set_index("transcript_id")["predicted_halflife"]

    # ------------------------------------------------------------------ #
    # 3. Merge predicted half-lives into the main TSV
    # ------------------------------------------------------------------ #
    print("\nMerging predictions ...")

    df["bov_pred_hl"]   = df["ensembl_gene_id"].map(bov_series)
    df["human_pred_hl"] = df["human_transcript_id"].map(human_pred_series)
    df["mouse_pred_hl"] = df["mouse_return_gene"].map(mouse_series)

    for label, col, pred_series, tsv_key in [
        ("Bovine", "bov_pred_hl",   bov_series,        "ensembl_gene_id"),
        ("Mouse",  "mouse_pred_hl", mouse_series,      "mouse_return_gene"),
        ("Human",  "human_pred_hl", human_pred_series, "human_transcript_id"),
    ]:
        n_matched    = df[col].notna().sum()
        n_total_tsv  = len(df)
        tsv_ids      = set(df[tsv_key].dropna())
        pred_ids     = set(pred_series.index)
        n_pred_total = len(pred_ids)
        n_in_both    = len(tsv_ids & pred_ids)
        n_pred_only  = len(pred_ids - tsv_ids)
        n_tsv_only   = len(tsv_ids - pred_ids)
        print(
            f"\n  {label}:"
            f"\n    TSV rows          : {n_total_tsv}"
            f"\n    Prediction genes  : {n_pred_total}"
            f"\n    In both           : {n_in_both}"
            f"\n    Merged (non-null) : {n_matched}"
            f"\n    Only in preds     : {n_pred_only}  (predicted but no TSV row)"
            f"\n    Only in TSV       : {n_tsv_only}  (TSV row but no prediction)"
        )

    # ------------------------------------------------------------------ #
    # 4. Save enriched TSV
    # ------------------------------------------------------------------ #
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    df.to_csv(output, index=False)
    print(f"\nEnriched table saved to: {output}")

    # ------------------------------------------------------------------ #
    # 5. Scatter plots
    # ------------------------------------------------------------------ #
    print("\nGenerating scatter plots ...")
    make_scatter_plots(df, output)


if __name__ == "__main__":
    args = parse_args()
    main(args.input_dir, args.tsv_file, args.output)