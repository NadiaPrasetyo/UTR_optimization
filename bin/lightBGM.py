import os
import sys
import pandas as pd
import lightgbm as lgb
import argparse

FEATURE_COLS_PATH = "lightGBM/halflife_lgbm_v4_feature_cols.csv"
MODEL_PATH = "lightGBM/halflife_lgbm_v4_lightgbm.txt"
JOIN_KEY = "transcript_id"

def engineer_features(df: pd.DataFrame) -> pd.DataFrame:
    """Compute derived features missing from raw TSVs."""

    # --- junction densities: normalise counts by region length ---
    # lengths come from sequence_basic (in nucleotides)
    for region, count_col, length_col in [
        ("5utr",  "n_5utr_junctions",  "length_5utr"),   # need to carry length through
        ("cds",   "n_cds_junctions",   "length_cds"),
        ("3utr",  "n_3utr_junctions",  "length_3utr"),
    ]:
        density_col = f"junctions_density_{region}"
        if count_col in df.columns and length_col in df.columns:
            df[density_col] = df[count_col] / df[length_col].replace(0, float("nan"))
        else:
            print(f"  Warning: cannot compute {density_col} — missing {count_col} or {length_col}")

    if all(c in df.columns for c in ["junctions_density_5utr", "junctions_density_cds", "junctions_density_3utr"]):
        df["junctions_density_mrna"] = (
            df["junctions_density_5utr"].fillna(0) +
            df["junctions_density_cds"].fillna(0) +
            df["junctions_density_3utr"].fillna(0)
        )

    # --- exon length fractions: first/last exon as fraction of total mRNA length ---
    if "exon_length_first_mrna" in df.columns and "length_mrna" in df.columns:
        df["exon_length_first_fraction_mrna"] = (
            df["exon_length_first_mrna"] / df["length_mrna"].replace(0, float("nan"))
        )
    if "exon_length_last_mrna" in df.columns and "length_mrna" in df.columns:
        df["exon_length_last_fraction_mrna"] = (
            df["exon_length_last_mrna"] / df["length_mrna"].replace(0, float("nan"))
        )

    return df

def normalise_columns(df: pd.DataFrame, fname: str) -> pd.DataFrame:
    """
    Rename TSV-specific columns to match the model's expected feature names.
    Handles each TSV's quirks individually based on filename.
    """
    name = os.path.splitext(fname)[0]  # e.g. "sequence_basic"

    # --- sequence_basic: wide → pivoted-style rename using 'region' column ---
    # Has columns: gc_content, frac_A, frac_C, ... with a 'region' column
    # These need to become gc_content_3utr, frac_a_3utr, etc.
    if name == "sequence_basic":
        region_map = {
            "3utr": "3utr", "5utr": "5utr", "cds": "cds",
            "last100": "last100", "mrna": "mrna", "start": "start", "stop": "stop"
        }
        metric_cols = ["gc_content", "frac_A", "frac_C", "frac_G", "frac_U",
                    "frac_other", "gc_skew", "at_skew", "purine_ratio", "amino_ratio"]

        frames = []
        for region, region_suffix in region_map.items():
            sub = df[df["region"].str.lower() == region][["transcript_id", "length"] + metric_cols].copy()
            sub = sub.rename(columns={
                "length": f"length_{region_suffix}",
                **{c: f"{c.lower()}_{region_suffix}" for c in metric_cols}
            })
            frames.append(sub.set_index("transcript_id"))

        df = pd.concat(frames, axis=1).reset_index()
        return df

    # --- stopfree: pivot by region (3utr, 5utr, mrna) ---
    if name == "stopfree":
        frames = []
        for region in ["3utr", "5utr", "mrna"]:
            sub = df[df["region"].str.lower() == region][["transcript_id", "stopfree_fraction"]].copy()
            sub = sub.rename(columns={"stopfree_fraction": f"stopfree_{region}"})
            frames.append(sub.set_index("transcript_id"))
        df = pd.concat(frames, axis=1).reset_index()
        return df

    # --- architecture ---
    if name == "architecture":
        df = df.rename(columns={
            "first_exon_length":  "exon_length_first_mrna",
            "last_exon_length":   "exon_length_last_mrna",
            "intron_mean":        "intron_length_mean_mrna",
            "intron_median":      "intron_median",
            "intron_sd":          "intron_sd",
        })
        return df

    # --- codon_aa_counts: lowercase codon/aa names, drop non-feature cols ---
    if name == "codon_aa_counts":
        rename = {}
        for col in df.columns:
            if col.startswith("codon_"):
                rename[col] = col.lower() + "_cds"        # codon_AAA → codon_aaa_cds
            elif col.startswith("aa_"):
                rename[col] = col.lower() + "_cds"        # aa_A → aa_a_cds
        df = df.rename(columns=rename)
        # drop helper cols not in feature list
        df = df.drop(columns=[c for c in ["cds_length_codons", "n_codons_scored",
                                           "n_stops", "aa_total"] if c in df.columns])
        return df

    # --- junctions ---
    if name == "junctions":
        df = df.rename(columns={
            # keep raw counts for density calculation post-merge
            "n_5utr_junctions":              "n_5utr_junctions",
            "n_cds_junctions":               "n_cds_junctions",
            "n_3utr_junctions":              "n_3utr_junctions",
            "stop_dist_closest_upstream":    "stop_dist_closest_upstream",
            "stop_dist_closest_downstream":  "stop_dist_closest_downstream",
            "stop_dist_last_downstream":     "stop_dist_last_downstream",
            "start_dist_closest_upstream":   "start_dist_closest_upstream",
            "start_dist_closest_downstream": "start_dist_closest_downstream",
        })
        return df

    # --- nmd_fragility ---
    if name == "nmd_fragility_full":
        df = df.rename(columns={
            "transition_fragile_codon_density":    "nmd_transition_fragile_codon_density_mrna",
            "transversion_fragile_codon_density":  "nmd_transversion_fragile_codon_density_mrna",
            "snv_fragile_codon_density":           "nmd_snv_fragile_codon_density_mrna",
            "alt_stop_codon_density":              "nmd_alt_stop_codon_density_mrna",
            "transition_fraction_of_snv_fragile":  "nmd_transition_fraction_of_snv_fragile_mrna",
        })
        return df

    # --- uorf ---
    if name == "uorf":
        df = df.rename(columns={
            "has_uorf":                      "uorf_present_mrna",
            "dist_cap_to_first_uatg":        "dist_cap_to_first_uatg_mrna",
            "n_overlapping_uorfs":           "n_overlapping_uorfs",
            "total_classical_uorf_codons":   "total_classical_uorf_codons",
            "max_classical_uorf_codons":     "max_classical_uorf_codons",
            "dist_last_uorf_stop_to_main_atg": "dist_last_uorf_stop_to_main_atg",
        })
        return df

    return df  # fallback: return as-is

def load_and_merge_tsvs(input_dir: str) -> pd.DataFrame:
    """Load all TSV files in input_dir and merge them on transcript_id."""
    tsv_files = sorted(f for f in os.listdir(input_dir) if f.endswith(".tsv"))
    if not tsv_files:
        sys.exit(f"Error: no .tsv files found in '{input_dir}'")

    merged = None
    shared_cols = None

    for fname in tsv_files:
        path = os.path.join(input_dir, fname)
        df = pd.read_csv(path, sep="\t")
        df = normalise_columns(df, fname)

        # normalise column names to snake_case
        df.columns = (
            df.columns
            .str.strip()
            .str.lower()
            .str.replace(r"[ \-]+", "_", regex=True)   # spaces/hyphens → underscore
            .str.replace(r"[^\w]", "", regex=True)      # strip any remaining non-word chars
        )

        if JOIN_KEY not in df.columns:
            sys.exit(f"Error: '{JOIN_KEY}' column not found in {fname} (after normalisation)")

        if merged is None:
            merged = df
            shared_cols = set(df.columns)
        else:
            duplicate_cols = [c for c in df.columns if c in shared_cols and c != JOIN_KEY]
            if duplicate_cols:
                print(f"  Dropping duplicate cols from {fname}: {duplicate_cols}")
            df = df.drop(columns=duplicate_cols)
            merged = merged.merge(df, on=JOIN_KEY, how="inner")
            shared_cols.update(df.columns)

        print(f"  Loaded {fname}  ({len(df):,} rows, {df.shape[1]} cols)")

    print(f"\nMerged dataset: {len(merged):,} transcripts, {merged.shape[1]} columns")
    return merged


def main(input_dir: str) -> None:
    print(f"Reading TSV files from: {input_dir}")
    df = load_and_merge_tsvs(input_dir)
    df = engineer_features(df) # Calculate derived features

    # --- load feature list ---
    if not os.path.exists(FEATURE_COLS_PATH):
        sys.exit(f"Error: feature column file not found: {FEATURE_COLS_PATH}")
    feature_cols = pd.read_csv(FEATURE_COLS_PATH)["feature"].tolist()

    missing = [c for c in feature_cols if c not in df.columns]
    if missing:
        sys.exit(f"Error: {len(missing)} feature(s) missing from merged data:\n  " + "\n  ".join(missing))

    X = df[feature_cols]

    # --- load model & predict ---
    if not os.path.exists(MODEL_PATH):
        sys.exit(f"Error: model file not found: {MODEL_PATH}")
    model = lgb.Booster(model_file=MODEL_PATH)

    model_features = model.feature_name()
    print(f"Model expects {len(model_features)} features:")
    print("\n".join(model_features))

    missing = [c for c in feature_cols if c not in model_features]
    if missing:
        sys.exit(f"Error: {len(missing)} feature(s) missing from model:\n  " + "\n  ".join(missing))

    print("Running predictions...")
    preds = model.predict(X)

    # --- save output ---
    out_df = pd.DataFrame({
        JOIN_KEY: df[JOIN_KEY],
        "predicted_halflife": preds
    })
    out_path = os.path.join(input_dir, "predictions.tsv")
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"Predictions saved to: {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run LightGBM half-life predictions")
    parser.add_argument("-i", "--input", required=True, type=str,
                        help="Directory containing input .tsv metric files")
                        
    args = parser.parse_args()
    main(args.input)