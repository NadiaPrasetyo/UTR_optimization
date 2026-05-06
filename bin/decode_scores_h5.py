#!/usr/bin/env python3

import argparse
import os
import h5py
import numpy as np
import pandas as pd


NUCS = ['A', 'C', 'G', 'T']

def safe_get(h5, key):    
    if key in h5:
        return h5[key][:]
    else:
        print(f"[WARNING] Missing dataset: '{key}'")
        return None


def parse_fasta_ids(fasta_file):
    """Extract sequence IDs from FASTA headers."""
    ids = []
    if fasta_file is None:
        return None

    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                ids.append(line[1:].strip().split()[0])
    return ids

def decode_preds(h5, seq_ids, out_dir):
    """
    Decode model predictions from HDF5 into a tabular CSV.

    OUTPUT COLUMNS:

    seq_id:
        Input sequence identifier from FASTA header.

    target:
        Index of model output head.
        (0 = primary half-life prediction in most setups)

    pred_half_life_mean:
        Mean prediction across ensemble models.
        NOTE: This is often in transformed space (e.g. log-scale or normalized),
        NOT raw physical half-life unless explicitly calibrated.

    pred_half_life_std:
        Standard deviation across ensemble models.
        Measures epistemic uncertainty:
        - low (~0.0–0.2): high confidence
        - moderate (~0.2–0.8): partial disagreement
        - high (>1.0): unstable prediction

    pred_half_life_min / max:
        Range of ensemble predictions (best/worst model).

    pred_half_life_ci_lower / ci_upper:
        Approximate 95% confidence interval assuming Gaussian ensemble spread:
            mean ± 1.96 * std

    ensemble_agreement:
        Heuristic confidence score derived from coefficient of variation:
            cv = std / (|mean| + 1e-8)
            agreement = 1 / (1 + cv)

        Higher = more consistent ensemble predictions.
        WARNING: depends on magnitude of mean; interpret carefully if mean ~ 0.

    stability_class:
        Heuristic binning of mean prediction:
            mean < 2   → very_unstable
            mean < 5   → unstable
            mean < 10  → moderate
            else       → stable

        NOTE: thresholds apply in MODEL SPACE, not necessarily real time units.
    """

    preds = safe_get(h5, 'preds')
    preds_mean = safe_get(h5, 'preds_mean')

    if preds is None and preds_mean is None:
        print("[WARNING] No prediction data found. Skipping predictions.")
        return

    rows = []

    if preds is not None:
        num_seqs, num_targets, num_models = preds.shape
    else:
        num_seqs, num_targets = preds_mean.shape
        num_models = 0

    for i in range(num_seqs):
        seq_id = seq_ids[i] if seq_ids and i < len(seq_ids) else f"seq_{i}"

        for t in range(num_targets):
            model_vals = None

            if preds is not None:
                model_vals = preds[i, t, :]
                mean = np.mean(model_vals)
                std = np.std(model_vals)
                min_v = np.min(model_vals)
                max_v = np.max(model_vals)
            else:
                mean = preds_mean[i, t]
                std = np.nan
                min_v = np.nan
                max_v = np.nan

            # fallback if preds_mean exists but preds missing
            if preds_mean is not None:
                mean = preds_mean[i, t]

            # CI only if std is valid
            if not np.isnan(std):
                ci_lower = mean - 1.96 * std
                ci_upper = mean + 1.96 * std
                cv = std / (np.abs(mean) + 1e-8)
                agreement = 1 / (1 + cv)
            else:
                ci_lower = ci_upper = np.nan
                agreement = np.nan

            # biological label (in model space)
            if mean < 2:
                stability_class = "very_unstable"
            elif mean < 5:
                stability_class = "unstable"
            elif mean < 10:
                stability_class = "moderate"
            else:
                stability_class = "stable"

            rows.append({
                "seq_id": seq_id,
                "target": t,
                "pred_log10_half_life_mean": mean,
                "pred_log10_half_life_std": std,
                "pred_log10_half_life_min": min_v,
                "pred_log10_half_life_max": max_v,
                "pred_log10_half_life_ci_lower": ci_lower,
                "pred_log10_half_life_ci_upper": ci_upper,
                "ensemble_agreement": agreement,
                "stability_class": stability_class
            })

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out_dir, "predictions.csv"), index=False)


def decode_grads(h5, seq_ids, out_dir, write_raw=False):
    seqs = safe_get(h5, 'seqs')
    def get_annotation(seq_i, pos):
        # CDS mask from channel 4
        cds_mask = seqs[seq_i, :, 4].astype(bool)

        cds_positions = np.where(cds_mask)[0]

        # no CDS detected → treat everything as UTR
        if len(cds_positions) == 0:
            return "utr"

        cds_start = cds_positions[0]
        cds_end = cds_positions[-1]

        if pos < cds_start:
            return "5'_utr"
        elif pos > cds_end:
            return "3'_utr"
        else:
            return "cds"

    grads = safe_get(h5, 'grads')

    if grads is None:
        print("[WARNING] No gradient data found. Skipping gradients.")
        return

    if grads.ndim != 4:
        print(f"[WARNING] Unexpected grads shape: {grads.shape}")
        return

    num_seqs, seq_len, num_channels, num_models = grads.shape

    if num_channels < 4:
        print("[WARNING] Not enough channels for nucleotide gradients.")
        return

    summary_rows = []
    raw_rows = []

    for i in range(num_seqs):
        seq_id = seq_ids[i] if seq_ids and i < len(seq_ids) else f"seq_{i}"

        for pos in range(seq_len):
            try:
                g = grads[i, pos, :, :].mean(axis=-1)
            except Exception as e:
                print(f"[WARNING] Failed at seq {i}, pos {pos}: {e}")
                continue

            # nucleotide from one-hot
            if seqs is not None:
                nuc_idx = np.argmax(seqs[i, pos, :4])
                nuc = NUCS[nuc_idx]
            else:
                nuc = None

            annotation = get_annotation(i, pos) if seqs is not None else None

            nuc_grads = g[:4]

            if np.all(np.isnan(nuc_grads)):
                best_nuc = worst_nuc = None
                best_grad = worst_grad = np.nan
            else:
                best_idx = np.nanargmax(nuc_grads)
                worst_idx = np.nanargmin(nuc_grads)
                best_nuc = NUCS[best_idx]
                worst_nuc = NUCS[worst_idx]
                best_grad = nuc_grads[best_idx]
                worst_grad = nuc_grads[worst_idx]

            row = {
                "seq_id": seq_id,
                "position": pos,
                "nuc": nuc,
                "annotation": annotation,
                "grad_mean": np.nanmean(nuc_grads),
                "grad_std": np.nanstd(nuc_grads),
                "best_nuc": best_nuc,
                "best_grad": best_grad,
                "worst_nuc": worst_nuc,
                "worst_grad": worst_grad,
                "delta_best_worst": (
                    best_grad - worst_grad
                    if not np.isnan(best_grad) and not np.isnan(worst_grad)
                    else np.nan
                ),
                "frame_flag": g[4] if num_channels > 4 else np.nan,
                "special_flag": g[5] if num_channels > 5 else np.nan
            }

            summary_rows.append(row)

            if write_raw:
                for m in range(num_models):
                    raw = grads[i, pos, :, m]

                    raw_rows.append({
                        "seq_id": seq_id,
                        "position": pos,
                        "nuc": nuc,
                        "annotation": annotation,
                        "model": m,
                        "A_grad": raw[0] if num_channels > 0 else np.nan,
                        "C_grad": raw[1] if num_channels > 1 else np.nan,
                        "G_grad": raw[2] if num_channels > 2 else np.nan,
                        "T_grad": raw[3] if num_channels > 3 else np.nan,
                        "frame_flag": raw[4] if num_channels > 4 else np.nan,
                        "special_flag": raw[5] if num_channels > 5 else np.nan
                    })

    pd.DataFrame(summary_rows).to_csv(
        os.path.join(out_dir, "gradients_summary.csv"), index=False
    )

    if write_raw:
        pd.DataFrame(raw_rows).to_csv(
            os.path.join(out_dir, "gradients_raw.csv"), index=False
        )


def main():
    parser = argparse.ArgumentParser(
        description="Decode Saluki/Basenji H5 output into human-readable CSVs"
    )
    parser.add_argument("h5_file", help="Input H5 file (scores.h5)")
    parser.add_argument(
        "-o", "--out_dir", default="csv_out",
        help="Output directory for CSVs [Default: csv_out]"
    )
    parser.add_argument(
        "-f", "--fasta_file", default=None,
        help="Original FASTA file to recover sequence IDs"
    )
    parser.add_argument(
        "--raw-grads", action="store_true",
        help="Write raw gradient CSV (can be very large)"
    )

    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    seq_ids = parse_fasta_ids(args.fasta_file)

    with h5py.File(args.h5_file, "r") as h5:
        print("Decoding predictions...")
        decode_preds(h5, seq_ids, args.out_dir)

        print("Decoding gradients...")
        decode_grads(h5, seq_ids, args.out_dir, args.raw_grads)

    print(f"Done. Output written to {args.out_dir}/")


if __name__ == "__main__":
    main()