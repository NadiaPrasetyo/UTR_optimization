# a genetic algorithm to optimize the sequence of a specific target: 3' UTR
import os
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import random
import logging
import subprocess
import shutil

reference_gfp_cds = SeqIO.read('data/GFP_CDS.fasta', 'fasta').seq
reference_pfizer_5UTR = SeqIO.read('data/Pfizer_5UTR.fasta', 'fasta').seq

SALUKI_PREDICT = 'bin/saluki_predict_fasta.py'
DECODE_SCORES  = 'bin/decode_scores_h5.py'
MODELS_DIR     = 'models/'   # adjust to your models directory
NUCS = list('ACGT')


def setup_logging(verbose=False, output_dir='data/'):
    log_level = logging.DEBUG if verbose else logging.INFO
    log_file = os.path.join(output_dir, "mutagenesis_saluki.log")
    os.makedirs(output_dir, exist_ok=True)
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode='w') if verbose else logging.NullHandler(),
            logging.StreamHandler(sys.stdout)
        ]
    )
    if verbose:
        logging.info("Logging initialized. Log file: %s", log_file)


def mutate(seq, i):
    """Yield all single-nucleotide substitutions at position i (excluding the original base)."""
    original = seq[i]
    for base in NUCS:
        if base != original:
            new_seq = list(seq)
            new_seq[i] = base
            yield ''.join(new_seq), f"{original}{i+1}{base}"


def write_fasta(sequences, fasta_path):
    """
    Write sequences to a FASTA file.
    sequences: list of (seq_id, seq_str) tuples
    """
    with open(fasta_path, 'w') as f:
        for seq_id, seq_str in sequences:
            f.write(f">{seq_id}\n{seq_str}\n")


def run_saluki_predict(fasta_path, out_dir, models_dir=MODELS_DIR):
    """Run saluki_predict_fasta.py on a multi-sequence FASTA and return the path to scores.h5."""
    h5_path = os.path.join(out_dir, 'scores.h5')
    cmd = [
        sys.executable, SALUKI_PREDICT,
        '-o', out_dir,
        '-d', '0',
        models_dir,
        fasta_path          # all candidates in one shot
    ]
    logging.debug("Running: %s", ' '.join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error("saluki_predict_fasta.py failed:\n%s", result.stderr)
        raise RuntimeError("Prediction step failed.")
    logging.debug(result.stdout)
    return h5_path


def run_decode_scores(h5_path, fasta_path, out_dir):
    """Run decode_scores_h5.py and return the predictions CSV path."""
    csv_path = os.path.join(out_dir, 'predictions.csv')
    cmd = [
        sys.executable, DECODE_SCORES,
        '-o', out_dir,
        '-f', fasta_path,
        h5_path
    ]
    logging.debug("Running: %s", ' '.join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error("decode_scores_h5.py failed:\n%s", result.stderr)
        raise RuntimeError("Decoding step failed.")
    logging.debug(result.stdout)
    return csv_path


def select(candidates, fasta_path, work_dir, models_dir=MODELS_DIR):
    """
    Score ALL candidate sequences in a single Saluki call and return the best one.

    candidates : list of (seq_id, full_seq_str) tuples — written as one FASTA
    fasta_path : path where the multi-sequence FASTA will be written
    work_dir   : scratch directory for this selection round (h5 + csvs land here)

    Returns: (best_seq_id, best_seq_str, best_score, scores_df)
    """
    os.makedirs(work_dir, exist_ok=True)
    write_fasta(candidates, fasta_path)          # one FASTA, N sequences

    h5_path  = run_saluki_predict(fasta_path, work_dir, models_dir)   # single call
    csv_path = run_decode_scores(h5_path, fasta_path, work_dir)

    df = pd.read_csv(csv_path)
    df_t0 = df[df['target'] == 0].copy().set_index('seq_id')

    best_id    = df_t0['pred_log10_half_life_mean'].idxmax()
    best_score = df_t0.loc[best_id, 'pred_log10_half_life_mean']
    best_seq   = next(seq for sid, seq in candidates if sid == best_id)

    logging.info("  Best candidate: %s  score=%.4f", best_id, best_score)
    return best_id, best_seq, best_score, df_t0


def genetic_algorithm(
    utr_seq,
    full_seq_builder,
    output_dir,
    models_dir=MODELS_DIR,
    num_generations=10,
    mutations_per_generation=5,
    save_cycle=False,
    cycle_dir=None,
):
    """
    Simple single-site greedy genetic algorithm.

    Each generation:
      1. Pick `mutations_per_generation` random positions in the 3' UTR.
      2. Generate all 3 possible substitutions at each picked position
         → up to 3 × mutations_per_generation candidates per generation.
      3. Score every candidate with Saluki.
      4. Keep the best-scoring sequence as the parent for the next generation.

    full_seq_builder : callable(utr_str) -> full_seq_str
                       Prepends CDS + 5'UTR to the 3' UTR variant.
    """
    current_utr   = str(utr_seq)
    current_score = None
    history = []

    for gen in range(num_generations):
        logging.info("=== Generation %d / %d ===", gen + 1, num_generations)

        utr_len   = len(current_utr)
        positions = random.sample(range(utr_len), min(mutations_per_generation, utr_len))
        logging.debug("Mutating positions: %s", positions)

        # Build one candidate list with current + all mutations
        candidates = [("current", full_seq_builder(current_utr))]
        for pos in positions:
            for mut_utr, mut_label in mutate(current_utr, pos):
                candidates.append((f"gen{gen+1}_{mut_label}", full_seq_builder(mut_utr)))

        logging.info("  Evaluating %d candidates in a single Saluki call …", len(candidates))

        gen_work   = os.path.join(output_dir, f"_gen_{gen+1:03d}_scratch")
        fasta_path = os.path.join(gen_work, "candidates.fasta")

        best_id, best_seq_full, best_score, scores_df = select(
            candidates, fasta_path, gen_work, models_dir
        )

        prefix_len = len(full_seq_builder(""))
        best_utr   = best_seq_full[prefix_len:]
        mutation_applied = best_id if best_id != "current" else "none"

        improved = current_score is None or best_score > current_score
        if improved:
            logging.info(
                "  Improvement: %.4f → %.4f  (%s)",
                current_score if current_score is not None else float('nan'),
                best_score, mutation_applied
            )
            current_utr   = best_utr
            current_score = best_score
        else:
            logging.info(
                "  No improvement (best=%.4f, current=%.4f). Keeping current sequence.",
                best_score, current_score
            )

        history.append({
            "generation":        gen + 1,
            "best_candidate_id": best_id,
            "mutation_applied":  mutation_applied,
            "score":             best_score,
            "improved":          improved,
            "current_utr":       current_utr,
        })

        if save_cycle and cycle_dir:
            dest = os.path.join(cycle_dir, f"gen_{gen+1:03d}")
            os.makedirs(dest, exist_ok=True)

            # Copy CSVs
            for fname in ("predictions.csv", "gradients_summary.csv"):
                src = os.path.join(gen_work, fname)
                if os.path.exists(src):
                    shutil.copy2(src, os.path.join(dest, fname))

            # Copy raw H5
            h5_src = os.path.join(gen_work, "scores.h5")
            if os.path.exists(h5_src):
                shutil.copy2(h5_src, os.path.join(dest, "scores.h5"))

            # Save winning 3' UTR for this generation
            write_fasta(
                [(f"gen{gen+1}_best_utr", current_utr)],
                os.path.join(dest, "best_utr.fasta")
            )
            logging.debug("  Saved cycle outputs to %s", dest)

        shutil.rmtree(gen_work, ignore_errors=True)   # always clean scratch

    return current_utr, current_score, pd.DataFrame(history)


def main(input_file, output_dir, models_dir, num_generations,
         mutations_per_generation, save_cycle, verbose):

    utr_seq = SeqIO.read(input_file, 'fasta').seq
    logging.info("Loaded 3' UTR: %d nt from %s", len(utr_seq), input_file)
    logging.info("GFP CDS     : %d nt", len(reference_gfp_cds))
    logging.info("Pfizer 5'UTR: %d nt", len(reference_pfizer_5UTR))

    prefix = str(reference_gfp_cds) + str(reference_pfizer_5UTR)

    def full_seq_builder(utr_str):
        return prefix + utr_str

    cycle_dir = None
    if save_cycle:
        cycle_dir = os.path.join(output_dir, "saluki_cycles")
        os.makedirs(cycle_dir, exist_ok=True)

    # Score the original sequence first as a reference
    logging.info("Scoring original sequence …")
    orig_work  = os.path.join(output_dir, "_original_scratch")
    orig_fasta = os.path.join(orig_work, "original.fasta")
    _, _, orig_score, _ = select(
        [("original", full_seq_builder(str(utr_seq)))],
        orig_fasta, orig_work, models_dir
    )
    logging.info("Original sequence score: %.4f", orig_score)
    if not save_cycle:
        shutil.rmtree(orig_work, ignore_errors=True)

    # Run the genetic algorithm
    best_utr, best_score, history_df = genetic_algorithm(
        utr_seq          = utr_seq,
        full_seq_builder = full_seq_builder,
        output_dir       = output_dir,
        models_dir       = models_dir,
        num_generations  = num_generations,
        mutations_per_generation = mutations_per_generation,
        save_cycle       = save_cycle,
        cycle_dir        = cycle_dir,
    )

    # ── Outputs ──────────────────────────────────────────────────────────────
    # 1. Optimised 3' UTR
    best_utr_fasta = os.path.join(output_dir, "best_utr_optimized.fasta")
    write_fasta([("optimized_3utr", best_utr)], best_utr_fasta)
    logging.info("Best 3' UTR written to %s", best_utr_fasta)

    # 2. Full optimised sequence
    best_full_fasta = os.path.join(output_dir, "best_full_sequence.fasta")
    write_fasta([("optimized_full", full_seq_builder(best_utr))], best_full_fasta)
    logging.info("Best full sequence written to %s", best_full_fasta)

    # 3. Generation history CSV
    history_csv = os.path.join(output_dir, "optimization_history.csv")
    history_df.to_csv(history_csv, index=False)
    logging.info("Optimization history written to %s", history_csv)

    # 4. Summary
    logging.info(
        "\n=== Optimization complete ===\n"
        "  Original score : %.4f\n"
        "  Best score     : %.4f\n"
        "  Improvement    : %+.4f\n"
        "  Generations    : %d",
        orig_score, best_score, best_score - orig_score, num_generations
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Optimize the sequence of the 3' UTR based on Saluki predictions"
    )
    parser.add_argument('-i', '--input-file', required=True, type=str,
                        help="Input FASTA file containing the 3' UTR sequence")
    parser.add_argument('-o', '--output-dir', default='data/', type=str,
                        help="Output directory")
    parser.add_argument('-m', '--models-dir', default='models/', type=str,
                        help="Directory containing Saluki fold/cross models")
    parser.add_argument('--num-generations', type=int, default=10,
                        help="Number of GA generations [Default: 10]")
    parser.add_argument('--mutations-per-generation', type=int, default=5,
                        help="Random positions mutated per generation [Default: 5]")
    parser.add_argument('--save-cycle', action='store_true',
                        help="Save intermediate cycle outputs and sequences")
    parser.add_argument('--verbose', action='store_true',
                        help="Enable verbose logging")
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    setup_logging(verbose=args.verbose, output_dir=args.output_dir)
    main(
        input_file               = args.input_file,
        output_dir               = args.output_dir,
        models_dir               = args.models_dir,
        num_generations          = args.num_generations,
        mutations_per_generation = args.mutations_per_generation,
        save_cycle               = args.save_cycle,
        verbose                  = args.verbose,
    )
