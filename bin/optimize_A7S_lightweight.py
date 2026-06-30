#!/usr/bin/env python3
"""
mutate_3utr_ga.py
─────────────────────────────────────────────────────────────────────────────
Genetic Algorithm to evolve 3' UTR sequences that maximise predicted mRNA
half-life, using the 01b_metrics.py + LightGBM prediction pipeline.

USAGE
─────
python mutate_3utr_ga.py \\
    --fasta-5utr   input_5utr.fa      \\   # single sequence
    --fasta-cds    input_cds.fa        \\   # single sequence
    --fasta-3utr   seed_3utr.fa        \\   # single sequence (seed)
    --species      "Homo sapiens"      \\
    --metrics-script  bin/01b_metrics.py \\
    --predict-script  predict_halflife.py \\
    --output-dir   ga_output           \\
    [--population  50]                 \\   # individuals per generation
    [--generations 30]                 \\   # number of GA iterations
    [--mutation-rate 0.02]             \\   # per-nucleotide mutation probability
    [--crossover-rate 0.5]             \\   # crossover probability
    [--elite-frac  0.1]                \\   # fraction kept as elites (no mutation)
    [--tournament-k 3]                 \\   # tournament selection size
    [--seed        42]                 \\
    [--no-plot]                        \\   # skip generating the fitness plot

WHAT IT DOES
────────────
1.  Reads the fixed 5' UTR and CDS (these are never mutated).
2.  Seeds a population of N copies of the initial 3' UTR, then applies
    random point mutations to every non-elite member each generation.
3.  Each generation:
      a.  Writes three temporary FASTA files (5UTR, CDS, 3UTR variants).
      b.  Calls 01b_metrics.py --fasta-* to compute feature TSVs.
      c.  Calls the LightGBM prediction script to get half-life scores.
      d.  Ranks by predicted half-life, keeps elites, breeds next generation
          via tournament selection + uniform crossover + point mutation.
4.  Writes results/  containing:
      - best_per_generation.tsv      — top score each generation
      - all_generations.tsv          — every evaluated individual
      - best_3utr.fa                 — FASTA of the all-time best sequence
      - evolution_summary.txt        — human-readable run summary
      - fitness_over_generations.png — plot of best/mean/worst half-life per gen
─────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('ga_3utr')

NUCLEOTIDES = ['A', 'T', 'G', 'C']


# ══════════════════════════════════════════════════════════════════════════════
# FASTA I/O
# ══════════════════════════════════════════════════════════════════════════════

def read_single_fasta(path: Path) -> Tuple[str, str]:
    """Read a FASTA with exactly one record → (seq_id, sequence)."""
    seq_id: Optional[str] = None
    parts: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if seq_id is not None:
                    raise ValueError(f"Multiple sequences found in {path}; expected exactly one.")
                seq_id = line[1:].split()[0]
            elif line:
                parts.append(line.upper().replace('U', 'T'))
    if seq_id is None:
        raise ValueError(f"No sequences found in {path}")
    return seq_id, ''.join(parts)


def write_fasta(path: Path, records: List[Tuple[str, str]], line_width: int = 60) -> None:
    """Write (seq_id, sequence) records to a FASTA file."""
    with open(path, 'w') as fh:
        for seq_id, seq in records:
            fh.write(f'>{seq_id}\n')
            for i in range(0, len(seq), line_width):
                fh.write(seq[i:i + line_width] + '\n')


# ══════════════════════════════════════════════════════════════════════════════
# Genetic operators
# ══════════════════════════════════════════════════════════════════════════════

def point_mutate(seq: str, mutation_rate: float, rng: random.Random) -> str:
    """Apply random point mutations at *mutation_rate* per nucleotide."""
    bases = list(seq)
    for i, base in enumerate(bases):
        if rng.random() < mutation_rate:
            alt = [n for n in NUCLEOTIDES if n != base]
            bases[i] = rng.choice(alt)
    return ''.join(bases)


def uniform_crossover(seq_a: str, seq_b: str, rng: random.Random) -> str:
    """
    Produce one offspring via uniform crossover.
    Each position is drawn independently from either parent with p=0.5.
    Both parents must be the same length.
    """
    if len(seq_a) != len(seq_b):
        # Length mismatch: just return one parent (safe fallback)
        return seq_a
    return ''.join(
        a if rng.random() < 0.5 else b
        for a, b in zip(seq_a, seq_b)
    )


def tournament_select(population: List[str], scores: List[float],
                       k: int, rng: random.Random) -> str:
    """Select one individual via k-way tournament (highest score wins)."""
    contestants = rng.sample(range(len(population)), min(k, len(population)))
    winner = max(contestants, key=lambda i: scores[i])
    return population[winner]


# ══════════════════════════════════════════════════════════════════════════════
# Pipeline: metrics → prediction
# ══════════════════════════════════════════════════════════════════════════════

def run_metrics(
    metrics_script: Path,
    fasta_5utr: Path,
    fasta_cds: Path,
    fasta_3utr: Path,
    species: str,
    output_dir: Path,
    force: bool = True,
) -> bool:
    """
    Call 01b_metrics.py in split-FASTA mode.

    Returns True if the script exited cleanly OR if it exited with a non-zero
    code only because some plugins failed (partial success).  The prediction
    script will catch any missing TSVs.  Returns False only when the metrics
    script crashes before producing any output at all.
    """
    cmd = [
        sys.executable, str(metrics_script),
        '--fasta-5utr', str(fasta_5utr),
        '--fasta-cds',  str(fasta_cds),
        '--fasta-3utr', str(fasta_3utr),
        '--species',    species,
        '--output-dir', str(output_dir),
    ]
    if force:
        cmd.append('--force')

    log.debug(f"metrics cmd: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    # Always surface the full stderr so failing plugin names are visible.
    if result.stderr:
        # Print each line with a prefix so it's easy to spot in terminal output.
        for line in result.stderr.splitlines():
            # Suppress blank lines and duplicate timestamp noise at DEBUG level.
            if line.strip():
                if '[ERROR]' in line or 'failed' in line.lower():
                    log.warning(f"[metrics] {line}")
                else:
                    log.debug(f"[metrics] {line}")

    if result.returncode != 0:
        # Check whether the metrics dir has *any* TSV output — if so, treat
        # this as a partial success (some plugins failed, rest succeeded).
        metrics_tsv_dir = output_dir / 'metrics'
        tsv_count = len(list(metrics_tsv_dir.glob('*.tsv'))) if metrics_tsv_dir.exists() else 0
        if tsv_count > 0:
            log.warning(
                f"metrics script exited with code {result.returncode} "
                f"but produced {tsv_count} TSV(s) — continuing (partial success). "
                f"Check [metrics] WARNING lines above for which plugins failed."
            )
            return True   # let the prediction script decide if it has enough
        else:
            log.error(
                f"metrics script failed (exit {result.returncode}) with no TSV output.\n"
                f"Full stderr:\n{result.stderr}"
            )
            return False
    return True


def run_prediction(predict_script: Path, metrics_dir: Path) -> Optional[pd.DataFrame]:
    """
    Call the LightGBM prediction script on the metrics TSV directory.
    Returns a DataFrame with columns [transcript_id, predicted_halflife]
    or None on failure.
    """
    cmd = [sys.executable, str(predict_script), '--input', str(metrics_dir)]
    log.debug(f"predict cmd: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(
            f"prediction script failed (exit {result.returncode}).\n"
            f"── stdout ──\n{result.stdout[-3000:]}\n"
            f"── stderr ──\n{result.stderr[-3000:]}"
        )
        return None

    out_path = metrics_dir / 'predictions.tsv'
    if not out_path.exists():
        log.error(
            f"Prediction output not found: {out_path}\n"
            f"── stdout ──\n{result.stdout[-2000:]}\n"
            f"── stderr ──\n{result.stderr[-2000:]}"
        )
        return None

    return pd.read_csv(out_path, sep='\t')


# ══════════════════════════════════════════════════════════════════════════════
# Fallback TSV synthesis
# Generates default-filled TSVs for plugins that cannot run without a GFF
# (nmd_fragility_*, junctions, architecture).  Triggered when a TSV is absent
# OR exists but contains only a header row (0 data rows) — the NMD plugins
# write an empty-but-headed file when no junction data is available.
# ══════════════════════════════════════════════════════════════════════════════

def _tsv_row_count(path: Path) -> int:
    """Return the number of data rows in a TSV (0 if missing or header-only)."""
    if not path.exists():
        return 0
    with open(path) as fh:
        lines = [l for l in fh if l.strip()]
    return max(0, len(lines) - 1)   # subtract header


# ── NMD shared columns (present in all three NMD variants) ───────────────────
_NMD_SHARED_COLS = [
    'transcript_id', 'gene_id', 'strand', 'model',
    'cds_length', 'zone_length',
    'n_transition_fragile_codons', 'n_transversion_fragile_codons',
    'n_snv_fragile_codons', 'n_alt_stop_codons',
]
_NMD_DENSITY_COLS = [
    'transition_fragile_codon_density',
    'transversion_fragile_codon_density',
    'snv_fragile_codon_density',
    'alt_stop_codon_density',
    'transition_fraction_of_snv_fragile',
]


def _nmd_base_row(transcript_id: str, cds_len: int, model_label: str) -> dict:
    """Shared NMD fields with zero counts."""
    return {
        'transcript_id':                transcript_id,
        'gene_id':                      transcript_id,
        'strand':                       '+',
        'model':                        model_label,
        'cds_length':                   cds_len,
        'zone_length':                  cds_len,
        'n_transition_fragile_codons':  0,
        'n_transversion_fragile_codons':0,
        'n_snv_fragile_codons':         0,
        'n_alt_stop_codons':            0,
    }


def _nmd_core_row(transcript_id: str, gene_id: str, cds_len: int) -> dict:
    """
    nmd_fragility_core.tsv — shared NMD counts only, no density columns.
    13 cols as seen in the prediction script output.
    """
    row = _nmd_base_row(transcript_id, cds_len, 'core')
    return row


def _nmd_full_row(transcript_id: str, gene_id: str, cds_len: int) -> dict:
    """
    nmd_fragility_full.tsv — shared counts + all five density/fraction cols.
    Fraction is NaN (0/0 undefined; LightGBM handles NaN natively).
    """
    row = _nmd_base_row(transcript_id, cds_len, 'full')
    row.update({
        'transition_fragile_codon_density':    0.0,
        'transversion_fragile_codon_density':  0.0,
        'snv_fragile_codon_density':           0.0,
        'alt_stop_codon_density':              0.0,
        'transition_fraction_of_snv_fragile':  '',   # blank → NaN on read
    })
    return row


def _nmd_window_row(transcript_id: str, gene_id: str, cds_len: int) -> dict:
    """
    nmd_fragility_window.tsv — same as full but model label = 'window'.
    The predict script drops density cols as duplicates of full, so only
    one unique column (e.g. a window-specific metric) remains.  We emit
    the full shared + density set so nothing is accidentally missing.
    """
    row = _nmd_base_row(transcript_id, cds_len, 'window')
    row.update({
        'transition_fragile_codon_density':    0.0,
        'transversion_fragile_codon_density':  0.0,
        'snv_fragile_codon_density':           0.0,
        'alt_stop_codon_density':              0.0,
        'transition_fraction_of_snv_fragile':  '',
    })
    return row


def _junctions_default_row(transcript_id: str, gene_id: str, cds_len: int) -> dict:
    """No introns → all junction counts 0, all distances NaN."""
    return {
        'transcript_id':                  transcript_id,
        'gene_id':                        transcript_id,
        'n_exons':                        1,
        'n_junctions':                    0,
        'strand':                         '+',
        'n_5utr_junctions':               0,
        'n_cds_junctions':                0,
        'n_3utr_junctions':               0,
        'stop_dist_closest_upstream':     '',
        'stop_dist_closest_downstream':   '',
        'stop_dist_last_downstream':      '',
        'start_dist_closest_upstream':    '',
        'start_dist_closest_downstream':  '',
    }


def _architecture_default_row(transcript_id: str, gene_id: str, cds_len: int) -> dict:
    """Single-exon transcript; intron stats all NaN."""
    return {
        'transcript_id':      transcript_id,
        'gene_id':            transcript_id,
        'n_exons':            1,
        'strand':             '+',
        'first_exon_length':  '',
        'last_exon_length':   '',
        'intron_mean':        '',
        'intron_median':      '',
        'intron_sd':          '',
    }


def _make_schema(factory) -> tuple:
    sample = factory('__sample__', '__sample__', 100)
    return (list(sample.keys()), factory)


# Maps TSV filename → (columns, row_factory).
# Only nmd_fragility_full is used by the prediction model; core and window
# are not needed and are intentionally excluded so the prediction script's
# inner-join does not require them.
_FALLBACK_SCHEMAS: dict = {
    'nmd_fragility_full.tsv': _make_schema(_nmd_full_row),
    'junctions.tsv':          _make_schema(_junctions_default_row),
    'architecture.tsv':       _make_schema(_architecture_default_row),
}


def _synthesise_fallback_tsvs(
    metrics_tsv_dir: Path,
    utr_ids: List[str],
    cds_len: int,
) -> None:
    """
    For each plugin TSV that is absent OR has 0 data rows (header-only),
    overwrite it with default-filled rows so the prediction inner-join succeeds.

    Header-only detection is important: NMD plugins write empty-headed files
    rather than no file at all when junction data is unavailable.
    """
    for tsv_name, (columns, row_factory) in _FALLBACK_SCHEMAS.items():
        tsv_path = metrics_tsv_dir / tsv_name
        n_rows = _tsv_row_count(tsv_path)
        if n_rows > 0:
            continue   # plugin produced real data — leave it alone

        action = "header-only (0 data rows)" if tsv_path.exists() else "missing"
        log.warning(
            f"[fallback] '{tsv_name}' is {action} — "
            f"writing default (zero/NaN) rows for {len(utr_ids)} transcripts."
        )
        metrics_tsv_dir.mkdir(parents=True, exist_ok=True)
        rows = [row_factory(sid, sid, cds_len) for sid in utr_ids]
        with open(tsv_path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=columns, delimiter='\t',
                                    lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        log.info(f"[fallback] Wrote {len(rows)} rows → {tsv_path}")


def evaluate_population(
    population: List[str],
    utr_ids: List[str],
    fixed_5utr_id: str,
    fixed_5utr_seq: str,
    fixed_cds_id: str,
    fixed_cds_seq: str,
    species: str,
    metrics_script: Path,
    predict_script: Path,
    work_dir: Path,
) -> Optional[List[float]]:
    """
    Write population sequences as temporary FASTAs, run the pipeline,
    return a list of predicted half-lives aligned with *population*.
    Returns None on pipeline failure.
    """
    # Write 5UTR FASTA (single fixed sequence — all samples share it)
    fasta_5 = work_dir / 'pop_5utr.fa'
    fasta_c = work_dir / 'pop_cds.fa'
    fasta_3 = work_dir / 'pop_3utr.fa'

    # 5UTR: replicate the fixed sequence once for every sample ID
    write_fasta(fasta_5, [(sid, fixed_5utr_seq) for sid in utr_ids])
    # CDS: single shared sequence (01b_metrics expects exactly 1)
    write_fasta(fasta_c, [(fixed_cds_id, fixed_cds_seq)])
    # 3UTR: one record per individual
    write_fasta(fasta_3, list(zip(utr_ids, population)))

    metrics_out = work_dir / 'metrics_run'
    metrics_tsv_dir = metrics_out / 'metrics'

    # Wipe the metrics TSV dir before each generation so that predictions.tsv
    # and stale TSVs from the previous generation don't pollute the inner-join.
    if metrics_tsv_dir.exists():
        import shutil as _shutil
        _shutil.rmtree(metrics_tsv_dir)
    metrics_out.mkdir(parents=True, exist_ok=True)

    ok = run_metrics(
        metrics_script, fasta_5, fasta_c, fasta_3,
        species, metrics_out, force=True,
    )
    if not ok:
        return None

    # Remove TSVs that the prediction model doesn't need but that the metrics
    # script may have written as empty/header-only files.  Leaving them in
    # causes the prediction script's inner-join to collapse to 0 rows.
    _UNUSED_TSVS = ('nmd_fragility_core.tsv', 'nmd_fragility_window.tsv')
    for _unused in _UNUSED_TSVS:
        _p = metrics_tsv_dir / _unused
        if _p.exists():
            _p.unlink()
            log.debug(f"[cleanup] Removed unused TSV: {_p.name}")

    # Synthesise default-filled TSVs for any GFF-dependent plugin that failed.
    # This is expected in split-FASTA / no-GFF mode for NMD, junctions, etc.
    _synthesise_fallback_tsvs(metrics_tsv_dir, utr_ids, len(fixed_cds_seq))

    preds = run_prediction(predict_script, metrics_tsv_dir)
    if preds is None:
        return None

    # Align predictions back to population order
    pred_map = dict(zip(preds['transcript_id'], preds['predicted_halflife']))
    scores = [pred_map.get(sid, float('-inf')) for sid in utr_ids]
    return scores


# ══════════════════════════════════════════════════════════════════════════════
# Plotting
# ══════════════════════════════════════════════════════════════════════════════

def plot_fitness_over_generations(best_per_gen: List[dict], out_path: Path) -> None:
    """
    Plot best / mean / worst predicted half-life per generation and save as PNG.
    Uses matplotlib with a non-interactive backend so it works headlessly.
    Failure to plot (e.g. matplotlib not installed) is logged as a warning,
    never a fatal error — the GA's TSV/FASTA outputs are unaffected.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        log.warning(
            "matplotlib is not installed — skipping fitness plot. "
            "Install it with `pip install matplotlib` to enable plotting."
        )
        return

    try:
        gens   = [row['generation'] for row in best_per_gen]
        best   = [row['best_score'] for row in best_per_gen]
        mean   = [row['mean_score'] for row in best_per_gen]
        worst  = [row['worst_score'] for row in best_per_gen]

        fig, ax = plt.subplots(figsize=(9, 5.5), dpi=150)

        ax.plot(gens, best, label='Best', color='#1b9e3e', linewidth=2, marker='o', markersize=3)
        ax.plot(gens, mean, label='Mean', color='#2166ac', linewidth=1.5, marker='o', markersize=3)
        ax.plot(gens, worst, label='Worst', color='#b2182b', linewidth=1.5, linestyle='--', alpha=0.7)

        # Shade the gap between best and worst each generation for visual context.
        ax.fill_between(gens, worst, best, color='#1b9e3e', alpha=0.07)

        # Mark the all-time best point.
        best_idx = max(range(len(best)), key=lambda i: best[i])
        ax.scatter([gens[best_idx]], [best[best_idx]], color='gold', edgecolor='black',
                   zorder=5, s=80, marker='*', label=f"All-time best (gen {gens[best_idx]})")

        ax.set_xlabel('Generation')
        ax.set_ylabel('Predicted half-life')
        ax.set_title("3' UTR GA — Predicted Half-Life Over Generations")
        ax.legend(loc='best', framealpha=0.9)
        ax.grid(True, alpha=0.3)
        if len(gens) <= 30:
            ax.set_xticks(gens)
        fig.tight_layout()
        fig.savefig(out_path)
        plt.close(fig)
        log.info(f"Wrote {out_path}")
    except Exception as exc:
        log.warning(f"Failed to generate fitness plot ({exc}); continuing without it.")


# ══════════════════════════════════════════════════════════════════════════════
# Genetic Algorithm
# ══════════════════════════════════════════════════════════════════════════════

def run_ga(
    seed_3utr_seq: str,
    fixed_5utr_id: str,
    fixed_5utr_seq: str,
    fixed_cds_id: str,
    fixed_cds_seq: str,
    species: str,
    metrics_script: Path,
    predict_script: Path,
    output_dir: Path,
    population_size: int,
    generations: int,
    mutation_rate: float,
    crossover_rate: float,
    elite_frac: float,
    tournament_k: int,
    rng_seed: int,
    make_plot: bool = True,
) -> None:
    rng = random.Random(rng_seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    results_dir = output_dir / 'results'
    results_dir.mkdir(exist_ok=True)

    n_elite = max(1, int(population_size * elite_frac))

    log.info("═" * 60)
    log.info("3' UTR Genetic Algorithm")
    log.info(f"  Population   : {population_size}")
    log.info(f"  Generations  : {generations}")
    log.info(f"  Mutation rate: {mutation_rate:.4f}")
    log.info(f"  Crossover    : {crossover_rate:.4f}")
    log.info(f"  Elites       : {n_elite}")
    log.info(f"  Seed 3UTR    : {len(seed_3utr_seq)} nt")
    log.info(f"  Output       : {output_dir}")
    log.info("═" * 60)

    population = [seed_3utr_seq] + [
        point_mutate(seed_3utr_seq, mutation_rate * 5, rng)
        for _ in range(population_size - 1)
    ]

    all_rows: List[dict] = []
    best_per_gen: List[dict] = []
    all_time_best_score = float('-inf')
    all_time_best_seq   = seed_3utr_seq

    work_dir = output_dir / '_workspace'
    work_dir.mkdir(exist_ok=True)

    for gen in range(1, generations + 1):
        gen_start = time.time()
        log.info(f"\n── Generation {gen}/{generations} ──")

        # Unique sample IDs so the merge in predict_script works
        utr_ids = [f"ind_{gen:04d}_{i:05d}" for i in range(population_size)]

        scores = evaluate_population(
            population, utr_ids,
            fixed_5utr_id, fixed_5utr_seq,
            fixed_cds_id, fixed_cds_seq,
            species, metrics_script, predict_script, work_dir,
        )

        if scores is None:
            log.error(f"Pipeline failed on generation {gen}. Aborting.")
            if best_per_gen and make_plot:
                plot_fitness_over_generations(
                    best_per_gen, results_dir / 'fitness_over_generations.png'
                )
            sys.exit(1)

        # ── Record results ────────────────────────────────────────────────
        for i, (sid, seq, score) in enumerate(zip(utr_ids, population, scores)):
            all_rows.append({
                'generation': gen,
                'rank_in_gen': 0,
                'sample_id': sid,
                'predicted_halflife': score,
                'sequence': seq,
                'length': len(seq),
                'mutations_from_seed': sum(a != b for a, b in zip(seq, seed_3utr_seq)
                                           if len(seq) == len(seed_3utr_seq)),
            })

        # Sort by score descending
        ranked = sorted(zip(scores, population, utr_ids),
                        key=lambda x: x[0], reverse=True)
        for rank, (sc, sq, sid) in enumerate(ranked):
            # update rank in all_rows
            for row in all_rows:
                if row['sample_id'] == sid:
                    row['rank_in_gen'] = rank + 1

        best_score, best_seq, best_id = ranked[0]
        gen_elapsed = time.time() - gen_start
        log.info(f"  Best this gen : {best_score:.4f}  ({best_id})")
        log.info(f"  Worst this gen: {ranked[-1][0]:.4f}")
        log.info(f"  Mean          : {sum(scores)/len(scores):.4f}")
        log.info(f"  Elapsed       : {gen_elapsed:.1f}s")

        best_per_gen.append({
            'generation': gen,
            'best_score': best_score,
            'mean_score': sum(scores) / len(scores),
            'worst_score': ranked[-1][0],
            'best_id': best_id,
            'best_sequence': best_seq,
        })

        if best_score > all_time_best_score:
            all_time_best_score = best_score
            all_time_best_seq   = best_seq
            log.info(f"  ★ New all-time best: {all_time_best_score:.4f}")

        if gen == generations:
            break

        # ── Breed next generation ─────────────────────────────────────────
        elite_seqs = [sq for _, sq, _ in ranked[:n_elite]]
        non_elite_scores = [sc for sc, _, _ in ranked]
        non_elite_seqs   = [sq for _, sq, _ in ranked]

        next_gen = list(elite_seqs)

        while len(next_gen) < population_size:
            parent_a = tournament_select(non_elite_seqs, non_elite_scores,
                                         tournament_k, rng)
            if rng.random() < crossover_rate:
                parent_b = tournament_select(non_elite_seqs, non_elite_scores,
                                              tournament_k, rng)
                child = uniform_crossover(parent_a, parent_b, rng)
            else:
                child = parent_a

            child = point_mutate(child, mutation_rate, rng)
            next_gen.append(child)

        population = next_gen

    # ══ Write outputs ══════════════════════════════════════════════════════
    # best_per_generation.tsv
    bpg_path = results_dir / 'best_per_generation.tsv'
    with open(bpg_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(best_per_gen[0].keys()),
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(best_per_gen)
    log.info(f"Wrote {bpg_path}")

    # all_generations.tsv
    all_path = results_dir / 'all_generations.tsv'
    with open(all_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()),
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(all_rows)
    log.info(f"Wrote {all_path}")

    # best_3utr.fa
    best_fa = results_dir / 'best_3utr.fa'
    write_fasta(best_fa, [('best_3utr_evolved', all_time_best_seq)])
    log.info(f"Wrote {best_fa}")

    # seed_3utr.fa  (for reference)
    seed_fa = results_dir / 'seed_3utr.fa'
    write_fasta(seed_fa, [('seed_3utr_original', seed_3utr_seq)])

    # fitness_over_generations.png
    if make_plot:
        plot_fitness_over_generations(
            best_per_gen, results_dir / 'fitness_over_generations.png'
        )

    summary_path = results_dir / 'evolution_summary.txt'
    with open(summary_path, 'w') as fh:
        fh.write("3' UTR Genetic Algorithm — Evolution Summary\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(f"Generations     : {generations}\n")
        fh.write(f"Population size : {population_size}\n")
        fh.write(f"Mutation rate   : {mutation_rate}\n")
        fh.write(f"Crossover rate  : {crossover_rate}\n")
        fh.write(f"Elite fraction  : {elite_frac}  ({n_elite} individuals)\n")
        fh.write(f"Tournament k    : {tournament_k}\n")
        fh.write(f"RNG seed        : {rng_seed}\n\n")
        fh.write(f"Seed 3UTR length : {len(seed_3utr_seq)} nt\n")
        fh.write(f"Best 3UTR length : {len(all_time_best_seq)} nt\n\n")
        fh.write(f"Seed predicted half-life  : {best_per_gen[0]['best_score']:.4f}  "
                 f"(generation 1 best)\n")
        fh.write(f"Final best half-life      : {all_time_best_score:.4f}\n\n")
        mutations = sum(a != b for a, b in zip(seed_3utr_seq, all_time_best_seq)
                        if len(seed_3utr_seq) == len(all_time_best_seq))
        fh.write(f"Point mutations from seed : {mutations} / {len(seed_3utr_seq)}\n\n")
        fh.write("Generation-by-generation best scores:\n")
        for row in best_per_gen:
            fh.write(f"  Gen {row['generation']:>3d} : {row['best_score']:.4f}  "
                     f"(mean {row['mean_score']:.4f}, worst {row['worst_score']:.4f})\n")
        fh.write("\nBest evolved 3UTR sequence:\n")
        for i in range(0, len(all_time_best_seq), 60):
            fh.write(all_time_best_seq[i:i + 60] + '\n')

    log.info(f"Wrote {summary_path}")
    log.info("\n" + "═" * 60)
    log.info(f"GA complete.  All-time best predicted half-life: {all_time_best_score:.4f}")
    log.info(f"Results in: {results_dir}")
    log.info("═" * 60)


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── Required inputs ───────────────────────────────────────────────────
    parser.add_argument('--fasta-5utr', required=True, dest='fasta_5utr',
                        metavar='FILE',
                        help="FASTA with a single 5' UTR sequence (fixed, not mutated).")
    parser.add_argument('--fasta-cds', required=True, dest='fasta_cds',
                        metavar='FILE',
                        help="FASTA with a single CDS sequence (fixed, not mutated).")
    parser.add_argument('--fasta-3utr', required=True, dest='fasta_3utr',
                        metavar='FILE',
                        help="FASTA with the seed 3' UTR sequence to evolve.")
    parser.add_argument('--species', required=True, metavar='BINOMIAL',
                        help="Species name passed to the metrics script (e.g. 'Homo sapiens').")

    # ── Script paths ──────────────────────────────────────────────────────
    parser.add_argument('--metrics-script', required=True, dest='metrics_script',
                        metavar='FILE',
                        help="Path to 01b_metrics.py.")
    parser.add_argument('--predict-script', required=True, dest='predict_script',
                        metavar='FILE',
                        help="Path to the LightGBM prediction script.")

    # ── GA parameters ─────────────────────────────────────────────────────
    parser.add_argument('--population', type=int, default=50, metavar='N',
                        help="Population size per generation (default: 50).")
    parser.add_argument('--generations', type=int, default=30, metavar='N',
                        help="Number of GA generations (default: 30).")
    parser.add_argument('--mutation-rate', type=float, default=0.02,
                        dest='mutation_rate', metavar='RATE',
                        help="Per-nucleotide point mutation probability (default: 0.02).")
    parser.add_argument('--crossover-rate', type=float, default=0.5,
                        dest='crossover_rate', metavar='RATE',
                        help="Probability of crossover vs. cloning (default: 0.5).")
    parser.add_argument('--elite-frac', type=float, default=0.1,
                        dest='elite_frac', metavar='FRAC',
                        help="Fraction of top individuals preserved unchanged (default: 0.1).")
    parser.add_argument('--tournament-k', type=int, default=3,
                        dest='tournament_k', metavar='K',
                        help="Tournament selection size (default: 3).")
    parser.add_argument('--seed', type=int, default=42, metavar='INT',
                        help="Random seed for reproducibility (default: 42).")

    # ── Output ────────────────────────────────────────────────────────────
    parser.add_argument('--output-dir', '-o', default='ga_output',
                        dest='output_dir', metavar='DIR',
                        help="Output directory (default: ga_output/).")
    parser.add_argument('--no-plot', action='store_true', dest='no_plot',
                        help="Skip generating fitness_over_generations.png.")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Enable DEBUG logging.")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # ── Validate paths ────────────────────────────────────────────────────
    for label, path_str in [
        ('--fasta-5utr', args.fasta_5utr),
        ('--fasta-cds',  args.fasta_cds),
        ('--fasta-3utr', args.fasta_3utr),
        ('--metrics-script', args.metrics_script),
        ('--predict-script', args.predict_script),
    ]:
        p = Path(path_str)
        if not p.exists():
            log.error(f"{label} not found: {p}")
            sys.exit(1)

    # ── Load fixed sequences ──────────────────────────────────────────────
    try:
        u5_id,  u5_seq  = read_single_fasta(Path(args.fasta_5utr))
        cds_id, cds_seq = read_single_fasta(Path(args.fasta_cds))
        u3_id,  u3_seq  = read_single_fasta(Path(args.fasta_3utr))
    except ValueError as exc:
        log.error(str(exc))
        sys.exit(1)

    log.info(f"5' UTR  : '{u5_id}'   {len(u5_seq)} nt  (fixed)")
    log.info(f"CDS     : '{cds_id}'  {len(cds_seq)} nt  (fixed)")
    log.info(f"3' UTR  : '{u3_id}'   {len(u3_seq)} nt  (seed — will be evolved)")

    # ── Run GA ────────────────────────────────────────────────────────────
    run_ga(
        seed_3utr_seq   = u3_seq,
        fixed_5utr_id   = u5_id,
        fixed_5utr_seq  = u5_seq,
        fixed_cds_id    = cds_id,
        fixed_cds_seq   = cds_seq,
        species         = args.species,
        metrics_script  = Path(args.metrics_script),
        predict_script  = Path(args.predict_script),
        output_dir      = Path(args.output_dir),
        population_size = args.population,
        generations     = args.generations,
        mutation_rate   = args.mutation_rate,
        crossover_rate  = args.crossover_rate,
        elite_frac      = args.elite_frac,
        tournament_k    = args.tournament_k,
        rng_seed        = args.seed,
        make_plot       = not args.no_plot,
    )


if __name__ == '__main__':
    main()