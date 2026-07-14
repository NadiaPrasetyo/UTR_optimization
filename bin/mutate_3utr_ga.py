#!/usr/bin/env python3
"""
mutate_3utr_ga.py
─────────────────────────────────────────────────────────────────────────────
Point-scored genetic algorithm for evolving a 3' UTR sequence, selecting
purely on:

  1. AlphaFold3-predicted structure RMSD vs. 3 representative reference
     structures (minimum over the 3) — DECREASED relative to the parent
     that produced this individual  → +1 point
  2. Infernal cmsearch bit score against a covariance model of the
     original structure — INCREASED relative to the parent → +1 point;
     NO HIT at all → --no-hit-penalty (default -10), overriding the +1
  3. Percent identity (PID) to a set of patent sequences (Smith-Waterman) —
     PID > --patent-pid-threshold (default 80%) to ANY patent →
     --patent-pid-penalty (default -10)
  4. LightGBM-predicted half-life is computed and reported every
     generation but carries NO weight in the score — it is not part of
     selection at all.

Per-individual integer scores are converted into a sampling probability
distribution (softmax, temperature-controlled), and --n-select individuals
(default 10) are drawn WITH REPLACEMENT from the full population (default
1,000) using that distribution. Those (possibly repeated) individuals seed
the next generation via point mutation only (uniform crossover is
available but OFF by default — see --crossover-rate).

There is no mutation-rate annealing/boosting of any kind: a single flat
--mutation-rate (per-nucleotide substitution probability) is used for
every generation.

AlphaFold3 structure prediction for the whole population is parallelized
on SLURM as a single array job per generation (reusing the exact JSON /
job-list conventions of prepare_af3_jobs.py), and RMSD is computed with
the Needleman-Wunsch (global-alignment) Biopython RMSD method only (no
PyMOL) from compare_af3_candidates.py. x3dna-dssr dot-bracket secondary
structure is additionally captured per candidate for reporting (it does
not affect scoring).

USAGE
─────
python mutate_3utr_ga.py \
    --fasta-5utr   input_5utr.fa \
    --fasta-cds    input_cds.fa \
    --fasta-3utr   seed_3utr.fa \
    --species      "Homo sapiens" \
    --metrics-script  bin/01b_metrics.py \
    --predict-script  predict_halflife.py \
    --ref-cif ref1.cif --ref-cif ref2.cif --ref-cif ref3.cif \
    --cm-model     original_3utr.cm \
    --af3-work-dir /scratch/me/af3_ga \
    --af3-models   /projects/.../alphafold3-weights \
    --output-dir   ga_output \
    [--population  1000] \
    [--n-select    10] \
    [--generations 30] \
    [--mutation-rate 0.02] \
    [--crossover-rate 0.0] \
    [--temperature 1.0] \
    [--cm-evalue-threshold 0.01] \
    [--no-hit-penalty -10] \
    [--patent-pid-threshold 80.0] \
    [--patent-pid-penalty -10] \
    [--seed 42] \
    [--no-plot]

REQUIRED COMPANION SCRIPTS
───────────────────────────
prepare_af3_jobs.py and compare_af3_candidates.py must sit next to this
script (they are imported directly as libraries for AF3 JSON/job-list
construction and for the Biopython RMSD + DSSR helpers).

WHAT IT DOES EACH GENERATION
──────────────────────────────
  a. Submits one SLURM array job covering the whole population's AF3
     monomer predictions (one task per individual), waits for it to
     finish, and collects each individual's top model .cif.
  b. Computes RMSD (Needleman-Wunsch Biopython RMSD, min over the 3
     reference .cif files) for every individual, and (optionally) runs
     x3dna-dssr on each candidate structure to log its dot-bracket
     string.
  c. Runs cmsearch for the whole population against --cm-model.
  d. Computes Smith-Waterman PID against the patent sequence set.
  e. Runs 01b_metrics.py + the LightGBM predictor for predicted
     half-life (reporting only).
  f. Scores each individual by comparing (1)-(3) to the metrics of the
     specific parent that produced it (for generation 1, the parent is
     the seed sequence).
  g. Converts scores → sampling probabilities (softmax) and draws
     --n-select individuals WITH REPLACEMENT.
  h. Breeds the next --population individuals from those n-select
     parents via point mutation (± optional crossover).
─────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import logging
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# ── Logging ──────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('ga_3utr')

NUCLEOTIDES = ['A', 'T', 'G', 'C']

# ── Patent sequences (from US/international patent applications) ─────────
PATENT_SEQUENCES = [
    {"id": "A7_30nt", "start": 9651, "end": 9680, "seq": "CAGACCCTGGTCCGGGGCAATGGGACCACT"},
    {"id": "A7_32nt", "start": 9650, "end": 9681, "seq": "TCAGACCCTGGTCCGGGGCAAATGGGACCACTG"},
    {"id": "A7_34nt", "start": 9649, "end": 9682, "seq": "GTCAGACCCTGGTCCGGGGCAATGGGACCACTGT"},
    {"id": "A7_36nt", "start": 9648, "end": 9683, "seq": "GGTCAGACCCTGGTCCGGGGCAAATGGGACCACTGTT"},
    {"id": "A7_40nt", "start": 9646, "end": 9685, "seq": "TGGGTCAGACCCTGGTCCGGGGCAAATGGGACCACTGTTTC"},
    {"id": "A7_43nt", "start": 9646, "end": 9688, "seq": "TGGGTCAGACCCTGGTCCGGGGCAATGGGACCACTGTTTCGCG"},
    {"id": "A7_47nt", "start": 9646, "end": 9692, "seq": "TGGGTCAGACCCTGGTCCGGGGCAATGGGACCACTGTTTCGCGTTTA"},
]


# ══════════════════════════════════════════════════════════════════════════
# Sibling-script imports (prepare_af3_jobs.py, compare_af3_candidates.py)
# ══════════════════════════════════════════════════════════════════════════

def _import_sibling(module_name: str, filename: str):
    here = Path(__file__).resolve().parent
    path = here / filename
    if not path.exists():
        raise FileNotFoundError(
            f"Required helper script '{filename}' was not found next to "
            f"{Path(__file__).name} (looked in {here}). Place "
            f"prepare_af3_jobs.py and compare_af3_candidates.py in the "
            f"same directory as this script."
        )
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


af3prep = None  # populated in main() once we know argv[0]'s directory
af3cmp = None


# ══════════════════════════════════════════════════════════════════════════
# FASTA I/O
# ══════════════════════════════════════════════════════════════════════════

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


# ══════════════════════════════════════════════════════════════════════════
# Genetic operators — flat mutation rate only, no annealing/boosting
# ══════════════════════════════════════════════════════════════════════════

def point_mutate(seq: str, mutation_rate: float, rng: random.Random) -> str:
    """Apply random point mutations at a constant *mutation_rate* per nt."""
    bases = list(seq)
    for i, base in enumerate(bases):
        if rng.random() < mutation_rate:
            alt = [n for n in NUCLEOTIDES if n != base]
            bases[i] = rng.choice(alt)
    return ''.join(bases)


def uniform_crossover(seq_a: str, seq_b: str, rng: random.Random) -> str:
    """Uniform crossover: each position drawn from either parent, p=0.5."""
    if len(seq_a) != len(seq_b):
        return seq_a
    return ''.join(a if rng.random() < 0.5 else b for a, b in zip(seq_a, seq_b))


# ══════════════════════════════════════════════════════════════════════════
# Patent sequence identity (PID) checking via Smith-Waterman alignment
# ══════════════════════════════════════════════════════════════════════════

def _sw_pure(seq_a: str, seq_b: str) -> dict:
    """Pure-Python Smith-Waterman local alignment fallback."""
    match_score, mismatch_score, gap_score = 2, -1, -1
    m, n = len(seq_a), len(seq_b)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    max_score, max_i, max_j = 0, 0, 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score = (H[i - 1][j - 1] + match_score) if seq_a[i - 1] == seq_b[j - 1] \
                else (H[i - 1][j - 1] + mismatch_score)
            score = max(score, H[i - 1][j] + gap_score, H[i][j - 1] + gap_score, 0)
            H[i][j] = score
            if score > max_score:
                max_score, max_i, max_j = score, i, j

    align_a, align_b = [], []
    i, j = max_i, max_j
    while i > 0 and j > 0 and H[i][j] > 0:
        if seq_a[i - 1] == seq_b[j - 1]:
            align_a.append(seq_a[i - 1]); align_b.append(seq_b[j - 1]); i -= 1; j -= 1
        elif H[i - 1][j - 1] >= H[i - 1][j] and H[i - 1][j - 1] >= H[i][j - 1]:
            align_a.append(seq_a[i - 1]); align_b.append(seq_b[j - 1]); i -= 1; j -= 1
        elif H[i - 1][j] > H[i][j - 1]:
            align_a.append(seq_a[i - 1]); align_b.append('-'); i -= 1
        else:
            align_a.append('-'); align_b.append(seq_b[j - 1]); j -= 1

    align_a = ''.join(reversed(align_a))
    align_b = ''.join(reversed(align_b))
    identities = sum(a == b and a != '-' for a, b in zip(align_a, align_b))
    return {'identities': identities, 'alignment_length': max(len(align_a), len(align_b))}


def _sw_biopython(seq_a: str, seq_b: str) -> Optional[dict]:
    try:
        from Bio import pairwise2
        from Bio.Seq import Seq
    except ImportError:
        return None

    alignments = pairwise2.align.localms(Seq(seq_a), Seq(seq_b), 2, -1, -1, -1)
    if not alignments:
        return None
    best = alignments[0]
    align_a, align_b = str(best[0]), str(best[1])
    identities = sum(a == b and a != '-' for a, b in zip(align_a, align_b))
    return {'identities': identities, 'alignment_length': max(len(align_a), len(align_b))}


def smith_waterman(seq_a: str, seq_b: str) -> dict:
    try:
        result = _sw_biopython(seq_a, seq_b)
        if result:
            return result
    except Exception:
        pass
    return _sw_pure(seq_a, seq_b)


def calculate_pid(seq_query: str, seq_subject: str) -> float:
    """PID = (identities / length_of_shorter_sequence) * 100."""
    alignment = smith_waterman(seq_query, seq_subject)
    shorter_len = min(len(seq_query), len(seq_subject))
    if shorter_len == 0:
        return 0.0
    return (alignment['identities'] / shorter_len) * 100.0


def check_patent_pid(sequence: str, patent_seqs: List[dict]) -> List[float]:
    """Return the list of PID values (%) aligned with patent_seqs."""
    return [calculate_pid(sequence, p['seq']) for p in patent_seqs]


# ══════════════════════════════════════════════════════════════════════════
# AlphaFold3: SLURM-parallelized structure prediction for a whole
# generation's population
# ══════════════════════════════════════════════════════════════════════════

def submit_af3_population(
    population: List[str],
    utr_ids: List[str],
    af3_work_dir: Path,
    af3_db: str,
    af3_models: str,
    af3_module: str,
    partition: str,
    time_limit: str,
    mem: str,
    cpus_per_task: int,
    gres: str,
    exclude: str,
    array_max_parallel: Optional[int],
    model_seeds: List[int],
    generation: int,
) -> Tuple[Path, Path, str]:
    """
    Build one AlphaFold3 JSON per individual (RNA monomer, chain A) using
    prepare_af3_jobs.py's own build_af3_json/sanitize_name/clean_RNA_seq,
    write a SLURM array script covering the whole population, and submit
    it with `sbatch`. Returns (jobs_tsv, output_root, slurm_job_id).
    """
    gen_dir = af3_work_dir / f"gen_{generation:04d}"
    jobs_dir = gen_dir / "af3_jobs"
    inputs_dir = jobs_dir / "inputs"
    logs_dir = jobs_dir / "logs"
    output_root = gen_dir / "af3_output"
    for d in (inputs_dir, logs_dir, output_root):
        d.mkdir(parents=True, exist_ok=True)

    jobs = []
    for sid, seq in zip(utr_ids, population):
        job_name = af3prep.sanitize_name(sid)
        chains = [{"id": "A", "sequence": af3prep.clean_RNA_seq(seq), "type": "rna"}]
        af3_json = af3prep.build_af3_json(job_name, chains, model_seeds)
        json_path = inputs_dir / f"{job_name}.json"
        json_path.write_text(json.dumps(af3_json, indent=2))
        jobs.append((job_name, json_path))

    jobs_tsv = jobs_dir / "jobs.tsv"
    with open(jobs_tsv, "w") as fh:
        for idx, (job_name, json_path) in enumerate(jobs):
            out_dir = output_root / job_name
            fh.write(f"{idx}\t{job_name}\t{json_path}\t{out_dir}\n")

    n_jobs = len(jobs)
    array_range = f"0-{n_jobs - 1}"
    if array_max_parallel:
        array_range += f"%{array_max_parallel}"
    exclude_line = f"#SBATCH --exclude={exclude}\n" if exclude else ""

    script = f"""#!/bin/bash
#SBATCH --job-name=af3_gen{generation:04d}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --mem={mem}
#SBATCH --time={time_limit}
#SBATCH --partition={partition}
#SBATCH --gres={gres}
#SBATCH --array={array_range}
#SBATCH --output={logs_dir}/af3_%A_%a.log
{exclude_line}
module purge
module load {af3_module}

export AF3_DB={af3_db}
export AF3_WD={gen_dir}
export AF3_MODELS={af3_models}

JOBS_TSV="{jobs_tsv}"
LINE=$(awk -F'\\t' -v idx="$SLURM_ARRAY_TASK_ID" '$1 == idx {{print}}' "$JOBS_TSV")

if [ -z "$LINE" ]; then
    echo "ERROR: no job found for array task ID $SLURM_ARRAY_TASK_ID in $JOBS_TSV"
    exit 1
fi

JOB_NAME=$(echo "$LINE" | cut -f2)
JSON_PATH=$(echo "$LINE" | cut -f3)
OUT_DIR=$(echo "$LINE" | cut -f4)
mkdir -p "$OUT_DIR"

echo "Running AlphaFold3 for job: $JOB_NAME"
af3 "$JSON_PATH" "$OUT_DIR" --run_data_pipeline=true --run_inference=true
echo "AlphaFold3 job '$JOB_NAME' finished at $(date)."
"""
    script_path = jobs_dir / "run_af3_array.sh"
    script_path.write_text(script)
    script_path.chmod(0o755)

    result = subprocess.run(["sbatch", "--parsable", str(script_path)],
                             capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"sbatch failed for generation {generation}: {result.stderr}")
    job_id = result.stdout.strip().split(';')[0]
    log.info(f"  Submitted AF3 array job {job_id} ({n_jobs} tasks) → {jobs_tsv}")
    return jobs_tsv, output_root, job_id


def wait_for_slurm_job(job_id: str, poll_interval: float = 30.0,
                        timeout: Optional[float] = None) -> None:
    """Poll `squeue` until *job_id* has no tasks left queued/running."""
    start = time.time()
    while True:
        result = subprocess.run(["squeue", "-j", job_id, "-h"],
                                 capture_output=True, text=True)
        # squeue returns nonzero once the job has fully left the queue on
        # some SLURM versions too, so treat "no rows" OR nonzero as done.
        if result.returncode != 0 or not result.stdout.strip():
            log.info(f"  SLURM array job {job_id} has left the queue.")
            return
        n_left = len(result.stdout.strip().splitlines())
        log.info(f"  Waiting on SLURM job {job_id}: {n_left} task(s) still queued/running "
                 f"(poll every {poll_interval:.0f}s)...")
        if timeout and (time.time() - start) > timeout:
            log.warning(f"  Timed out waiting for SLURM job {job_id} after {timeout}s "
                        f"— continuing with whatever output exists.")
            return
        time.sleep(poll_interval)


def collect_af3_cifs(output_root: Path, utr_ids: List[str],
                      model_glob: str) -> Dict[str, Optional[Path]]:
    """
    Locate each individual's top-ranked AF3 model .cif under
    output_root/<sanitized_job_name>/ using *model_glob* (relative,
    e.g. "**/*_model.cif"). AF3's exact output layout/naming can vary by
    version/site — adjust --af3-model-glob if nothing is found.
    """
    result: Dict[str, Optional[Path]] = {}
    for sid in utr_ids:
        job_name = af3prep.sanitize_name(sid)
        job_out_dir = output_root / job_name
        matches = sorted(job_out_dir.glob(model_glob))
        result[sid] = matches[0] if matches else None
        if not matches:
            log.warning(f"  No AF3 output structure found for {sid} in "
                        f"{job_out_dir} (glob '{model_glob}') — treated as a "
                        f"failed prediction for this generation.")
    return result


def evaluate_af3_rmsd(
    cif_paths: Dict[str, Optional[Path]],
    utr_ids: List[str],
    ref_cifs: List[Path],
    atom_name: str,
    run_dssr: bool,
    dssr_path: str,
    dssr_outdir: Path,
) -> Tuple[Dict[str, Optional[float]], Dict[str, Optional[str]]]:
    """
    For every individual's predicted structure, compute the Needleman-
    Wunsch (global-alignment) Biopython RMSD against each reference .cif
    and keep the minimum. If *run_dssr*, also fetch the candidate's
    x3dna-dssr dot-bracket string for reporting only (never affects
    scoring; failures here are logged and ignored).
    """
    rmsd_min: Dict[str, Optional[float]] = {}
    dbn_map: Dict[str, Optional[str]] = {}

    for sid in utr_ids:
        cif = cif_paths.get(sid)
        if cif is None or not Path(cif).exists():
            rmsd_min[sid] = None
            dbn_map[sid] = None
            continue

        best = None
        for ref in ref_cifs:
            try:
                summary, _sup, _s2 = af3cmp.biopython_rmsd(
                    str(ref), str(cif), None, None, atom_name)
                if best is None or summary.rmsd < best:
                    best = summary.rmsd
            except Exception as exc:
                log.debug(f"  RMSD failed for {sid} vs {ref}: {exc}")
        rmsd_min[sid] = best

        if run_dssr:
            try:
                dres = af3cmp.get_dssr_summary(str(cif), str(dssr_outdir), dssr_path)
                dbn_map[sid] = dres.dbn
            except Exception as exc:
                log.debug(f"  DSSR failed for {sid}: {exc}")
                dbn_map[sid] = None
        else:
            dbn_map[sid] = None

    n_ok = sum(1 for v in rmsd_min.values() if v is not None)
    log.info(f"  AF3/RMSD: {n_ok}/{len(utr_ids)} individuals have a usable "
             f"predicted structure and RMSD.")
    return rmsd_min, dbn_map


# ══════════════════════════════════════════════════════════════════════════
# Pipeline: metrics → predicted half-life (reported only, never weighted)
# ══════════════════════════════════════════════════════════════════════════

def run_metrics(metrics_script: Path, fasta_5utr: Path, fasta_cds: Path,
                 fasta_3utr: Path, species: str, output_dir: Path,
                 force: bool = True) -> bool:
    """Call 01b_metrics.py in split-FASTA mode."""
    cmd = [
        sys.executable, str(metrics_script),
        '--fasta-5utr', str(fasta_5utr), '--fasta-cds', str(fasta_cds),
        '--fasta-3utr', str(fasta_3utr), '--species', species,
        '--output-dir', str(output_dir),
    ]
    if force:
        cmd.append('--force')

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stderr:
        for line in result.stderr.splitlines():
            if line.strip():
                (log.warning if ('[ERROR]' in line or 'failed' in line.lower())
                 else log.debug)(f"[metrics] {line}")

    if result.returncode != 0:
        metrics_tsv_dir = output_dir / 'metrics'
        tsv_count = len(list(metrics_tsv_dir.glob('*.tsv'))) if metrics_tsv_dir.exists() else 0
        if tsv_count > 0:
            log.warning(f"metrics script exited {result.returncode} but produced "
                        f"{tsv_count} TSV(s) — continuing (partial success).")
            return True
        log.error(f"metrics script failed (exit {result.returncode}) with no TSV output.\n"
                  f"Full stderr:\n{result.stderr}")
        return False
    return True


def run_prediction(predict_script: Path, metrics_dir: Path) -> Optional[pd.DataFrame]:
    """Call the LightGBM prediction script; returns [transcript_id, predicted_halflife]."""
    cmd = [sys.executable, str(predict_script), '--input', str(metrics_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"prediction script failed (exit {result.returncode}).\n"
                  f"stdout:\n{result.stdout[-3000:]}\nstderr:\n{result.stderr[-3000:]}")
        return None
    out_path = metrics_dir / 'predictions.tsv'
    if not out_path.exists():
        log.error(f"Prediction output not found: {out_path}")
        return None
    return pd.read_csv(out_path, sep='\t')


_NMD_FULL_COLS = [
    'transcript_id', 'gene_id', 'strand', 'model', 'cds_length', 'zone_length',
    'n_transition_fragile_codons', 'n_transversion_fragile_codons',
    'n_snv_fragile_codons', 'n_alt_stop_codons',
    'transition_fragile_codon_density', 'transversion_fragile_codon_density',
    'snv_fragile_codon_density', 'alt_stop_codon_density',
    'transition_fraction_of_snv_fragile',
]
_JUNCTIONS_COLS = [
    'transcript_id', 'gene_id', 'n_exons', 'n_junctions', 'strand',
    'n_5utr_junctions', 'n_cds_junctions', 'n_3utr_junctions',
    'stop_dist_closest_upstream', 'stop_dist_closest_downstream',
    'stop_dist_last_downstream', 'start_dist_closest_upstream',
    'start_dist_closest_downstream',
]
_ARCHITECTURE_COLS = [
    'transcript_id', 'gene_id', 'strand', 'n_exons', 'n_internal_exons',
    'first_exon_length', 'last_exon_length', 'internal_exon_mean',
    'internal_exon_median', 'internal_exon_sd', 'intron_mean',
    'intron_median', 'intron_sd',
]


def _nmd_full_row(sid: str, cds_len: int) -> dict:
    return {
        'transcript_id': sid, 'gene_id': sid, 'strand': '+', 'model': 'full',
        'cds_length': cds_len, 'zone_length': cds_len,
        'n_transition_fragile_codons': 0, 'n_transversion_fragile_codons': 0,
        'n_snv_fragile_codons': 0, 'n_alt_stop_codons': 0,
        'transition_fragile_codon_density': 0.0, 'transversion_fragile_codon_density': 0.0,
        'snv_fragile_codon_density': 0.0, 'alt_stop_codon_density': 0.0,
        'transition_fraction_of_snv_fragile': '',
    }


def _junctions_row(sid: str) -> dict:
    return {
        'transcript_id': sid, 'gene_id': sid, 'n_exons': 1, 'n_junctions': 0,
        'strand': '+', 'n_5utr_junctions': 0, 'n_cds_junctions': 0,
        'n_3utr_junctions': 0, 'stop_dist_closest_upstream': '',
        'stop_dist_closest_downstream': '', 'stop_dist_last_downstream': '',
        'start_dist_closest_upstream': '', 'start_dist_closest_downstream': '',
    }


def _architecture_row(sid: str) -> dict:
    return {
        'transcript_id': sid, 'gene_id': sid, 'strand': '+', 'n_exons': 1,
        'n_internal_exons': 0, 'first_exon_length': '', 'last_exon_length': '',
        'internal_exon_mean': '', 'internal_exon_median': '', 'internal_exon_sd': '',
        'intron_mean': '', 'intron_median': '', 'intron_sd': '',
    }


def _tsv_row_count(path: Path) -> int:
    if not path.exists():
        return 0
    with open(path) as fh:
        lines = [l for l in fh if l.strip()]
    return max(0, len(lines) - 1)


def _synthesise_fallback_tsvs(metrics_tsv_dir: Path, utr_ids: List[str], cds_len: int) -> None:
    """
    Fill in default (zero/NaN) rows for GFF-dependent plugin TSVs that are
    absent or header-only, so the prediction script's inner-join succeeds
    in split-FASTA / no-GFF mode.
    """
    schemas = {
        'nmd_fragility_full.tsv': (_NMD_FULL_COLS, lambda sid: _nmd_full_row(sid, cds_len)),
        'junctions.tsv':          (_JUNCTIONS_COLS, _junctions_row),
        'architecture.tsv':       (_ARCHITECTURE_COLS, _architecture_row),
    }
    for tsv_name, (columns, row_factory) in schemas.items():
        tsv_path = metrics_tsv_dir / tsv_name
        if _tsv_row_count(tsv_path) > 0:
            continue
        log.debug(f"[fallback] Writing default rows for '{tsv_name}' "
                  f"({len(utr_ids)} transcripts).")
        metrics_tsv_dir.mkdir(parents=True, exist_ok=True)
        rows = [row_factory(sid) for sid in utr_ids]
        with open(tsv_path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=columns, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)


def evaluate_halflife(
    population: List[str], utr_ids: List[str],
    fixed_5utr_id: str, fixed_5utr_seq: str,
    fixed_cds_id: str, fixed_cds_seq: str,
    species: str, metrics_script: Path, predict_script: Path, work_dir: Path,
) -> Optional[Dict[str, float]]:
    """Run the metrics+LightGBM pipeline; return {sample_id: predicted_halflife}."""
    fasta_5, fasta_c, fasta_3 = work_dir / 'pop_5utr.fa', work_dir / 'pop_cds.fa', work_dir / 'pop_3utr.fa'
    write_fasta(fasta_5, [(sid, fixed_5utr_seq) for sid in utr_ids])
    write_fasta(fasta_c, [(fixed_cds_id, fixed_cds_seq)])
    write_fasta(fasta_3, list(zip(utr_ids, population)))

    metrics_out = work_dir / 'metrics_run'
    metrics_tsv_dir = metrics_out / 'metrics'
    if metrics_tsv_dir.exists():
        shutil.rmtree(metrics_tsv_dir)
    metrics_out.mkdir(parents=True, exist_ok=True)

    if not run_metrics(metrics_script, fasta_5, fasta_c, fasta_3, species, metrics_out, force=True):
        return None

    for unused in ('nmd_fragility_core.tsv', 'nmd_fragility_window.tsv'):
        p = metrics_tsv_dir / unused
        if p.exists():
            p.unlink()

    _synthesise_fallback_tsvs(metrics_tsv_dir, utr_ids, len(fixed_cds_seq))

    preds = run_prediction(predict_script, metrics_tsv_dir)
    if preds is None:
        return None
    pred_map = dict(zip(preds['transcript_id'], preds['predicted_halflife']))
    return {sid: pred_map.get(sid) for sid in utr_ids}


# ══════════════════════════════════════════════════════════════════════════
# Structural constraint: cmsearch against the original 3' UTR CM
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class CmsearchResult:
    query_name: str
    target_name: str
    score: str
    evalue: str


def _parse_tblout(tblout_file: str) -> List[dict]:
    hits: List[dict] = []
    if not os.path.exists(tblout_file):
        return hits
    with open(tblout_file) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split(maxsplit=17)
            if len(fields) < 17:
                continue
            hits.append({
                "target_name": fields[0], "query_name": fields[2],
                "score": fields[14], "evalue": fields[15],
            })
    return hits


def run_cmsearch(queries, cm_file, output_dir=".") -> List[CmsearchResult]:
    """Run cmsearch permissively (-E 1000); filtering to a hit/no-hit call happens later."""
    if shutil.which("cmsearch") is None:
        log.warning("cmsearch not found on PATH — skipping CM scoring for this generation.")
        return []
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, "cmsearch.out")
    tblout_file = os.path.join(output_dir, "cmsearch.tblout")
    sto_file = os.path.join(output_dir, "cmsearch.sto")

    results: List[CmsearchResult] = []
    with tempfile.TemporaryDirectory() as tmpdir:
        query_fa = os.path.join(tmpdir, "queries.fa")
        with open(query_fa, "w") as fh:
            for qname, qseq in queries:
                fh.write(f">{qname}\n{qseq}\n")
        cmd = ["cmsearch", "--notextw", "-A", sto_file, "-o", out_file,
               "--tblout", tblout_file, "-E", "1000", cm_file, query_fa]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            log.warning(f"cmsearch error:\n{r.stderr}")
        for h in _parse_tblout(tblout_file):
            results.append(CmsearchResult(query_name=h["query_name"],
                                           target_name=h["target_name"],
                                           score=h["score"], evalue=h["evalue"]))
    return results


def evaluate_cmscores(
    cm_model: Path, utr_ids: List[str], population: List[str],
    work_dir: Path, evalue_threshold: float, generation: int,
) -> Tuple[Dict[str, Optional[float]], Dict[str, bool]]:
    """Run cmsearch for the population; return ({sid: bitscore|None}, {sid: hit_flag})."""
    queries = list(zip(utr_ids, population))
    cmsearch_out_dir = work_dir / 'cmsearch_logs' / f"gen_{generation:04d}"
    results = run_cmsearch(queries, str(cm_model), str(cmsearch_out_dir))

    best_by_target: Dict[str, CmsearchResult] = {}
    for r in results:
        try:
            ev = float(r.evalue)
        except (TypeError, ValueError):
            continue
        current = best_by_target.get(r.target_name)
        if current is None or ev < float(current.evalue):
            best_by_target[r.target_name] = r

    cmscores: Dict[str, Optional[float]] = {}
    hit_flags: Dict[str, bool] = {}
    for sid in utr_ids:
        r = best_by_target.get(sid)
        if r is not None and float(r.evalue) <= evalue_threshold:
            cmscores[sid] = float(r.score)
            hit_flags[sid] = True
        else:
            cmscores[sid] = None
            hit_flags[sid] = False

    n_hits = sum(hit_flags.values())
    log.info(f"  cmsearch: {n_hits}/{len(utr_ids)} sequences hit the original "
             f"structure CM (E <= {evalue_threshold})")
    return cmscores, hit_flags


def evaluate_cmscore_single(seq_id: str, seq: str, cm_model: Path,
                             work_dir: Path, evalue_threshold: float) -> Tuple[Optional[float], bool]:
    """Convenience wrapper for scoring a single sequence (used for the seed baseline)."""
    cmscores, hit_flags = evaluate_cmscores(cm_model, [seq_id], [seq], work_dir, evalue_threshold, generation=0)
    return cmscores[seq_id], hit_flags[seq_id]


# ══════════════════════════════════════════════════════════════════════════
# Point scoring, softmax probabilities, weighted sampling with replacement
# ══════════════════════════════════════════════════════════════════════════

def compute_scores(
    utr_ids: List[str],
    rmsd_now: Dict[str, Optional[float]],
    rmsd_parent: Dict[str, Optional[float]],
    cmscore_now: Dict[str, Optional[float]],
    cmscore_parent: Dict[str, Optional[float]],
    hit_flags: Dict[str, bool],
    max_pids: Dict[str, Optional[float]],
    patent_pid_threshold: float,
    no_hit_penalty: float,
    patent_pid_penalty: float,
) -> Dict[str, int]:
    """
    Point score per individual:
      +1  RMSD (min over 3 refs) decreased vs. the parent that produced it
      +1  cmsearch bit score increased vs. the parent (only if there IS a hit)
      no_hit_penalty (default -10) if there is NO cmsearch hit at all,
          overriding the possible RMSD/cmscore +1's from cmscore's side
      patent_pid_penalty (default -10), additionally, if PID > threshold
          to any patent sequence
    Half-life is intentionally absent — reported elsewhere, never scored.
    """
    scores: Dict[str, int] = {}
    for sid in utr_ids:
        s = 0
        r_now, r_par = rmsd_now.get(sid), rmsd_parent.get(sid)
        if r_now is not None and r_par is not None and r_now < r_par:
            s += 1

        if not hit_flags.get(sid, False):
            s += no_hit_penalty
        else:
            c_now, c_par = cmscore_now.get(sid), cmscore_parent.get(sid)
            if c_now is not None and c_par is not None and c_now > c_par:
                s += 1

        pid = max_pids.get(sid)
        if pid is not None and pid > patent_pid_threshold:
            s += patent_pid_penalty

        scores[sid] = s
    return scores


def scores_to_probabilities(scores: List[int], temperature: float = 1.0) -> List[float]:
    """Softmax conversion of integer point scores → sampling probabilities."""
    if temperature <= 0:
        raise ValueError("--temperature must be > 0")
    m = max(scores)
    exps = [math.exp((s - m) / temperature) for s in scores]
    total = sum(exps)
    return [e / total for e in exps]


def sample_with_replacement(items: List[str], probabilities: List[float],
                             n: int, rng: random.Random) -> List[str]:
    """Weighted sampling WITH REPLACEMENT — duplicates are expected/allowed."""
    return rng.choices(items, weights=probabilities, k=n)


def breed_next_generation(
    selected_ids: List[str],
    seq_by_id: Dict[str, str],
    rmsd_by_id: Dict[str, Optional[float]],
    cmscore_by_id: Dict[str, Optional[float]],
    population_size: int,
    mutation_rate: float,
    crossover_rate: float,
    rng: random.Random,
) -> Tuple[List[str], List[str], List[Optional[float]], List[Optional[float]]]:
    """
    Fill the next generation (population_size individuals) by repeatedly
    drawing a parent (uniformly) from the n-select pool and mutating it —
    only a constant point-mutation rate is applied, no annealing/boosting.
    Returns (children_seqs, parent_of_child, parent_rmsd_of_child,
    parent_cmscore_of_child) aligned by index, so next generation's
    scoring has a baseline to compare against.
    """
    children, parent_of, parent_rmsd, parent_cm = [], [], [], []
    for _ in range(population_size):
        pa = rng.choice(selected_ids)
        seq = seq_by_id[pa]
        parent_id = pa
        if crossover_rate > 0 and rng.random() < crossover_rate and len(selected_ids) > 1:
            pb = rng.choice(selected_ids)
            seq = uniform_crossover(seq, seq_by_id[pb], rng)
        child = point_mutate(seq, mutation_rate, rng)
        children.append(child)
        parent_of.append(parent_id)
        parent_rmsd.append(rmsd_by_id.get(parent_id))
        parent_cm.append(cmscore_by_id.get(parent_id))
    return children, parent_of, parent_rmsd, parent_cm


# ══════════════════════════════════════════════════════════════════════════
# Plotting
# ══════════════════════════════════════════════════════════════════════════

def plot_scores_over_generations(best_per_gen: List[dict], out_path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        log.warning("matplotlib not installed — skipping plot.")
        return
    try:
        gens = [r['generation'] for r in best_per_gen]
        best_score = [r['best_score'] for r in best_per_gen]
        mean_score = [r['mean_score'] for r in best_per_gen]
        worst_score = [r['worst_score'] for r in best_per_gen]
        mean_rmsd = [r['mean_rmsd'] if r['mean_rmsd'] is not None else float('nan') for r in best_per_gen]
        mean_cm = [r['mean_cmscore'] if r['mean_cmscore'] is not None else float('nan') for r in best_per_gen]
        mean_halflife = [r['mean_halflife'] if r['mean_halflife'] is not None else float('nan') for r in best_per_gen]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8.5), dpi=150, sharex=True)

        ax1.plot(gens, best_score, label='Best score', color='#1b9e3e', marker='o', markersize=3)
        ax1.plot(gens, mean_score, label='Mean score', color='#2166ac', marker='o', markersize=3)
        ax1.plot(gens, worst_score, label='Worst score', color='#b2182b', linestyle='--', alpha=0.7)
        ax1.set_ylabel('Point score')
        ax1.set_title("3' UTR GA — Point Score (RMSD + cmscore + patent PID) Over Generations")
        ax1.legend(loc='best', framealpha=0.9)
        ax1.grid(True, alpha=0.3)

        ax2.plot(gens, mean_rmsd, label='Mean min-RMSD (Å)', color='#762a83', marker='o', markersize=3)
        ax2.set_ylabel('RMSD (Å)')
        ax2.tick_params(axis='y', labelcolor='#762a83')

        ax2b = ax2.twinx()
        ax2b.plot(gens, mean_cm, label='Mean cmscore (bits)', color='#e08214', marker='s', markersize=3)
        ax2b.set_ylabel('cmscore (bits)')

        ax2c = ax2.twinx()
        ax2c.spines['right'].set_position(('outward', 60))
        ax2c.plot(gens, mean_halflife, label='Mean predicted half-life (unweighted)',
                  color='#4d4d4d', linestyle=':', marker='^', markersize=3)
        ax2c.set_ylabel('Predicted half-life')

        l1, lb1 = ax2.get_legend_handles_labels()
        l2, lb2 = ax2b.get_legend_handles_labels()
        l3, lb3 = ax2c.get_legend_handles_labels()
        ax2.legend(l1 + l2 + l3, lb1 + lb2 + lb3, loc='best', framealpha=0.9, fontsize=8)
        ax2.set_xlabel('Generation')
        ax2.grid(True, alpha=0.3)
        if len(gens) <= 30:
            ax2.set_xticks(gens)

        fig.tight_layout()
        fig.savefig(out_path)
        plt.close(fig)
        log.info(f"Wrote {out_path}")
    except Exception as exc:
        log.warning(f"Failed to generate plot ({exc}); continuing without it.")


# ══════════════════════════════════════════════════════════════════════════
# Genetic Algorithm
# ══════════════════════════════════════════════════════════════════════════

def run_ga(
    seed_3utr_seq: str,
    fixed_5utr_id: str, fixed_5utr_seq: str,
    fixed_cds_id: str, fixed_cds_seq: str,
    species: str,
    metrics_script: Path, predict_script: Path,
    ref_cifs: List[Path],
    cm_model: Path,
    cm_evalue_threshold: float,
    no_hit_penalty: int,
    patent_seqs: List[dict],
    patent_pid_threshold: float,
    patent_pid_penalty: int,
    af3_work_dir: Path, af3_db: str, af3_models: str, af3_module: str,
    af3_partition: str, af3_time: str, af3_mem: str, af3_cpus: int, af3_gres: str,
    af3_exclude: str, af3_array_max_parallel: Optional[int], af3_model_seeds: List[int],
    af3_model_glob: str, af3_poll_interval: float, af3_timeout: Optional[float],
    atom_name: str, run_dssr: bool, dssr_path: str,
    output_dir: Path,
    population_size: int, n_select: int, generations: int,
    mutation_rate: float, crossover_rate: float, temperature: float,
    rng_seed: int, make_plot: bool,
) -> None:
    rng = random.Random(rng_seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    results_dir = output_dir / 'results'
    results_dir.mkdir(exist_ok=True)
    work_dir = output_dir / '_workspace'
    work_dir.mkdir(exist_ok=True)
    dssr_outdir = work_dir / 'dssr_logs'

    log.info("═" * 60)
    log.info("3' UTR Genetic Algorithm — point-scored, AF3-in-the-loop")
    log.info(f"  Population        : {population_size}")
    log.info(f"  n-select          : {n_select}  (weighted sampling WITH replacement)")
    log.info(f"  Generations       : {generations}")
    log.info(f"  Mutation rate     : {mutation_rate}  (flat — no annealing/boosting)")
    log.info(f"  Crossover rate    : {crossover_rate}")
    log.info(f"  Temperature       : {temperature}")
    log.info(f"  Reference CIFs    : {len(ref_cifs)}")
    log.info(f"  CM model          : {cm_model}")
    log.info(f"  cmsearch E-value  : <= {cm_evalue_threshold}")
    log.info(f"  No-hit penalty    : {no_hit_penalty}")
    log.info(f"  Patent sequences  : {len(patent_seqs)} motifs")
    log.info(f"  Patent PID thresh : > {patent_pid_threshold}% → penalty {patent_pid_penalty}")
    log.info(f"  Seed 3UTR         : {len(seed_3utr_seq)} nt")
    log.info(f"  Output            : {output_dir}")
    log.info("═" * 60)

    # ── Baseline (generation-0 "parent") metrics for the seed sequence ──
    log.info("Computing baseline (seed) metrics for generation-1 comparisons...")
    seed_id = "seed_3utr"
    jobs_tsv, out_root, job_id = submit_af3_population(
        [seed_3utr_seq], [seed_id], af3_work_dir, af3_db, af3_models, af3_module,
        af3_partition, af3_time, af3_mem, af3_cpus, af3_gres, af3_exclude,
        af3_array_max_parallel, af3_model_seeds, generation=0,
    )
    wait_for_slurm_job(job_id, af3_poll_interval, af3_timeout)
    seed_cif = collect_af3_cifs(out_root, [seed_id], af3_model_glob)
    seed_rmsd_map, _seed_dbn = evaluate_af3_rmsd(
        seed_cif, [seed_id], ref_cifs, atom_name, run_dssr, dssr_path, dssr_outdir)
    seed_cmscore_map, seed_hit_map = evaluate_cmscores(
        cm_model, [seed_id], [seed_3utr_seq], work_dir, cm_evalue_threshold, generation=0)

    seed_rmsd = seed_rmsd_map.get(seed_id)
    seed_cmscore = seed_cmscore_map.get(seed_id)
    log.info(f"  Seed baseline: RMSD={seed_rmsd}, cmscore={seed_cmscore}, "
             f"cm_hit={seed_hit_map.get(seed_id)}")

    # ── Generation-1 population: point mutations of the seed only ──────
    population = [point_mutate(seed_3utr_seq, mutation_rate, rng) for _ in range(population_size)]
    parent_of = [seed_id] * population_size
    parent_rmsd = [seed_rmsd] * population_size
    parent_cmscore = [seed_cmscore] * population_size

    all_rows: List[dict] = []
    selected_rows: List[dict] = []
    best_per_gen: List[dict] = []
    all_time_best_score = None
    all_time_best_seq = seed_3utr_seq
    all_time_best_id = seed_id
    all_time_best_gen = 0

    for gen in range(1, generations + 1):
        gen_start = time.time()
        log.info(f"\n── Generation {gen}/{generations} "
                 f"(mutation rate = {mutation_rate}) ──")

        utr_ids = [f"ind_{gen:04d}_{i:05d}" for i in range(population_size)]
        seq_by_id = dict(zip(utr_ids, population))
        parent_of_by_id = dict(zip(utr_ids, parent_of))
        parent_rmsd_by_id = dict(zip(utr_ids, parent_rmsd))
        parent_cm_by_id = dict(zip(utr_ids, parent_cmscore))

        # (1) AlphaFold3 structures, parallelized on SLURM as one array job
        jobs_tsv, out_root, job_id = submit_af3_population(
            population, utr_ids, af3_work_dir, af3_db, af3_models, af3_module,
            af3_partition, af3_time, af3_mem, af3_cpus, af3_gres, af3_exclude,
            af3_array_max_parallel, af3_model_seeds, generation=gen,
        )
        wait_for_slurm_job(job_id, af3_poll_interval, af3_timeout)
        cif_paths = collect_af3_cifs(out_root, utr_ids, af3_model_glob)

        # (2) RMSD vs. the 3 reference structures (min), + DSSR dot-bracket
        rmsd_now, dbn_now = evaluate_af3_rmsd(
            cif_paths, utr_ids, ref_cifs, atom_name, run_dssr, dssr_path, dssr_outdir)

        # (3) cmsearch bit score
        cmscore_now, hit_flags = evaluate_cmscores(
            cm_model, utr_ids, population, work_dir, cm_evalue_threshold, generation=gen)

        # (4) Patent PID (Smith-Waterman)
        max_pids: Dict[str, Optional[float]] = {}
        all_pids: Dict[str, List[float]] = {}
        for sid, seq in zip(utr_ids, population):
            pids = check_patent_pid(seq, patent_seqs) if patent_seqs else []
            all_pids[sid] = pids
            max_pids[sid] = max(pids) if pids else None

        # (5) Predicted half-life — reported only, never weighted into score
        halflife_map = evaluate_halflife(
            population, utr_ids, fixed_5utr_id, fixed_5utr_seq,
            fixed_cds_id, fixed_cds_seq, species, metrics_script, predict_script, work_dir,
        ) or {sid: None for sid in utr_ids}

        # ── Score, relative to each individual's own parent ────────────
        scores = compute_scores(
            utr_ids, rmsd_now, parent_rmsd_by_id, cmscore_now, parent_cm_by_id,
            hit_flags, max_pids, patent_pid_threshold, no_hit_penalty, patent_pid_penalty,
        )
        score_list = [scores[sid] for sid in utr_ids]
        probabilities = scores_to_probabilities(score_list, temperature)
        prob_by_id = dict(zip(utr_ids, probabilities))

        # ── Record every individual this generation ────────────────────
        for sid in utr_ids:
            all_rows.append({
                'generation': gen,
                'sample_id': sid,
                'parent_id': parent_of_by_id[sid],
                'score': scores[sid],
                'probability': prob_by_id[sid],
                'rmsd_min': rmsd_now.get(sid),
                'parent_rmsd': parent_rmsd_by_id[sid],
                'cmscore': cmscore_now.get(sid),
                'parent_cmscore': parent_cm_by_id[sid],
                'cm_hit': hit_flags.get(sid, False),
                'max_patent_pid': max_pids.get(sid),
                'predicted_halflife': halflife_map.get(sid),
                'dssr_dbn': dbn_now.get(sid),
                'sequence': seq_by_id[sid],
                'length': len(seq_by_id[sid]),
                'mutations_from_seed': sum(
                    a != b for a, b in zip(seq_by_id[sid], seed_3utr_seq)
                    if len(seq_by_id[sid]) == len(seed_3utr_seq)),
            })

        gen_elapsed = time.time() - gen_start
        valid_rmsd = [v for v in rmsd_now.values() if v is not None]
        valid_cm = [v for v in cmscore_now.values() if v is not None]
        valid_hl = [v for v in halflife_map.values() if v is not None]
        best_sid = max(utr_ids, key=lambda s: scores[s])

        log.info(f"  Best score this gen   : {scores[best_sid]}  ({best_sid})")
        log.info(f"  Score  best/mean/worst: {max(score_list)}/"
                 f"{sum(score_list)/len(score_list):.2f}/{min(score_list)}")
        log.info(f"  CM hits               : {sum(hit_flags.values())}/{population_size}")
        if patent_seqs:
            n_patent_pass = sum(1 for p in max_pids.values() if p is None or p <= patent_pid_threshold)
            log.info(f"  Patent PID pass       : {n_patent_pass}/{population_size}")
        log.info(f"  Elapsed               : {gen_elapsed:.1f}s")

        best_per_gen.append({
            'generation': gen,
            'best_score': max(score_list),
            'mean_score': sum(score_list) / len(score_list),
            'worst_score': min(score_list),
            'mean_rmsd': (sum(valid_rmsd) / len(valid_rmsd)) if valid_rmsd else None,
            'mean_cmscore': (sum(valid_cm) / len(valid_cm)) if valid_cm else None,
            'mean_halflife': (sum(valid_hl) / len(valid_hl)) if valid_hl else None,
            'best_id': best_sid,
            'best_sequence': seq_by_id[best_sid],
        })

        if all_time_best_score is None or scores[best_sid] > all_time_best_score:
            all_time_best_score = scores[best_sid]
            all_time_best_seq = seq_by_id[best_sid]
            all_time_best_id = best_sid
            all_time_best_gen = gen
            log.info(f"  ★ New all-time best score: {all_time_best_score}")

        # ── Select n_select individuals WITH REPLACEMENT, by probability ──
        selected_ids = sample_with_replacement(utr_ids, probabilities, n_select, rng)
        for rank, sid in enumerate(selected_ids):
            selected_rows.append({
                'generation': gen, 'rank': rank + 1, 'sample_id': sid,
                'score': scores[sid], 'probability': prob_by_id[sid],
                'sequence': seq_by_id[sid],
            })
        log.info(f"  Selected (with replacement): "
                 f"{', '.join(f'{s}(score={scores[s]})' for s in selected_ids)}")

        if gen == generations:
            break

        # ── Breed the next population purely by mutation (+ optional
        #    crossover) from the selected pool ─────────────────────────
        population, parent_of, parent_rmsd, parent_cmscore = breed_next_generation(
            selected_ids, seq_by_id, rmsd_now, cmscore_now,
            population_size, mutation_rate, crossover_rate, rng,
        )

    # ══ Write outputs ═══════════════════════════════════════════════════
    bpg_path = results_dir / 'best_per_generation.tsv'
    with open(bpg_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(best_per_gen[0].keys()),
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(best_per_gen)
    log.info(f"Wrote {bpg_path}")

    all_path = results_dir / 'all_generations.tsv'
    with open(all_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()),
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(all_rows)
    log.info(f"Wrote {all_path}")

    sel_path = results_dir / 'selected_per_generation.tsv'
    with open(sel_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(selected_rows[0].keys()),
                                delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(selected_rows)
    log.info(f"Wrote {sel_path}")

    best_fa = results_dir / 'best_3utr.fa'
    write_fasta(best_fa, [('best_3utr_evolved', all_time_best_seq)])
    log.info(f"Wrote {best_fa}")

    seed_fa = results_dir / 'seed_3utr.fa'
    write_fasta(seed_fa, [('seed_3utr_original', seed_3utr_seq)])

    if make_plot:
        plot_scores_over_generations(best_per_gen, results_dir / 'score_over_generations.png')

    summary_path = results_dir / 'evolution_summary.txt'
    with open(summary_path, 'w') as fh:
        fh.write("3' UTR Genetic Algorithm — Evolution Summary\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(f"Generations          : {generations}\n")
        fh.write(f"Population size      : {population_size}\n")
        fh.write(f"n-select             : {n_select}\n")
        fh.write(f"Mutation rate (flat) : {mutation_rate}\n")
        fh.write(f"Crossover rate       : {crossover_rate}\n")
        fh.write(f"Softmax temperature  : {temperature}\n")
        fh.write(f"RNG seed             : {rng_seed}\n")
        fh.write(f"Reference CIFs       : {[str(c) for c in ref_cifs]}\n")
        fh.write(f"CM model             : {cm_model}\n")
        fh.write(f"cmsearch E-value     : <= {cm_evalue_threshold}\n")
        fh.write(f"No-hit penalty       : {no_hit_penalty}\n")
        fh.write(f"Patent sequences     : {len(patent_seqs)} motifs\n")
        fh.write(f"Patent PID threshold : > {patent_pid_threshold}%\n")
        fh.write(f"Patent penalty       : {patent_pid_penalty}\n\n")
        fh.write(f"Seed 3UTR length : {len(seed_3utr_seq)} nt\n")
        fh.write(f"Seed RMSD/cmscore: {seed_rmsd} / {seed_cmscore}\n")
        fh.write(f"Best 3UTR length : {len(all_time_best_seq)} nt\n\n")
        fh.write(f"All-time best score : {all_time_best_score} "
                 f"(id={all_time_best_id}, generation={all_time_best_gen})\n\n")
        fh.write("Generation-by-generation:\n")
        for row in best_per_gen:
            fh.write(
                f"  Gen {row['generation']:>3d} : score best {row['best_score']} "
                f"(mean {row['mean_score']:.2f}, worst {row['worst_score']}) | "
                f"mean_rmsd {row['mean_rmsd']} | mean_cmscore {row['mean_cmscore']} | "
                f"mean_halflife(unweighted) {row['mean_halflife']}\n"
            )
        fh.write("\nBest evolved 3UTR sequence:\n")
        for i in range(0, len(all_time_best_seq), 60):
            fh.write(all_time_best_seq[i:i + 60] + '\n')
    log.info(f"Wrote {summary_path}")

    log.info("\n" + "═" * 60)
    log.info(f"GA complete. All-time best score: {all_time_best_score} "
             f"(id={all_time_best_id}, generation={all_time_best_gen})")
    log.info(f"Results in: {results_dir}")
    log.info("═" * 60)


# ══════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════

def main():
    global af3prep, af3cmp
    af3prep = _import_sibling('af3prep', 'prepare_af3_jobs.py')
    af3cmp = _import_sibling('af3cmp', 'compare_af3_candidates.py')

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    # ── Required biological inputs ──────────────────────────────────────
    parser.add_argument('--fasta-5utr', required=True, dest='fasta_5utr', metavar='FILE',
                        help="FASTA with a single 5' UTR sequence (fixed, not mutated).")
    parser.add_argument('--fasta-cds', required=True, dest='fasta_cds', metavar='FILE',
                        help="FASTA with a single CDS sequence (fixed, not mutated).")
    parser.add_argument('--fasta-3utr', required=True, dest='fasta_3utr', metavar='FILE',
                        help="FASTA with the seed 3' UTR sequence to evolve.")
    parser.add_argument('--species', required=True, metavar='BINOMIAL',
                        help="Species name passed to the metrics script.")
    parser.add_argument('--metrics-script', required=True, dest='metrics_script', metavar='FILE',
                        help="Path to 01b_metrics.py.")
    parser.add_argument('--predict-script', required=True, dest='predict_script', metavar='FILE',
                        help="Path to the LightGBM half-life prediction script.")

    # ── Structural scoring inputs ───────────────────────────────────────
    parser.add_argument('--ref-cif', required=True, action='append', dest='ref_cifs', metavar='FILE',
                        help="A representative reference .cif. Repeat exactly 3 times "
                             "(one per --ref-cif); RMSD for each individual is the "
                             "minimum Needleman-Wunsch Biopython RMSD across all of them.")
    parser.add_argument('--atom', default="C1'", metavar='NAME',
                        help="Representative per-residue atom for RMSD (default: %(default)s).")
    parser.add_argument('--skip-dssr', action='store_true',
                        help="Don't run x3dna-dssr on candidate structures (dot-bracket "
                             "reporting only — scoring is unaffected either way).")
    parser.add_argument('--dssr-path', default='x3dna-dssr', metavar='PATH',
                        help="Path to the x3dna-dssr executable (default: %(default)s).")
    parser.add_argument('--cm-model', required=True, dest='cm_model', metavar='FILE',
                        help="Calibrated Infernal covariance model (.cm) of the original "
                             "3' UTR structure.")
    parser.add_argument('--cm-evalue-threshold', type=float, default=0.01,
                        dest='cm_evalue_threshold', metavar='EVALUE',
                        help="E-value cutoff for a cmsearch hit (default: 0.01).")
    parser.add_argument('--no-hit-penalty', type=int, default=-10, dest='no_hit_penalty', metavar='N',
                        help="Score penalty when there is no cmsearch hit (default: -10).")

    # ── Patent PID inputs ────────────────────────────────────────────────
    parser.add_argument('--patent-pid-threshold', type=float, default=80.0,
                        dest='patent_pid_threshold', metavar='PCT',
                        help="PID (%%) above which the patent penalty applies (default: 80.0).")
    parser.add_argument('--patent-pid-penalty', type=int, default=-10,
                        dest='patent_pid_penalty', metavar='N',
                        help="Score penalty for PID > threshold to any patent sequence (default: -10).")
    parser.add_argument('--no-default-patents', action='store_true', dest='no_default_patents',
                        help="Disable the built-in A7 patent sequence set.")

    # ── AlphaFold3 / SLURM ───────────────────────────────────────────────
    parser.add_argument('--af3-work-dir', required=True, dest='af3_work_dir', metavar='DIR',
                        help="AF3 working directory (AF3_WD); one gen_NNNN/ subdir per generation.")
    parser.add_argument('--af3-db', default='/opt/alphafold3/databases', dest='af3_db',
                        help="Path for AF3_DB (default: %(default)s).")
    parser.add_argument('--af3-models', required=True, dest='af3_models',
                        help="Path for AF3_MODELS (model weights).")
    parser.add_argument('--af3-module', default='alphafold3/3.0.3', dest='af3_module',
                        help="Environment module to load (default: %(default)s).")
    parser.add_argument('--af3-model-seed', type=int, action='append', default=None,
                        dest='af3_model_seeds', help="AF3 model seed(s); repeatable (default: [1]).")
    parser.add_argument('--af3-model-glob', default='**/*_model.cif', dest='af3_model_glob',
                        metavar='GLOB',
                        help="Glob (relative to each job's output dir) used to find the top "
                             "model .cif (default: %(default)s). Adjust to match your AF3 "
                             "version's output layout if nothing is found.")
    parser.add_argument('--af3-partition', default='aoraki_gpu_L40', dest='af3_partition',
                        help="SLURM partition (default: %(default)s).")
    parser.add_argument('--af3-time', default='02:00:00', dest='af3_time',
                        help="SLURM time (default: %(default)s).")
    parser.add_argument('--af3-mem', default='32G', dest='af3_mem',
                        help="SLURM memory (default: %(default)s).")
    parser.add_argument('--af3-cpus-per-task', type=int, default=8, dest='af3_cpus',
                        help="SLURM cpus-per-task (default: %(default)s).")
    parser.add_argument('--af3-gres', default='gpu:1', dest='af3_gres',
                        help="SLURM gres (default: %(default)s).")
    parser.add_argument('--af3-exclude', default='', dest='af3_exclude',
                        help="Nodes to exclude, comma-separated (default: none).")
    parser.add_argument('--af3-array-max-parallel', type=int, default=None,
                        dest='af3_array_max_parallel', metavar='N',
                        help="Cap on simultaneously running array tasks (default: no cap).")
    parser.add_argument('--af3-poll-interval', type=float, default=30.0, dest='af3_poll_interval',
                        help="Seconds between squeue polls while waiting for AF3 (default: 30).")
    parser.add_argument('--af3-timeout', type=float, default=None, dest='af3_timeout',
                        help="Max seconds to wait for a generation's AF3 array job "
                             "(default: no timeout).")

    # ── GA parameters ────────────────────────────────────────────────────
    parser.add_argument('--population', type=int, default=1000, metavar='N',
                        help="Population size per generation (default: 1000).")
    parser.add_argument('--n-select', type=int, default=10, dest='n_select', metavar='N',
                        help="Individuals drawn WITH REPLACEMENT each generation, "
                             "weighted by softmax(score), to seed the next generation "
                             "(default: 10).")
    parser.add_argument('--generations', type=int, default=30, metavar='N',
                        help="Number of GA generations (default: 30).")
    parser.add_argument('--mutation-rate', type=float, default=0.02, dest='mutation_rate', metavar='RATE',
                        help="Constant per-nucleotide mutation probability, applied "
                             "identically every generation — no annealing or boosting "
                             "of any kind (default: 0.02).")
    parser.add_argument('--crossover-rate', type=float, default=0.0, dest='crossover_rate', metavar='RATE',
                        help="Probability of uniform crossover between two selected "
                             "parents before mutation (default: 0.0, i.e. mutation only).")
    parser.add_argument('--temperature', type=float, default=1.0, metavar='T',
                        help="Softmax temperature for score→probability conversion; "
                             "lower = more concentrated on high scorers (default: 1.0).")
    parser.add_argument('--seed', type=int, default=42, metavar='INT',
                        help="Random seed for reproducibility (default: 42).")

    # ── Output ────────────────────────────────────────────────────────────
    parser.add_argument('--output-dir', '-o', default='ga_output', dest='output_dir', metavar='DIR')
    parser.add_argument('--no-plot', action='store_true', dest='no_plot')
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if len(args.ref_cifs) != 3:
        log.error(f"--ref-cif must be given exactly 3 times (got {len(args.ref_cifs)}).")
        sys.exit(1)

    required_paths = [
        ('--fasta-5utr', args.fasta_5utr), ('--fasta-cds', args.fasta_cds),
        ('--fasta-3utr', args.fasta_3utr), ('--metrics-script', args.metrics_script),
        ('--predict-script', args.predict_script), ('--cm-model', args.cm_model),
    ] + [('--ref-cif', c) for c in args.ref_cifs]
    for label, path_str in required_paths:
        if not Path(path_str).exists():
            log.error(f"{label} not found: {path_str}")
            sys.exit(1)

    if shutil.which("cmsearch") is None:
        log.error("'cmsearch' (Infernal) not found on PATH — required for CM scoring.")
        sys.exit(1)
    if shutil.which("sbatch") is None:
        log.error("'sbatch' not found on PATH — required to submit AF3 SLURM array jobs.")
        sys.exit(1)
    if not args.skip_dssr and shutil.which(args.dssr_path) is None:
        log.warning(f"'{args.dssr_path}' not found on PATH — DSSR dot-bracket reporting "
                    f"will be skipped (scoring is unaffected).")
        args.skip_dssr = True

    try:
        u5_id, u5_seq = read_single_fasta(Path(args.fasta_5utr))
        cds_id, cds_seq = read_single_fasta(Path(args.fasta_cds))
        u3_id, u3_seq = read_single_fasta(Path(args.fasta_3utr))
    except ValueError as exc:
        log.error(str(exc))
        sys.exit(1)

    log.info(f"5' UTR : '{u5_id}'  {len(u5_seq)} nt (fixed)")
    log.info(f"CDS    : '{cds_id}' {len(cds_seq)} nt (fixed)")
    log.info(f"3' UTR : '{u3_id}'  {len(u3_seq)} nt (seed — will be evolved)")

    run_ga(
        seed_3utr_seq=u3_seq,
        fixed_5utr_id=u5_id, fixed_5utr_seq=u5_seq,
        fixed_cds_id=cds_id, fixed_cds_seq=cds_seq,
        species=args.species,
        metrics_script=Path(args.metrics_script), predict_script=Path(args.predict_script),
        ref_cifs=[Path(c) for c in args.ref_cifs],
        cm_model=Path(args.cm_model),
        cm_evalue_threshold=args.cm_evalue_threshold,
        no_hit_penalty=args.no_hit_penalty,
        patent_seqs=[] if args.no_default_patents else PATENT_SEQUENCES,
        patent_pid_threshold=args.patent_pid_threshold,
        patent_pid_penalty=args.patent_pid_penalty,
        af3_work_dir=Path(args.af3_work_dir), af3_db=args.af3_db, af3_models=args.af3_models,
        af3_module=args.af3_module, af3_partition=args.af3_partition, af3_time=args.af3_time,
        af3_mem=args.af3_mem, af3_cpus=args.af3_cpus, af3_gres=args.af3_gres,
        af3_exclude=args.af3_exclude, af3_array_max_parallel=args.af3_array_max_parallel,
        af3_model_seeds=args.af3_model_seeds or [1], af3_model_glob=args.af3_model_glob,
        af3_poll_interval=args.af3_poll_interval, af3_timeout=args.af3_timeout,
        atom_name=args.atom, run_dssr=not args.skip_dssr, dssr_path=args.dssr_path,
        output_dir=Path(args.output_dir),
        population_size=args.population, n_select=args.n_select, generations=args.generations,
        mutation_rate=args.mutation_rate, crossover_rate=args.crossover_rate,
        temperature=args.temperature, rng_seed=args.seed, make_plot=not args.no_plot,
    )


if __name__ == '__main__':
    main()