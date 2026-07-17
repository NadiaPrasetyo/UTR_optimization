#!/usr/bin/env python3
"""
mutate_3utr_ga.py — minimal edition
─────────────────────────────────────────────────────────────────────────────
Evolves a 3' UTR sequence over generations by point mutation (+ optional
crossover), scoring each individual on:
  +1  AlphaFold3 structure RMSD improved vs. its parent
  +1  cmsearch bit score improved vs. its parent (only if it hits the CM)
  no_hit_penalty  if it does NOT hit the structural CM at all
  patent_pid_penalty  if it's too similar (Smith-Waterman PID) to a patent motif
Predicted half-life is computed and reported every generation but never
included in the score.

RESUMABILITY: a checkpoint is written after every generation to
<output-dir>/_workspace/checkpoints/. Re-running the same command resumes
from the last checkpoint automatically. `rm -rf <output-dir>/_workspace/checkpoints`
to start over.

AF3 SCHEDULING: each generation's population is split into SLURM array
batches (`--af3-batch-size`, default 100) and up to `--af3-max-concurrent-
batches` (default 10) batches are submitted at once, each via a blocking
`sbatch --wait` call run in its own thread. There is no manual squeue/sacct
polling loop — a batch's thread simply blocks until SLURM says it's done,
then every individual in that batch is checked (exit code + presence of the
expected `.cif`) and logged immediately as OK or FAILED. Structures already
present on disk from a previous run are skipped and reused.
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
import pickle
import random
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

logging.basicConfig(level=logging.INFO,
                     format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
                     datefmt='%H:%M:%S')
log = logging.getLogger('ga_3utr')

NUCLEOTIDES = ['A', 'T', 'G', 'C']
RMSD_ATOM = "C1'"

PATENT_SEQUENCES = [
    {"id": "A7_30nt", "seq": "CAGACCCTGGTCCGGGGCAATGGGACCACT"},
    {"id": "A7_32nt", "seq": "TCAGACCCTGGTCCGGGGCAAATGGGACCACTG"},
    {"id": "A7_34nt", "seq": "GTCAGACCCTGGTCCGGGGCAATGGGACCACTGT"},
    {"id": "A7_36nt", "seq": "GGTCAGACCCTGGTCCGGGGCAAATGGGACCACTGTT"},
    {"id": "A7_40nt", "seq": "TGGGTCAGACCCTGGTCCGGGGCAAATGGGACCACTGTTTC"},
    {"id": "A7_43nt", "seq": "TGGGTCAGACCCTGGTCCGGGGCAATGGGACCACTGTTTCGCG"},
    {"id": "A7_47nt", "seq": "TGGGTCAGACCCTGGTCCGGGGCAATGGGACCACTGTTTCGCGTTTA"},
]


# ══════════════════════════════════════════════════════════════════════════
# Checkpoint/Resume
# ══════════════════════════════════════════════════════════════════════════

@dataclass
class GACheckpoint:
    generation: int
    population: List[str]
    parent_of: List[str]
    parent_rmsd: List[Optional[float]]
    parent_cmscore: List[Optional[float]]
    seed_rmsd: Optional[float]
    seed_cmscore: Optional[float]
    all_rows: List[dict]
    selected_rows: List[dict]
    best_per_gen: List[dict]
    all_time_best_score: Optional[int]
    all_time_best_seq: str
    all_time_best_id: str
    all_time_best_gen: int
    rng_state: tuple


def checkpoint_dir(output_dir: Path) -> Path:
    return output_dir / '_workspace' / 'checkpoints'


def find_latest_checkpoint(output_dir: Path) -> Optional[Path]:
    d = checkpoint_dir(output_dir)
    if not d.exists():
        return None
    ckpts = sorted(d.glob('checkpoint_gen_*.pkl'))
    return ckpts[-1] if ckpts else None


def save_checkpoint(ckpt: GACheckpoint, output_dir: Path) -> None:
    d = checkpoint_dir(output_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"checkpoint_gen_{ckpt.generation:04d}.pkl"
    with open(path, 'wb') as fh:
        pickle.dump(ckpt, fh)
    log.info(f"  Saved checkpoint: {path}")


def load_checkpoint(path: Path) -> GACheckpoint:
    with open(path, 'rb') as fh:
        ckpt = pickle.load(fh)
    log.info(f"  Loaded checkpoint from generation {ckpt.generation}: {path}")
    return ckpt


# ══════════════════════════════════════════════════════════════════════════
# Sibling-script imports
# ══════════════════════════════════════════════════════════════════════════

def _import_sibling(module_name: str, filename: str):
    here = Path(__file__).resolve().parent
    path = here / filename
    if not path.exists():
        raise FileNotFoundError(f"Required helper '{filename}' not found next to this script.")
    spec = importlib.util.spec_from_file_location(module_name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # type: ignore[union-attr]
    return mod


af3prep = None
af3cmp = None


# ══════════════════════════════════════════════════════════════════════════
# FASTA I/O
# ══════════════════════════════════════════════════════════════════════════

def read_single_fasta(path: Path) -> Tuple[str, str]:
    seq_id, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if seq_id is not None:
                    raise ValueError(f"Multiple sequences in {path}; expected exactly one.")
                seq_id = line[1:].split()[0]
            elif line:
                parts.append(line.upper().replace('U', 'T'))
    if seq_id is None:
        raise ValueError(f"No sequences found in {path}")
    return seq_id, ''.join(parts)


def write_fasta(path: Path, records: List[Tuple[str, str]], width: int = 60) -> None:
    with open(path, 'w') as fh:
        for seq_id, seq in records:
            fh.write(f'>{seq_id}\n')
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + '\n')


# ══════════════════════════════════════════════════════════════════════════
# Genetic operators (flat mutation rate — no annealing)
# ══════════════════════════════════════════════════════════════════════════

def point_mutate(seq: str, rate: float, rng: random.Random) -> str:
    bases = list(seq)
    for i, base in enumerate(bases):
        if rng.random() < rate:
            bases[i] = rng.choice([n for n in NUCLEOTIDES if n != base])
    return ''.join(bases)


def uniform_crossover(a: str, b: str, rng: random.Random) -> str:
    if len(a) != len(b):
        return a
    return ''.join(x if rng.random() < 0.5 else y for x, y in zip(a, b))


# ══════════════════════════════════════════════════════════════════════════
# Patent sequence identity via Smith-Waterman
# ══════════════════════════════════════════════════════════════════════════

def _sw_pure(a: str, b: str) -> dict:
    match, mismatch, gap = 2, -1, -1
    m, n = len(a), len(b)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    best, bi, bj = 0, 0, 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = H[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            s = max(s, H[i - 1][j] + gap, H[i][j - 1] + gap, 0)
            H[i][j] = s
            if s > best:
                best, bi, bj = s, i, j
    ai, aj = [], []
    i, j = bi, bj
    while i > 0 and j > 0 and H[i][j] > 0:
        if a[i - 1] == b[j - 1] or (H[i - 1][j - 1] >= H[i - 1][j] and H[i - 1][j - 1] >= H[i][j - 1]):
            ai.append(a[i - 1]); aj.append(b[j - 1]); i -= 1; j -= 1
        elif H[i - 1][j] > H[i][j - 1]:
            ai.append(a[i - 1]); aj.append('-'); i -= 1
        else:
            ai.append('-'); aj.append(b[j - 1]); j -= 1
    ai, aj = ''.join(reversed(ai)), ''.join(reversed(aj))
    identities = sum(x == y and x != '-' for x, y in zip(ai, aj))
    return {'identities': identities}


def calculate_pid(query: str, subject: str) -> float:
    shorter = min(len(query), len(subject))
    if shorter == 0:
        return 0.0
    try:
        from Bio import pairwise2
        from Bio.Seq import Seq
        aln = pairwise2.align.localms(Seq(query), Seq(subject), 2, -1, -1, -1)
        if aln:
            a, b = str(aln[0][0]), str(aln[0][1])
            ident = sum(x == y and x != '-' for x, y in zip(a, b))
            return (ident / shorter) * 100.0
    except Exception:
        pass
    return (_sw_pure(query, subject)['identities'] / shorter) * 100.0


def check_patent_pid(seq: str, patents: List[dict]) -> List[float]:
    return [calculate_pid(seq, p['seq']) for p in patents]


# ══════════════════════════════════════════════════════════════════════════
# AlphaFold3 on SLURM — batched, blocking (sbatch --wait), thread-pooled
# ══════════════════════════════════════════════════════════════════════════

def _build_af3_batch_script(job_name_paths: List[Tuple[str, Path]], jobs_dir: Path,
                             logs_dir: Path, output_root: Path, af3_db: str, af3_models: str,
                             af3_module: str, partition: str, time_limit: str, mem: str,
                             cpus: int, gres: str, exclude: str, generation: int,
                             batch_idx: int, gen_dir: Path, model_glob: str,
                             af3_per_gpu: int) -> Tuple[Path, Path]:
    """
    Write jobs.tsv + an sbatch array script for a batch. Every `af3_per_gpu`
    individuals share one array task (one GPU allocation) and run
    concurrently in the background on it — set this based on how many
    predictions your input size actually fits in one GPU's VRAM (check
    `nvidia-smi` while a single af3 run is in flight). Default 1 reproduces
    one-individual-per-GPU behavior exactly.
    """
    jobs_tsv = jobs_dir / "jobs.tsv"
    with open(jobs_tsv, "w") as fh:
        for idx, (job_name, json_path) in enumerate(job_name_paths):
            fh.write(f"{idx}\t{job_name}\t{json_path}\t{output_root / job_name}\n")

    n_slots = math.ceil(len(job_name_paths) / af3_per_gpu)
    exclude_line = f"#SBATCH --exclude={exclude}\n" if exclude else ""
    script = f"""#!/bin/bash
#SBATCH --job-name=af3_gen{generation:04d}_b{batch_idx:03d}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time_limit}
#SBATCH --partition={partition}
#SBATCH --gres={gres}
#SBATCH --array=0-{n_slots - 1}
#SBATCH --output={logs_dir}/af3_%A_%a.log
{exclude_line}
module purge
module load {af3_module}

export AF3_DB={af3_db}
export AF3_WD={gen_dir}
export AF3_MODELS={af3_models}

PER_GPU={af3_per_gpu}
SLOT_START=$((SLURM_ARRAY_TASK_ID * PER_GPU))
SLOT_END=$((SLOT_START + PER_GPU - 1))
LINES=$(awk -F'\\t' -v s="$SLOT_START" -v e="$SLOT_END" '$1 >= s && $1 <= e' "{jobs_tsv}")

# JAX/TF each preallocate a large default share of the GPU per process, so
# running PER_GPU>1 processes on one GPU needs each capped to ~1/PER_GPU of
# it (with a small safety margin) or they OOM before ever touching the
# actual input — independent of how small the sequence is.
MEM_FRACTION=$(awk -v n="$PER_GPU" 'BEGIN{{printf "%.4f", 0.9 / n}}')
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION="$MEM_FRACTION"
export TF_FORCE_GPU_ALLOW_GROWTH=true

echo "Slot $SLURM_ARRAY_TASK_ID: $(echo "$LINES" | wc -l) job(s) running concurrently on this GPU (XLA_PYTHON_CLIENT_MEM_FRACTION=$MEM_FRACTION)."

run_one() {{
    local job_name="$1" json_path="$2" out_dir="$3"
    mkdir -p "$out_dir"
    rm -f "$out_dir/.failed"
    echo "[$job_name] starting at $(date)"
    af3 "$json_path" "$out_dir" --run_data_pipeline=true --run_inference=true \\
        > "{logs_dir}/af3_task_${{job_name}}.log" 2>&1
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "[$job_name] FAILED at $(date) (exit $rc) — see {logs_dir}/af3_task_${{job_name}}.log"
        echo "$rc" > "$out_dir/.failed"
        return $rc
    fi
    local n_cifs
    n_cifs=$(find "$out_dir" -path "$out_dir/{model_glob}" 2>/dev/null | wc -l)
    if [ "$n_cifs" -eq 0 ]; then
        echo "[$job_name] exited 0 but no output matching '{model_glob}' found — treating as FAILED (likely OOM)."
        echo "no_output_despite_exit_0" > "$out_dir/.failed"
        return 1
    fi
    echo "[$job_name] finished at $(date), $n_cifs structure file(s) confirmed."
}}

slot_failed=0
while IFS=$'\\t' read -r idx job_name json_path out_dir; do
    [ -z "$job_name" ] && continue
    run_one "$job_name" "$json_path" "$out_dir" &
done <<< "$LINES"

for pid in $(jobs -p); do
    wait "$pid" || slot_failed=1
done

echo "Slot $SLURM_ARRAY_TASK_ID: all assigned jobs finished at $(date)."
[ "$slot_failed" -ne 0 ] && exit 1
exit 0
"""
    script_path = jobs_dir / "run_af3_array.sh"
    script_path.write_text(script)
    script_path.chmod(0o755)
    return jobs_tsv, script_path


def _find_cached_cif(job_out_dir: Path, model_glob: str) -> Optional[Path]:
    if not job_out_dir.exists():
        return None
    matches = sorted(job_out_dir.glob(model_glob))
    return matches[0] if matches else None


def _run_batch_blocking(batch_idx: int, batch: List[Tuple[str, Path]], jobs_root: Path,
                         output_root: Path, gen_dir: Path, af3_db: str, af3_models: str,
                         af3_module: str, partition: str, time_limit: str, mem: str,
                         cpus: int, gres: str, exclude: str, generation: int,
                         model_glob: str, af3_per_gpu: int) -> List[Tuple[str, bool, str]]:
    """Submit one batch with `sbatch --wait` (blocks this thread until it's done),
    then check + log every individual's outcome. Returns [(job_name, ok, reason)]."""
    jobs_dir = jobs_root / f"batch_{batch_idx:04d}"
    logs_dir = jobs_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    _, script_path = _build_af3_batch_script(
        batch, jobs_dir, logs_dir, output_root, af3_db, af3_models, af3_module,
        partition, time_limit, mem, cpus, gres, exclude, generation, batch_idx,
        gen_dir, model_glob, af3_per_gpu)

    log.info(f"  [gen {generation}] submitting batch {batch_idx} ({len(batch)} jobs)...")
    result = subprocess.run(["sbatch", "--wait", "--parsable", str(script_path)],
                             capture_output=True, text=True)
    job_id = result.stdout.strip().split(';')[0] or "?"
    if result.returncode != 0:
        log.error(f"  [gen {generation} batch {batch_idx}] sbatch failed: {result.stderr.strip()}")

    outcomes = []
    for job_name, _ in batch:
        out_dir = output_root / job_name
        failed_marker = out_dir / ".failed"
        cif = _find_cached_cif(out_dir, model_glob)
        ok = cif is not None and not failed_marker.exists()
        reason = "" if ok else (failed_marker.read_text().strip() if failed_marker.exists() else "no output found")
        outcomes.append((job_name, ok, reason))
        log.info(f"  [gen {generation} batch {batch_idx} job {job_id}] {job_name}: "
                 f"{'OK' if ok else f'FAILED ({reason})'}")
    return outcomes


def submit_af3_population(population: List[str], utr_ids: List[str], af3_work_dir: Path,
                           af3_db: str, af3_models: str, af3_module: str, partition: str,
                           time_limit: str, mem: str, cpus: int, gres: str, exclude: str,
                           model_seeds: List[int], generation: int, model_glob: str,
                           batch_size: int, max_concurrent_batches: int, af3_per_gpu: int) -> Path:
    """
    Build one AF3 JSON per individual, skip anything already cached on disk
    (matching sidecar sequence), and submit the rest as SLURM array batches,
    up to `max_concurrent_batches` running concurrently (each batch blocks
    its own thread via `sbatch --wait`; results are logged as each batch
    finishes — no separate polling loop).
    """
    gen_dir = af3_work_dir / f"gen_{generation:04d}"
    jobs_root = gen_dir / "af3_jobs"
    output_root = gen_dir / "af3_output"
    output_root.mkdir(parents=True, exist_ok=True)
    inputs_dir = jobs_root / "inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)

    job_name_paths: List[Tuple[str, Path]] = []
    n_cached = 0
    for sid, seq in zip(utr_ids, population):
        job_name = af3prep.sanitize_name(sid)
        cleaned = af3prep.clean_RNA_seq(seq)
        out_dir = output_root / job_name
        sidecar = inputs_dir / f"{job_name}.seq"
        cached = _find_cached_cif(out_dir, model_glob)
        matches = sidecar.exists() and sidecar.read_text().strip() == cleaned
        if cached is not None and matches:
            n_cached += 1
            continue
        if cached is not None and not matches:
            log.warning(f"    Cached AF3 output for {sid} doesn't match current sequence; re-predicting.")
            n_cached += 1 # only for 1st gen
            continue
        sidecar.write_text(cleaned)
        af3_json = af3prep.build_af3_json(job_name, [{"id": "A", "sequence": cleaned, "type": "rna"}], model_seeds)
        json_path = inputs_dir / f"{job_name}.json"
        json_path.write_text(json.dumps(af3_json, indent=2))
        job_name_paths.append((job_name, json_path))

    if n_cached:
        log.info(f"  Generation {generation}: {n_cached}/{len(utr_ids)} already cached, skipping.")
    if not job_name_paths:
        return output_root

    batches = [job_name_paths[i:i + batch_size] for i in range(0, len(job_name_paths), batch_size)]
    log.info(f"  Generation {generation}: {len(job_name_paths)} to predict in "
             f"{len(batches)} batch(es), {max_concurrent_batches} concurrent, "
             f"{af3_per_gpu} individual(s)/GPU.")

    all_outcomes: List[Tuple[str, bool, str]] = []
    with ThreadPoolExecutor(max_workers=max_concurrent_batches) as ex:
        futures = [
            ex.submit(_run_batch_blocking, idx, batch, jobs_root, output_root, gen_dir,
                      af3_db, af3_models, af3_module, partition, time_limit, mem, cpus,
                      gres, exclude, generation, model_glob, af3_per_gpu)
            for idx, batch in enumerate(batches)
        ]
        for f in as_completed(futures):
            all_outcomes.extend(f.result())

    n_ok = sum(1 for _, ok, _ in all_outcomes if ok)
    log.info(f"  Generation {generation}: AF3 complete — {n_ok}/{len(all_outcomes)} succeeded "
             f"({n_cached} reused from cache).")
    return output_root


def collect_af3_cifs(output_root: Path, utr_ids: List[str], model_glob: str) -> Dict[str, Optional[Path]]:
    result = {}
    for sid in utr_ids:
        out_dir = output_root / af3prep.sanitize_name(sid)
        matches = sorted(out_dir.glob(model_glob))
        result[sid] = matches[0] if matches else None
    return result


def evaluate_af3_rmsd(cif_paths: Dict[str, Optional[Path]], utr_ids: List[str], ref_cifs: List[Path],
                       run_dssr: bool, dssr_path: str, dssr_outdir: Path
                       ) -> Tuple[Dict[str, Optional[float]], Dict[str, Optional[str]]]:
    rmsd_min, dbn_map = {}, {}
    for sid in utr_ids:
        cif = cif_paths.get(sid)
        if cif is None or not Path(cif).exists():
            rmsd_min[sid], dbn_map[sid] = None, None
            continue
        best = None
        for ref in ref_cifs:
            try:
                summary, _, _ = af3cmp.biopython_rmsd(str(ref), str(cif), None, None, RMSD_ATOM)
                if best is None or summary.rmsd < best:
                    best = summary.rmsd
            except Exception as exc:
                log.debug(f"  RMSD failed for {sid} vs {ref}: {exc}")
        rmsd_min[sid] = best
        dbn_map[sid] = None
        if run_dssr:
            try:
                dbn_map[sid] = af3cmp.get_dssr_summary(str(cif), str(dssr_outdir), dssr_path).dbn
            except Exception as exc:
                log.debug(f"  DSSR failed for {sid}: {exc}")
    n_ok = sum(1 for v in rmsd_min.values() if v is not None)
    log.info(f"  AF3/RMSD: {n_ok}/{len(utr_ids)} individuals have a usable structure.")
    return rmsd_min, dbn_map


# ══════════════════════════════════════════════════════════════════════════
# Predicted half-life (reported only, never scored)
# ══════════════════════════════════════════════════════════════════════════

def run_metrics(metrics_script: Path, fasta_5utr: Path, fasta_cds: Path, fasta_3utr: Path,
                 species: str, output_dir: Path) -> bool:
    cmd = [sys.executable, str(metrics_script), '--fasta-5utr', str(fasta_5utr),
           '--fasta-cds', str(fasta_cds), '--fasta-3utr', str(fasta_3utr),
           '--species', species, '--output-dir', str(output_dir), '--force']
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        tsv_dir = output_dir / 'metrics'
        if tsv_dir.exists() and list(tsv_dir.glob('*.tsv')):
            log.warning(f"metrics script exited {result.returncode} but produced TSVs — continuing.")
            return True
        log.error(f"metrics script failed (exit {result.returncode}):\n{result.stderr}")
        return False
    return True


def run_prediction(predict_script: Path, metrics_dir: Path) -> Optional[pd.DataFrame]:
    cmd = [sys.executable, str(predict_script), '--input', str(metrics_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error(f"prediction script failed (exit {result.returncode}): {result.stderr[-2000:]}")
        return None
    out_path = metrics_dir / 'predictions.tsv'
    return pd.read_csv(out_path, sep='\t') if out_path.exists() else None


_FALLBACK_SCHEMAS = {
    'nmd_fragility_full.tsv': [
        'transcript_id', 'gene_id', 'strand', 'model', 'cds_length', 'zone_length',
        'n_transition_fragile_codons', 'n_transversion_fragile_codons', 'n_snv_fragile_codons',
        'n_alt_stop_codons', 'transition_fragile_codon_density', 'transversion_fragile_codon_density',
        'snv_fragile_codon_density', 'alt_stop_codon_density', 'transition_fraction_of_snv_fragile',
    ],
    'junctions.tsv': [
        'transcript_id', 'gene_id', 'n_exons', 'n_junctions', 'strand', 'n_5utr_junctions',
        'n_cds_junctions', 'n_3utr_junctions', 'stop_dist_closest_upstream',
        'stop_dist_closest_downstream', 'stop_dist_last_downstream',
        'start_dist_closest_upstream', 'start_dist_closest_downstream',
    ],
    'architecture.tsv': [
        'transcript_id', 'gene_id', 'strand', 'n_exons', 'n_internal_exons', 'first_exon_length',
        'last_exon_length', 'internal_exon_mean', 'internal_exon_median', 'internal_exon_sd',
        'intron_mean', 'intron_median', 'intron_sd',
    ],
}


def _synthesise_fallback_tsvs(metrics_tsv_dir: Path, utr_ids: List[str], cds_len: int) -> None:
    """Fill default rows for GFF-dependent plugin TSVs so the predictor's join succeeds."""
    for name, cols in _FALLBACK_SCHEMAS.items():
        path = metrics_tsv_dir / name
        if path.exists() and sum(1 for _ in open(path)) > 1:
            continue
        metrics_tsv_dir.mkdir(parents=True, exist_ok=True)
        with open(path, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for sid in utr_ids:
                row = {c: 0 for c in cols}
                row.update(transcript_id=sid, gene_id=sid, strand='+')
                if name == 'nmd_fragility_full.tsv':
                    row.update(model='full', cds_length=cds_len, zone_length=cds_len)
                elif name == 'junctions.tsv':
                    row.update(n_exons=1)
                elif name == 'architecture.tsv':
                    row.update(n_exons=1)
                writer.writerow(row)


def evaluate_halflife(population: List[str], utr_ids: List[str], u5_id: str, u5_seq: str,
                       cds_id: str, cds_seq: str, species: str, metrics_script: Path,
                       predict_script: Path, work_dir: Path) -> Dict[str, Optional[float]]:
    fasta_5, fasta_c, fasta_3 = work_dir / 'pop_5utr.fa', work_dir / 'pop_cds.fa', work_dir / 'pop_3utr.fa'
    write_fasta(fasta_5, [(sid, u5_seq) for sid in utr_ids])
    write_fasta(fasta_c, [(cds_id, cds_seq)])
    write_fasta(fasta_3, list(zip(utr_ids, population)))

    metrics_out = work_dir / 'metrics_run'
    tsv_dir = metrics_out / 'metrics'
    if tsv_dir.exists():
        shutil.rmtree(tsv_dir)
    metrics_out.mkdir(parents=True, exist_ok=True)

    if not run_metrics(metrics_script, fasta_5, fasta_c, fasta_3, species, metrics_out):
        return {sid: None for sid in utr_ids}

    for unused in ('nmd_fragility_core.tsv', 'nmd_fragility_window.tsv'):
        (tsv_dir / unused).unlink(missing_ok=True)
    _synthesise_fallback_tsvs(tsv_dir, utr_ids, len(cds_seq))

    preds = run_prediction(predict_script, tsv_dir)
    if preds is None:
        return {sid: None for sid in utr_ids}
    pred_map = dict(zip(preds['transcript_id'], preds['predicted_halflife']))
    return {sid: pred_map.get(sid) for sid in utr_ids}


# ══════════════════════════════════════════════════════════════════════════
# cmsearch (structural constraint)
# ══════════════════════════════════════════════════════════════════════════

def _parse_tblout(path: str) -> List[dict]:
    hits = []
    if not os.path.exists(path):
        return hits
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.split(maxsplit=17)
            if len(f) < 17:
                continue
            hits.append({"target_name": f[0], "query_name": f[2], "score": f[14], "evalue": f[15]})
    return hits


def run_cmsearch(queries: List[Tuple[str, str]], cm_file: Path, out_dir: Path) -> List[dict]:
    if shutil.which("cmsearch") is None:
        log.warning("cmsearch not on PATH — skipping CM scoring for this generation.")
        return []
    out_dir.mkdir(parents=True, exist_ok=True)
    tblout = out_dir / "cmsearch.tblout"
    with tempfile.TemporaryDirectory() as tmp:
        query_fa = Path(tmp) / "queries.fa"
        write_fasta(query_fa, queries)
        cmd = ["cmsearch", "--notextw", "-o", str(out_dir / "cmsearch.out"),
               "--tblout", str(tblout), "-E", "1000", str(cm_file), str(query_fa)]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            log.warning(f"cmsearch error:\n{r.stderr}")
    return _parse_tblout(str(tblout))


def evaluate_cmscores(cm_model: Path, utr_ids: List[str], population: List[str], work_dir: Path,
                       evalue_threshold: float, generation: int
                       ) -> Tuple[Dict[str, Optional[float]], Dict[str, bool]]:
    hits = run_cmsearch(list(zip(utr_ids, population)), cm_model,
                         work_dir / 'cmsearch_logs' / f"gen_{generation:04d}")
    best_by_target = {}
    for h in hits:
        try:
            ev = float(h['evalue'])
        except (TypeError, ValueError):
            continue
        cur = best_by_target.get(h['target_name'])
        if cur is None or ev < float(cur['evalue']):
            best_by_target[h['target_name']] = h

    cmscores, hit_flags = {}, {}
    for sid in utr_ids:
        h = best_by_target.get(sid)
        if h is not None and float(h['evalue']) <= evalue_threshold:
            cmscores[sid], hit_flags[sid] = float(h['score']), True
        else:
            cmscores[sid], hit_flags[sid] = None, False

    log.info(f"  cmsearch: {sum(hit_flags.values())}/{len(utr_ids)} sequences hit the CM "
             f"(E <= {evalue_threshold})")
    return cmscores, hit_flags


# ══════════════════════════════════════════════════════════════════════════
# Scoring, selection, breeding
# ══════════════════════════════════════════════════════════════════════════

def compute_scores(utr_ids: List[str], rmsd_now, rmsd_parent, cmscore_now, cmscore_parent,
                    hit_flags, max_pids, patent_pid_threshold, no_hit_penalty, patent_pid_penalty
                    ) -> Dict[str, int]:
    scores = {}
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


def scores_to_probabilities(scores: List[int], temperature: float) -> List[float]:
    m = max(scores)
    exps = [math.exp((s - m) / temperature) for s in scores]
    total = sum(exps)
    return [e / total for e in exps]


def breed_next_generation(selected_ids, seq_by_id, rmsd_by_id, cmscore_by_id, population_size,
                           mutation_rate, crossover_rate, rng
                           ) -> Tuple[List[str], List[str], List[Optional[float]], List[Optional[float]]]:
    children, parent_of, parent_rmsd, parent_cm = [], [], [], []
    for _ in range(population_size):
        pa = rng.choice(selected_ids)
        seq = seq_by_id[pa]
        if crossover_rate > 0 and rng.random() < crossover_rate and len(selected_ids) > 1:
            pb = rng.choice(selected_ids)
            seq = uniform_crossover(seq, seq_by_id[pb], rng)
        children.append(point_mutate(seq, mutation_rate, rng))
        parent_of.append(pa)
        parent_rmsd.append(rmsd_by_id.get(pa))
        parent_cm.append(cmscore_by_id.get(pa))
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
    gens = [r['generation'] for r in best_per_gen]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), dpi=150, sharex=True)
    ax1.plot(gens, [r['best_score'] for r in best_per_gen], label='Best', marker='o', markersize=3)
    ax1.plot(gens, [r['mean_score'] for r in best_per_gen], label='Mean', marker='o', markersize=3)
    ax1.plot(gens, [r['worst_score'] for r in best_per_gen], label='Worst', linestyle='--', alpha=0.7)
    ax1.set_ylabel('Point score'); ax1.legend(); ax1.grid(alpha=0.3)
    ax1.set_title("3' UTR GA — score over generations")

    ax2.plot(gens, [r['mean_rmsd'] or float('nan') for r in best_per_gen], color='#762a83', marker='o', markersize=3)
    ax2.set_ylabel('Mean min-RMSD (Å)', color='#762a83')
    ax2b = ax2.twinx()
    ax2b.plot(gens, [r['mean_cmscore'] or float('nan') for r in best_per_gen], color='#e08214', marker='s', markersize=3)
    ax2b.set_ylabel('Mean cmscore (bits)', color='#e08214')
    ax2.set_xlabel('Generation'); ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    log.info(f"Wrote {out_path}")


# ══════════════════════════════════════════════════════════════════════════
# Main GA loop
# ══════════════════════════════════════════════════════════════════════════

def run_ga(seed_3utr_seq: str, u5_id: str, u5_seq: str, cds_id: str, cds_seq: str, species: str,
           metrics_script: Path, predict_script: Path, ref_cifs: List[Path], cm_model: Path,
           cm_evalue_threshold: float, no_hit_penalty: int, patent_seqs: List[dict],
           patent_pid_threshold: float, patent_pid_penalty: int, af3_work_dir: Path, af3_db: str,
           af3_models: str, af3_module: str, af3_partition: str, af3_time: str, af3_mem: str,
           af3_cpus: int, af3_gres: str, af3_exclude: str, af3_model_seeds: List[int],
           af3_model_glob: str, af3_batch_size: int, af3_max_concurrent_batches: int,
           af3_per_gpu: int, run_dssr: bool, dssr_path: str, output_dir: Path, population_size: int, n_select: int,
           generations: int, mutation_rate: float, crossover_rate: float, temperature: float,
           rng_seed: int, make_plot: bool) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    results_dir = output_dir / 'results'
    results_dir.mkdir(exist_ok=True)
    work_dir = output_dir / '_workspace'
    work_dir.mkdir(exist_ok=True)
    dssr_outdir = work_dir / 'dssr_logs'

    log.info("═" * 60)
    log.info("3' UTR Genetic Algorithm — point-scored, AF3-in-the-loop")
    log.info(f"  Population {population_size} | n-select {n_select} | generations {generations}")
    log.info(f"  Mutation rate {mutation_rate} (flat) | crossover rate {crossover_rate} | temp {temperature}")
    log.info(f"  cmsearch E<= {cm_evalue_threshold} | no-hit penalty {no_hit_penalty}")
    log.info(f"  Patent motifs {len(patent_seqs)} | PID> {patent_pid_threshold}% -> {patent_pid_penalty}")
    log.info(f"  AF3 batch size {af3_batch_size} | max concurrent batches {af3_max_concurrent_batches} | "
             f"{af3_per_gpu} individual(s)/GPU")
    log.info(f"  Output: {output_dir}")

    ckpt_path = find_latest_checkpoint(output_dir)
    if ckpt_path:
        log.info("═" * 60 + "\nRESUMING FROM CHECKPOINT")
        ckpt = load_checkpoint(ckpt_path)
        rng = random.Random(rng_seed)
        rng.setstate(ckpt.rng_state)
        population, parent_of = ckpt.population, ckpt.parent_of
        parent_rmsd, parent_cmscore = ckpt.parent_rmsd, ckpt.parent_cmscore
        seed_rmsd, seed_cmscore = ckpt.seed_rmsd, ckpt.seed_cmscore
        all_rows, selected_rows, best_per_gen = ckpt.all_rows, ckpt.selected_rows, ckpt.best_per_gen
        all_time_best_score = ckpt.all_time_best_score
        all_time_best_seq, all_time_best_id, all_time_best_gen = (
            ckpt.all_time_best_seq, ckpt.all_time_best_id, ckpt.all_time_best_gen)
        start_generation = ckpt.generation + 1
        log.info(f"  Resuming at generation {start_generation}; best score so far {all_time_best_score}")

        if start_generation <= generations:
            # ckpt.population was scored to produce ckpt.generation's results but never
            # bred forward — redo that breeding step now with the restored RNG so the
            # result is identical to an uninterrupted run.
            gen_selected = sorted((r for r in selected_rows if r['generation'] == ckpt.generation),
                                   key=lambda r: r['rank'])
            selected_ids = [r['sample_id'] for r in gen_selected]
            gen_rows = [r for r in all_rows if r['generation'] == ckpt.generation]
            seq_by_id = {r['sample_id']: r['sequence'] for r in gen_rows}
            rmsd_by_id = {r['sample_id']: r['rmsd_min'] for r in gen_rows}
            cm_by_id = {r['sample_id']: r['cmscore'] for r in gen_rows}
            if not selected_ids:
                raise RuntimeError(f"Checkpoint gen {ckpt.generation} has no selected_rows to breed from.")
            population, parent_of, parent_rmsd, parent_cmscore = breed_next_generation(
                selected_ids, seq_by_id, rmsd_by_id, cm_by_id, population_size,
                mutation_rate, crossover_rate, rng)
    else:
        log.info("═" * 60 + "\nNO CHECKPOINT — starting fresh")
        rng = random.Random(rng_seed)
        start_generation = 1
        seed_id = "seed_3utr"

        log.info("Computing baseline (seed) metrics for generation-1 comparisons...")
        out_root = submit_af3_population(
            [seed_3utr_seq], [seed_id], af3_work_dir, af3_db, af3_models, af3_module,
            af3_partition, af3_time, af3_mem, af3_cpus, af3_gres, af3_exclude,
            af3_model_seeds, 0, af3_model_glob, af3_batch_size, af3_max_concurrent_batches, af3_per_gpu)
        seed_cif = collect_af3_cifs(out_root, [seed_id], af3_model_glob)
        seed_rmsd_map, _ = evaluate_af3_rmsd(seed_cif, [seed_id], ref_cifs, run_dssr, dssr_path, dssr_outdir)
        seed_cm_map, seed_hit_map = evaluate_cmscores(cm_model, [seed_id], [seed_3utr_seq], work_dir,
                                                        cm_evalue_threshold, 0)
        seed_rmsd, seed_cmscore = seed_rmsd_map.get(seed_id), seed_cm_map.get(seed_id)
        log.info(f"  Seed baseline: RMSD={seed_rmsd}, cmscore={seed_cmscore}, cm_hit={seed_hit_map.get(seed_id)}")

        population = [point_mutate(seed_3utr_seq, mutation_rate, rng) for _ in range(population_size)]
        parent_of = [seed_id] * population_size
        parent_rmsd = [seed_rmsd] * population_size
        parent_cmscore = [seed_cmscore] * population_size
        all_rows, selected_rows, best_per_gen = [], [], []
        all_time_best_score = None
        all_time_best_seq, all_time_best_id, all_time_best_gen = seed_3utr_seq, seed_id, 0

    log.info("═" * 60)

    for gen in range(start_generation, generations + 1):
        gen_start = time.time()
        log.info(f"\n── Generation {gen}/{generations} (mutation rate = {mutation_rate}) ──")

        utr_ids = [f"ind_{gen:04d}_{i:05d}" for i in range(population_size)]
        seq_by_id = dict(zip(utr_ids, population))
        parent_of_by_id = dict(zip(utr_ids, parent_of))
        parent_rmsd_by_id = dict(zip(utr_ids, parent_rmsd))
        parent_cm_by_id = dict(zip(utr_ids, parent_cmscore))

        out_root = submit_af3_population(
            population, utr_ids, af3_work_dir, af3_db, af3_models, af3_module,
            af3_partition, af3_time, af3_mem, af3_cpus, af3_gres, af3_exclude,
            af3_model_seeds, gen, af3_model_glob, af3_batch_size, af3_max_concurrent_batches, af3_per_gpu)
        cif_paths = collect_af3_cifs(out_root, utr_ids, af3_model_glob)
        rmsd_now, dbn_now = evaluate_af3_rmsd(cif_paths, utr_ids, ref_cifs, run_dssr, dssr_path, dssr_outdir)
        cmscore_now, hit_flags = evaluate_cmscores(cm_model, utr_ids, population, work_dir,
                                                     cm_evalue_threshold, gen)

        max_pids = {}
        for sid, seq in zip(utr_ids, population):
            pids = check_patent_pid(seq, patent_seqs) if patent_seqs else []
            max_pids[sid] = max(pids) if pids else None

        halflife_map = evaluate_halflife(population, utr_ids, u5_id, u5_seq, cds_id, cds_seq,
                                          species, metrics_script, predict_script, work_dir)

        scores = compute_scores(utr_ids, rmsd_now, parent_rmsd_by_id, cmscore_now, parent_cm_by_id,
                                 hit_flags, max_pids, patent_pid_threshold, no_hit_penalty, patent_pid_penalty)
        score_list = [scores[sid] for sid in utr_ids]
        probabilities = scores_to_probabilities(score_list, temperature)
        prob_by_id = dict(zip(utr_ids, probabilities))

        for sid in utr_ids:
            all_rows.append({
                'generation': gen, 'sample_id': sid, 'parent_id': parent_of_by_id[sid],
                'score': scores[sid], 'probability': prob_by_id[sid],
                'rmsd_min': rmsd_now.get(sid), 'parent_rmsd': parent_rmsd_by_id[sid],
                'cmscore': cmscore_now.get(sid), 'parent_cmscore': parent_cm_by_id[sid],
                'cm_hit': hit_flags.get(sid, False), 'max_patent_pid': max_pids.get(sid),
                'predicted_halflife': halflife_map.get(sid), 'dssr_dbn': dbn_now.get(sid),
                'sequence': seq_by_id[sid], 'length': len(seq_by_id[sid]),
                'mutations_from_seed': sum(a != b for a, b in zip(seq_by_id[sid], seed_3utr_seq)
                                            if len(seq_by_id[sid]) == len(seed_3utr_seq)),
            })

        valid_rmsd = [v for v in rmsd_now.values() if v is not None]
        valid_cm = [v for v in cmscore_now.values() if v is not None]
        valid_hl = [v for v in halflife_map.values() if v is not None]
        best_sid = max(utr_ids, key=lambda s: scores[s])

        log.info(f"  Best score this gen: {scores[best_sid]} ({best_sid}) | "
                 f"mean {sum(score_list)/len(score_list):.2f} | worst {min(score_list)}")
        log.info(f"  CM hits: {sum(hit_flags.values())}/{population_size} | "
                 f"elapsed {time.time() - gen_start:.1f}s")

        best_per_gen.append({
            'generation': gen, 'best_score': max(score_list), 'mean_score': sum(score_list) / len(score_list),
            'worst_score': min(score_list),
            'mean_rmsd': (sum(valid_rmsd) / len(valid_rmsd)) if valid_rmsd else None,
            'mean_cmscore': (sum(valid_cm) / len(valid_cm)) if valid_cm else None,
            'mean_halflife': (sum(valid_hl) / len(valid_hl)) if valid_hl else None,
            'best_id': best_sid, 'best_sequence': seq_by_id[best_sid],
        })

        if all_time_best_score is None or scores[best_sid] > all_time_best_score:
            all_time_best_score = scores[best_sid]
            all_time_best_seq, all_time_best_id, all_time_best_gen = seq_by_id[best_sid], best_sid, gen
            log.info(f"  \u2605 New all-time best score: {all_time_best_score}")

        selected_ids = rng.choices(utr_ids, weights=probabilities, k=n_select)
        for rank, sid in enumerate(selected_ids):
            selected_rows.append({'generation': gen, 'rank': rank + 1, 'sample_id': sid,
                                   'score': scores[sid], 'probability': prob_by_id[sid],
                                   'sequence': seq_by_id[sid]})

        save_checkpoint(GACheckpoint(
            generation=gen, population=population, parent_of=parent_of, parent_rmsd=parent_rmsd,
            parent_cmscore=parent_cmscore, seed_rmsd=seed_rmsd, seed_cmscore=seed_cmscore,
            all_rows=all_rows, selected_rows=selected_rows, best_per_gen=best_per_gen,
            all_time_best_score=all_time_best_score, all_time_best_seq=all_time_best_seq,
            all_time_best_id=all_time_best_id, all_time_best_gen=all_time_best_gen,
            rng_state=rng.getstate()), output_dir)

        if gen == generations:
            break
        population, parent_of, parent_rmsd, parent_cmscore = breed_next_generation(
            selected_ids, seq_by_id, rmsd_now, cmscore_now, population_size,
            mutation_rate, crossover_rate, rng)

    # ── Outputs ──────────────────────────────────────────────────────────
    def _write_tsv(path, rows):
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(rows)
        log.info(f"Wrote {path}")

    _write_tsv(results_dir / 'best_per_generation.tsv', best_per_gen)
    _write_tsv(results_dir / 'all_generations.tsv', all_rows)
    _write_tsv(results_dir / 'selected_per_generation.tsv', selected_rows)

    write_fasta(results_dir / 'best_3utr.fa', [('best_3utr_evolved', all_time_best_seq)])
    write_fasta(results_dir / 'seed_3utr.fa', [('seed_3utr_original', seed_3utr_seq)])
    log.info(f"Wrote {results_dir / 'best_3utr.fa'}")

    if make_plot:
        plot_scores_over_generations(best_per_gen, results_dir / 'score_over_generations.png')

    with open(results_dir / 'evolution_summary.txt', 'w') as fh:
        fh.write("3' UTR Genetic Algorithm — Evolution Summary\n" + "=" * 60 + "\n\n")
        fh.write(f"Generations: {generations} | Population: {population_size} | n-select: {n_select}\n")
        fh.write(f"Mutation rate: {mutation_rate} | Crossover rate: {crossover_rate} | Temp: {temperature}\n")
        fh.write(f"RNG seed: {rng_seed}\nRef CIFs: {[str(c) for c in ref_cifs]}\nCM model: {cm_model}\n")
        fh.write(f"cmsearch E<= {cm_evalue_threshold} | No-hit penalty: {no_hit_penalty}\n")
        fh.write(f"Patent motifs: {len(patent_seqs)} | PID> {patent_pid_threshold}% -> {patent_pid_penalty}\n\n")
        fh.write(f"Seed length: {len(seed_3utr_seq)} nt | Seed RMSD/cmscore: {seed_rmsd} / {seed_cmscore}\n")
        fh.write(f"Best length: {len(all_time_best_seq)} nt\n")
        fh.write(f"All-time best score: {all_time_best_score} (id={all_time_best_id}, gen={all_time_best_gen})\n\n")
        for row in best_per_gen:
            fh.write(f"  Gen {row['generation']:>3d}: best {row['best_score']} "
                     f"(mean {row['mean_score']:.2f}, worst {row['worst_score']}) | "
                     f"rmsd {row['mean_rmsd']} | cmscore {row['mean_cmscore']} | "
                     f"halflife {row['mean_halflife']}\n")
        fh.write("\nBest evolved 3' UTR sequence:\n")
        for i in range(0, len(all_time_best_seq), 60):
            fh.write(all_time_best_seq[i:i + 60] + '\n')
    log.info(f"Wrote {results_dir / 'evolution_summary.txt'}")

    log.info("\n" + "═" * 60)
    log.info(f"GA complete. All-time best score: {all_time_best_score} "
             f"(id={all_time_best_id}, gen={all_time_best_gen})")
    log.info(f"Results in: {results_dir}\n" + "═" * 60)


# ══════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════

def main():
    global af3prep, af3cmp
    af3prep = _import_sibling('af3prep', 'prepare_af3_jobs.py')
    af3cmp = _import_sibling('af3cmp', 'compare_af3_candidates.py')

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('--fasta-5utr', required=True, dest='fasta_5utr')
    p.add_argument('--fasta-cds', required=True, dest='fasta_cds')
    p.add_argument('--fasta-3utr', required=True, dest='fasta_3utr')
    p.add_argument('--species', required=True)
    p.add_argument('--metrics-script', required=True, dest='metrics_script')
    p.add_argument('--predict-script', required=True, dest='predict_script')

    p.add_argument('--ref-cif', required=True, action='append', dest='ref_cifs',
                    help="Repeat exactly 3 times; RMSD = min over all three.")
    p.add_argument('--skip-dssr', action='store_true')
    p.add_argument('--dssr-path', default='x3dna-dssr')
    p.add_argument('--cm-model', required=True, dest='cm_model')
    p.add_argument('--cm-evalue-threshold', type=float, default=0.01, dest='cm_evalue_threshold')
    p.add_argument('--no-hit-penalty', type=int, default=-10, dest='no_hit_penalty')

    p.add_argument('--patent-pid-threshold', type=float, default=80.0, dest='patent_pid_threshold')
    p.add_argument('--patent-pid-penalty', type=int, default=-10, dest='patent_pid_penalty')
    p.add_argument('--no-default-patents', action='store_true', dest='no_default_patents')

    p.add_argument('--af3-work-dir', required=True, dest='af3_work_dir')
    p.add_argument('--af3-db', default='/opt/alphafold3/databases', dest='af3_db')
    p.add_argument('--af3-models', required=True, dest='af3_models')
    p.add_argument('--af3-module', default='alphafold3/3.0.3', dest='af3_module')
    p.add_argument('--af3-model-seed', type=int, action='append', dest='af3_model_seeds')
    p.add_argument('--af3-model-glob', default='**/*_model.cif', dest='af3_model_glob',
                    help="Glob (relative to each job's output dir) for the top model .cif.")
    p.add_argument('--af3-partition', default='aoraki_gpu_L40', dest='af3_partition')
    p.add_argument('--af3-time', default='02:00:00', dest='af3_time')
    p.add_argument('--af3-mem', default='32G', dest='af3_mem')
    p.add_argument('--af3-cpus-per-task', type=int, default=8, dest='af3_cpus')
    p.add_argument('--af3-gres', default='gpu:1', dest='af3_gres')
    p.add_argument('--af3-exclude', default='', dest='af3_exclude')
    p.add_argument('--af3-batch-size', type=int, default=100, dest='af3_batch_size',
                    help="Individuals per SLURM array batch (default: 100).")
    p.add_argument('--af3-max-concurrent-batches', type=int, default=10, dest='af3_max_concurrent_batches',
                    help="Batches submitted/blocking at once (default: 10). Tune both of "
                         "these to your cluster's per-user job/task QOS limits.")
    p.add_argument('--af3-per-gpu', type=int, default=1, dest='af3_per_gpu',
                    help="Individuals grouped onto one GPU allocation, run concurrently "
                         "in the background on it (default: 1). Raise this for small "
                         "inputs that don't use a full GPU on their own — check "
                         "nvidia-smi VRAM usage for one prediction first, and cuts the "
                         "number of SLURM array tasks/GPU-node waits needed too.")

    p.add_argument('--population', type=int, default=1000)
    p.add_argument('--n-select', type=int, default=10, dest='n_select')
    p.add_argument('--generations', type=int, default=30)
    p.add_argument('--mutation-rate', type=float, default=0.02, dest='mutation_rate')
    p.add_argument('--crossover-rate', type=float, default=0.0, dest='crossover_rate')
    p.add_argument('--temperature', type=float, default=1.0)
    p.add_argument('--seed', type=int, default=42)

    p.add_argument('--output-dir', '-o', default='ga_output', dest='output_dir')
    p.add_argument('--no-plot', action='store_true', dest='no_plot')
    p.add_argument('-v', '--verbose', action='store_true')

    args = p.parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if len(args.ref_cifs) != 3:
        log.error(f"--ref-cif must be given exactly 3 times (got {len(args.ref_cifs)}).")
        sys.exit(1)

    required = [('--fasta-5utr', args.fasta_5utr), ('--fasta-cds', args.fasta_cds),
                ('--fasta-3utr', args.fasta_3utr), ('--metrics-script', args.metrics_script),
                ('--predict-script', args.predict_script), ('--cm-model', args.cm_model),
                ] + [('--ref-cif', c) for c in args.ref_cifs]
    for label, path_str in required:
        if not Path(path_str).exists():
            log.error(f"{label} not found: {path_str}")
            sys.exit(1)

    if shutil.which("cmsearch") is None:
        log.error("'cmsearch' (Infernal) not found on PATH.")
        sys.exit(1)
    if shutil.which("sbatch") is None:
        log.error("'sbatch' not found on PATH — required to submit AF3 SLURM array jobs.")
        sys.exit(1)
    if not args.skip_dssr and shutil.which(args.dssr_path) is None:
        log.warning(f"'{args.dssr_path}' not found — DSSR reporting skipped (scoring unaffected).")
        args.skip_dssr = True
    if args.af3_batch_size < 1 or args.af3_max_concurrent_batches < 1 or args.af3_per_gpu < 1:
        log.error("--af3-batch-size, --af3-max-concurrent-batches, and --af3-per-gpu must all be >= 1.")
        sys.exit(1)

    try:
        u5_id, u5_seq = read_single_fasta(Path(args.fasta_5utr))
        cds_id, cds_seq = read_single_fasta(Path(args.fasta_cds))
        u3_id, u3_seq = read_single_fasta(Path(args.fasta_3utr))
    except ValueError as exc:
        log.error(str(exc))
        sys.exit(1)

    log.info(f"5' UTR: '{u5_id}' {len(u5_seq)} nt (fixed) | CDS: '{cds_id}' {len(cds_seq)} nt (fixed) | "
             f"3' UTR: '{u3_id}' {len(u3_seq)} nt (seed)")

    run_ga(
        seed_3utr_seq=u3_seq, u5_id=u5_id, u5_seq=u5_seq, cds_id=cds_id, cds_seq=cds_seq,
        species=args.species, metrics_script=Path(args.metrics_script), predict_script=Path(args.predict_script),
        ref_cifs=[Path(c) for c in args.ref_cifs], cm_model=Path(args.cm_model),
        cm_evalue_threshold=args.cm_evalue_threshold, no_hit_penalty=args.no_hit_penalty,
        patent_seqs=[] if args.no_default_patents else PATENT_SEQUENCES,
        patent_pid_threshold=args.patent_pid_threshold, patent_pid_penalty=args.patent_pid_penalty,
        af3_work_dir=Path(args.af3_work_dir), af3_db=args.af3_db, af3_models=args.af3_models,
        af3_module=args.af3_module, af3_partition=args.af3_partition, af3_time=args.af3_time,
        af3_mem=args.af3_mem, af3_cpus=args.af3_cpus, af3_gres=args.af3_gres, af3_exclude=args.af3_exclude,
        af3_model_seeds=args.af3_model_seeds or [1], af3_model_glob=args.af3_model_glob,
        af3_batch_size=args.af3_batch_size, af3_max_concurrent_batches=args.af3_max_concurrent_batches,
        af3_per_gpu=args.af3_per_gpu,
        run_dssr=not args.skip_dssr, dssr_path=args.dssr_path, output_dir=Path(args.output_dir),
        population_size=args.population, n_select=args.n_select, generations=args.generations,
        mutation_rate=args.mutation_rate, crossover_rate=args.crossover_rate, temperature=args.temperature,
        rng_seed=args.seed, make_plot=not args.no_plot,
    )


if __name__ == '__main__':
    main()