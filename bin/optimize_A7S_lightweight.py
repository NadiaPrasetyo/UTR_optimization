#!/usr/bin/env python3
"""
mutate_3utr_ga.py
─────────────────────────────────────────────────────────────────────────────
Genetic Algorithm to evolve 3' UTR sequences that maximise predicted mRNA
half-life while preserving structural homology to the original 3' UTR,
using the 01b_metrics.py + LightGBM prediction pipeline plus an Infernal
covariance-model (CM) search against the original secondary structure.

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
    [--mutation-rate-start 0.05]       \\   # per-nt mutation prob, generation 1
    [--mutation-rate-end   0.005]      \\   # per-nt mutation prob, final gen
    [--cooling-schedule exponential]   \\   # 'exponential' (SA-style) or 'linear'
    [--cooling-rate 3.0]               \\   # decay constant for exponential cooling
    [--crossover-rate 0.5]             \\   # crossover probability
    [--elite-frac  0.1]                \\   # fraction kept as elites (no mutation)
    [--tournament-k 3]                 \\   # tournament selection size
    [--cm-model    original_3utr.cm]   \\   # calibrated Infernal CM of the seed structure
    [--cm-evalue-threshold 0.01]       \\   # E-value cutoff defining a "hit"
    [--halflife-weight 0.5]            \\   # weight of half-life in combined fitness
    [--cm-weight 0.5]                  \\   # weight of cmscore in combined fitness
    [--no-hit-penalty -1000.0]         \\   # fitness assigned when there is no CM hit
    [--use-default-patents]            \\   # enable A7 patent filtering
    [--patent-pid-threshold 80.0]      \\   # PID threshold (%) for patent filtering
    [--patent-weight 0.0]              \\   # weight of low-PID reward in combined fitness
    [--patent-pid-penalty -1000.0]     \\   # fitness assigned when PID >= threshold
    [--seed        42]                 \\
    [--no-plot]                        \\   # skip generating the fitness plot

WHAT IT DOES
────────────
1.  Reads the fixed 5' UTR and CDS (these are never mutated).
2.  Seeds a population of N copies of the initial 3' UTR, then applies
    random point mutations to every non-elite member each generation, with
    the per-nucleotide mutation rate annealed (simulated-annealing style)
    from --mutation-rate-start down to --mutation-rate-end over the run.
3.  Each generation:
      a.  Writes three temporary FASTA files (5UTR, CDS, 3UTR variants).
      b.  Calls 01b_metrics.py --fasta-* to compute feature TSVs.
      c.  Calls the LightGBM prediction script to get predicted half-life.
      d.  If --cm-model is given, runs cmsearch against the covariance
          model built from the original 3' UTR structure. Sequences with
          no hit above --cm-evalue-threshold receive --no-hit-penalty as
          their fitness, regardless of predicted half-life.
      e.  Combines (normalised) predicted half-life and (normalised)
          cmscore into a single fitness value using --halflife-weight and
          --cm-weight, and ranks/selects on that combined fitness.
      f.  Keeps elites, breeds next generation via tournament selection +
          uniform crossover + annealed point mutation.
4.  Writes results/  containing:
      - best_per_generation.tsv      — top fitness/half-life/cmscore per gen
      - all_generations.tsv          — every evaluated individual
      - best_3utr.fa                 — FASTA of the all-time best sequence
      - evolution_summary.txt        — human-readable run summary
      - fitness_over_generations.png — plot of fitness (and components) per gen

NOTE ON THE CM MODEL
─────────────────────
--cm-model expects a *calibrated* Infernal covariance model (the output of
`cmbuild` + `cmcalibrate` run on a Stockholm alignment/secondary-structure
annotation of the original 3' UTR). Building/calibrating that model is
outside the scope of this script — it is assumed to already exist. If
--cm-model is omitted, the structural constraint is skipped entirely and
fitness reduces to (normalised) predicted half-life alone.
─────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import argparse
import csv
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

# ── Logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('ga_3utr')

NUCLEOTIDES = ['A', 'T', 'G', 'C']

# ── Patent sequences (from US/international patent applications) ──────────
# Used for PID (percent identity) filtering to avoid re-patenting existing motifs
PATENT_SEQUENCES = [
    {"id": "A7_30nt", "start": 9651, "end": 9680, "seq": "CAGACCCTGGTCCGGGGCAATGGGACCACT"},
    {"id": "A7_32nt", "start": 9650, "end": 9681, "seq": "TCAGACCCTGGTCCGGGGCAAATGGGACCACTG"},
    {"id": "A7_34nt", "start": 9649, "end": 9682, "seq": "GTCAGACCCTGGTCCGGGGCAATGGGACCACTGT"},
    {"id": "A7_36nt", "start": 9648, "end": 9683, "seq": "GGTCAGACCCTGGTCCGGGGCAAATGGGACCACTGTT"},
    {"id": "A7_40nt", "start": 9646, "end": 9685, "seq": "TGGGTCAGACCCTGGTCCGGGGCAAATGGGACCACTGTTTC"},
    {"id": "A7_43nt", "start": 9646, "end": 9688, "seq": "TGGGTCAGACCCTGGTCCGGGGCAATGGGACCACTGTTTCGCG"},
    {"id": "A7_47nt", "start": 9646, "end": 9692, "seq": "TGGGTCAGACCCTGGTCCGGGGCAATGGGACCACTGTTTCGCGTTTA"},
]


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
# Simulated-annealing mutation-rate schedule
# ══════════════════════════════════════════════════════════════════════════════

def compute_mutation_rate(
    generation: int,
    total_generations: int,
    rate_start: float,
    rate_end: float,
    schedule: str = 'exponential',
    cooling_rate: float = 3.0,
) -> float:
    """
    Return the per-nucleotide mutation rate for *generation* (1-indexed),
    annealed from rate_start (generation 1) down to rate_end (final
    generation).

    schedule='exponential' follows a classic simulated-annealing cooling
    curve:  rate(frac) = rate_end + (rate_start - rate_end) * exp(-k * frac)
    where frac goes from 0 (first generation) to 1 (last generation) and
    k = cooling_rate controls how quickly the rate drops early on.

    schedule='linear' interpolates rate_start → rate_end evenly across
    generations.
    """
    if total_generations <= 1:
        return rate_end

    frac = (generation - 1) / (total_generations - 1)
    frac = min(max(frac, 0.0), 1.0)

    if schedule == 'linear':
        return rate_start + (rate_end - rate_start) * frac

    # exponential / simulated-annealing cooling
    decay = math.exp(-cooling_rate * frac)
    return rate_end + (rate_start - rate_end) * decay


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


def tournament_select(population: List[str], fitness: List[float],
                       k: int, rng: random.Random) -> str:
    """Select one individual via k-way tournament (highest fitness wins)."""
    contestants = rng.sample(range(len(population)), min(k, len(population)))
    winner = max(contestants, key=lambda i: fitness[i])
    return population[winner]


# ══════════════════════════════════════════════════════════════════════════════
# Patent sequence identity (PID) checking via Smith-Waterman alignment
# ══════════════════════════════════════════════════════════════════════════════

def _sw_pure(seq_a: str, seq_b: str) -> dict:
    """
    Pure-Python Smith-Waterman local alignment.
    Returns {'query': ..., 'subject': ..., 'identities': int, 'alignment_length': int}
    for the highest-scoring local alignment found.
    """
    match_score = 2
    mismatch_score = -1
    gap_score = -1

    # DP matrix
    m, n = len(seq_a), len(seq_b)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0
    max_i, max_j = 0, 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq_a[i - 1] == seq_b[j - 1]:
                score = H[i - 1][j - 1] + match_score
            else:
                score = H[i - 1][j - 1] + mismatch_score
            score = max(score, H[i - 1][j] + gap_score, H[i][j - 1] + gap_score, 0)
            H[i][j] = score
            if score > max_score:
                max_score = score
                max_i, max_j = i, j

    # Traceback from (max_i, max_j) back to start of local alignment
    align_a, align_b = [], []
    i, j = max_i, max_j
    while i > 0 and j > 0 and H[i][j] > 0:
        if seq_a[i - 1] == seq_b[j - 1]:
            align_a.append(seq_a[i - 1])
            align_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif H[i - 1][j - 1] >= H[i - 1][j] and H[i - 1][j - 1] >= H[i][j - 1]:
            align_a.append(seq_a[i - 1])
            align_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif H[i - 1][j] > H[i][j - 1]:
            align_a.append(seq_a[i - 1])
            align_b.append('-')
            i -= 1
        else:
            align_a.append('-')
            align_b.append(seq_b[j - 1])
            j -= 1

    align_a = ''.join(reversed(align_a))
    align_b = ''.join(reversed(align_b))
    identities = sum(a == b and a != '-' for a, b in zip(align_a, align_b))
    alignment_length = max(len(align_a), len(align_b))

    return {
        'query': align_a,
        'subject': align_b,
        'identities': identities,
        'alignment_length': alignment_length,
    }


def _sw_biopython(seq_a: str, seq_b: str) -> Optional[dict]:
    """
    Smith-Waterman alignment using BioPython's pairwise2 module.
    Returns a dict with 'query', 'subject', 'identities', 'alignment_length'
    or None if BioPython is not available.
    """
    try:
        from Bio import pairwise2
        from Bio.Seq import Seq
    except ImportError:
        return None

    seq_a_bio = Seq(seq_a)
    seq_b_bio = Seq(seq_b)

    # Use BLOSUM62-like scoring adapted for nucleotides
    match_score = 2
    mismatch_score = -1
    gap_open = -1
    gap_extend = -1

    alignments = pairwise2.align.localms(
        seq_a_bio, seq_b_bio,
        match_score, mismatch_score,
        gap_open, gap_extend,
    )
    if not alignments:
        return None

    best = alignments[0]
    align_a = str(best[0])
    align_b = str(best[1])
    identities = sum(a == b and a != '-' for a, b in zip(align_a, align_b))
    alignment_length = max(len(align_a), len(align_b))

    return {
        'query': align_a,
        'subject': align_b,
        'identities': identities,
        'alignment_length': alignment_length,
    }


def smith_waterman(seq_a: str, seq_b: str) -> dict:
    """
    Smith-Waterman local alignment; tries BioPython first, falls back to pure Python.
    Returns {'query': ..., 'subject': ..., 'identities': int, 'alignment_length': int}.
    """
    try:
        result = _sw_biopython(seq_a, seq_b)
        if result:
            return result
    except Exception:
        pass
    return _sw_pure(seq_a, seq_b)


def calculate_pid(seq_query: str, seq_subject: str) -> float:
    """
    Calculate percent identity (PID) between two sequences using Smith-Waterman.
    PID = (identities / length_of_shorter_sequence) * 100
    """
    alignment = smith_waterman(seq_query, seq_subject)
    shorter_len = min(len(seq_query), len(seq_subject))
    if shorter_len == 0:
        return 0.0
    pid = (alignment['identities'] / shorter_len) * 100.0
    return pid


def check_patent_pid(
    sequence: str,
    patent_seqs: List[dict],
    pid_threshold: float = 80.0,
) -> Tuple[List[float], bool]:
    """
    Check the query sequence against all patent sequences.

    Returns (pid_list, passes_filter) where:
      - pid_list: list of PID values (%) aligned with patent_seqs
      - passes_filter: True if the sequence has < pid_threshold PID to *all* patents
                       False if it has >= pid_threshold PID to *any* patent
    """
    pids = []
    for patent in patent_seqs:
        pid = calculate_pid(sequence, patent['seq'])
        pids.append(pid)

    # Passes filter if ALL PIDs are below threshold
    passes = all(pid < pid_threshold for pid in pids)
    return pids, passes


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
# Structural constraint: cmsearch against the original 3' UTR covariance model
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class CmsearchResult:
    query_name: str          # CM model name (from the .cm file)
    target_name: str         # sequence id from the query FASTA (a GA individual)
    cm_sig_name: Optional[str]
    target_accession: str
    score: str
    evalue: str
    bias: str
    strand: str
    seq_from: str
    seq_to: str
    gc: str
    description: str


def _parse_tblout(tblout_file: str) -> List[dict]:
    """
    Parse an Infernal 1.1.x --tblout file into a list of dicts.

    Column layout (whitespace-delimited):
      0 target name, 1 target accession, 2 query name, 3 query accession,
      4 mdl, 5 mdl from, 6 mdl to, 7 seq from, 8 seq to, 9 strand,
      10 trunc, 11 pass, 12 gc, 13 bias, 14 score, 15 E-value, 16 inc,
      17+ description of target
    """
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
                "target_name":      fields[0],
                "target_accession": fields[1],
                "query_name":       fields[2],
                "query_accession":  fields[3],
                "mdl":              fields[4],
                "mdl_from":         fields[5],
                "mdl_to":           fields[6],
                "seq_from":         fields[7],
                "seq_to":           fields[8],
                "strand":           fields[9],
                "trunc":            fields[10],
                "pass":             fields[11],
                "gc":               fields[12],
                "bias":             fields[13],
                "score":            fields[14],
                "evalue":           fields[15],
                "inc":              fields[16],
                "description":      fields[17].rstrip('\n') if len(fields) > 17 else '',
            })
    return hits


def run_cmsearch(queries, cm_file, output_dir="."):
    """
    Run cmsearch for each query (a GA individual's 3' UTR) against the
    original-structure covariance model.

    Parameters
    ----------
    queries    : list of (name, seq) tuples
    cm_file    : path to Infernal .cm file
    output_dir : directory where persistent output files are written

    Persistent output files (in output_dir):
      patent_cmsearch.out     — full human-readable cmsearch output  (-o)
      patent_cmsearch.tblout  — tabular hits                         (--tblout)
      patent_cmsearch.sto     — Stockholm alignment of all hits      (-A)
    """
    if shutil.which("cmsearch") is None:
        print("Warning: cmsearch not found on PATH — skipping CM alignment.", file=sys.stderr)
        return [], "", None, None
    os.makedirs(output_dir, exist_ok=True)
    out_file    = os.path.join(output_dir, "patent_cmsearch.out")
    tblout_file = os.path.join(output_dir, "patent_cmsearch.tblout")
    sto_file    = os.path.join(output_dir, "patent_cmsearch.sto")
    results: List[CmsearchResult] = []
    with tempfile.TemporaryDirectory() as tmpdir:
        query_fa = os.path.join(tmpdir, "queries.fa")
        with open(query_fa, "w") as fh:
            for qname, qseq in queries:
                fh.write(f">{qname}\n{qseq}\n")
        cmd = [
            "cmsearch",
            "--notextw",
            "-A",      sto_file,
            "-o",      out_file,
            "--tblout", tblout_file,
            "-E",      "1000",
            cm_file,
            query_fa,
        ]
        print(f"  cmsearch command: {' '.join(cmd)}", file=sys.stderr)
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if r.returncode != 0:
            print(f"cmsearch error:\n{r.stderr}", file=sys.stderr)
        else:
            print(f"  cmsearch output → {out_file}", file=sys.stderr)
            print(f"  tblout          → {tblout_file}", file=sys.stderr)
            print(f"  Stockholm aln   → {sto_file}", file=sys.stderr)
        hits = _parse_tblout(tblout_file)
        for h in hits:
            results.append(CmsearchResult(
                query_name=h["query_name"],
                target_name=h["target_name"],
                cm_sig_name=h["target_name"] if float(h["evalue"]) <= 0.1 else None,
                target_accession=h["target_accession"],
                score=h["score"],
                evalue=h["evalue"],
                bias=h["bias"],
                strand=h["strand"],
                seq_from=h["seq_from"],
                seq_to=h["seq_to"],
                gc=h["gc"],
                description=h["description"],
            ))
        raw_out = ""
        if os.path.exists(out_file):
            try:
                with open(out_file) as fh:
                    raw_out = fh.read()
            except Exception:
                pass
    return results, raw_out, tblout_file, sto_file


def evaluate_cmscores(
    cm_model: Path,
    utr_ids: List[str],
    population: List[str],
    work_dir: Path,
    evalue_threshold: float,
    generation: Optional[int] = None,
) -> Tuple[List[Optional[float]], List[bool]]:
    """
    Run cmsearch for the whole population's 3' UTR sequences against
    *cm_model* and return (cmscores, hit_flags) aligned with utr_ids.

    cmsearch itself is run permissively (-E 1000, i.e. "report everything"),
    and hits are then filtered in Python against *evalue_threshold* so the
    raw tblout/Stockholm output always contains the full picture even when
    the GA's hit/no-hit cutoff is stricter. cmscores[i] is None and
    hit_flags[i] is False when sequence i had no hit at or below
    evalue_threshold (including when it had no cmsearch hit at all).
    """
    queries = list(zip(utr_ids, population))
    gen_tag = f"gen_{generation:04d}" if generation is not None else "run"
    cmsearch_out_dir = work_dir / 'cmsearch_logs' / gen_tag

    results, _raw_out, _tblout_file, _sto_file = run_cmsearch(
        queries, str(cm_model), str(cmsearch_out_dir)
    )

    # Keep the most significant (lowest E-value) hit per individual.
    best_by_target: Dict[str, CmsearchResult] = {}
    for r in results:
        try:
            ev = float(r.evalue)
        except (TypeError, ValueError):
            continue
        current = best_by_target.get(r.target_name)
        if current is None or ev < float(current.evalue):
            best_by_target[r.target_name] = r

    cmscores: List[Optional[float]] = []
    hit_flags: List[bool] = []
    for sid in utr_ids:
        r = best_by_target.get(sid)
        if r is not None and float(r.evalue) <= evalue_threshold:
            cmscores.append(float(r.score))
            hit_flags.append(True)
        else:
            cmscores.append(None)
            hit_flags.append(False)

    n_hits = sum(hit_flags)
    log.info(
        f"  cmsearch: {n_hits}/{len(utr_ids)} sequences hit the original "
        f"structure CM (E <= {evalue_threshold})"
    )
    return cmscores, hit_flags


# ══════════════════════════════════════════════════════════════════════════════
# Combined fitness: predicted half-life + structural cmscore
# ══════════════════════════════════════════════════════════════════════════════

def compute_fitness(
    halflife_scores: List[float],
    cmscores: List[Optional[float]],
    hit_flags: List[bool],
    patent_pids: List[List[float]],
    patent_pid_flags: List[bool],
    halflife_weight: float,
    cm_weight: float,
    patent_weight: float,
    no_hit_penalty: float,
    patent_pid_penalty: float,
) -> Tuple[List[float], List[float], List[float], List[float]]:
    """
    Combine predicted half-life, cmscore, and patent PID into a single fitness value.

    All three components are min-max normalised to [0, 1] across the *current*
    population before being combined, so the three weights are comparable
    regardless of the raw score scales.

    Patent PID filtering:
      - Individuals with >= 80% PID to any patent receive *patent_pid_penalty*
        (overriding the weighted combination), guaranteeing they lose selection.
      - Individuals with < 80% PID to all patents participate in normal
        weighted fitness.

    CM hit filtering (existing):
      - Individuals with no CM hit receive *no_hit_penalty*.

    Returns (fitness, normalised_halflife, normalised_cmscore, normalised_patent_pid).
    """
    valid_hl = [s for s in halflife_scores if s is not None and s != float('-inf')]
    hl_min, hl_max = (min(valid_hl), max(valid_hl)) if valid_hl else (0.0, 1.0)
    hl_range = (hl_max - hl_min) or 1.0

    valid_cm = [c for c, hit in zip(cmscores, hit_flags) if hit and c is not None]
    cm_min, cm_max = (min(valid_cm), max(valid_cm)) if valid_cm else (0.0, 1.0)
    cm_range = (cm_max - cm_min) or 1.0

    # Patent PID: lower is better (we want < 80%), so we normalize as (1 - normalized_pid)
    # Extract the max PID (worst match to any patent) for each sequence
    max_pids = [max(pids) if pids else 0.0 for pids in patent_pids]
    valid_pp = [p for p in max_pids if p is not None]
    pp_min, pp_max = (min(valid_pp), max(valid_pp)) if valid_pp else (0.0, 100.0)
    pp_range = (pp_max - pp_min) or 1.0

    fitness: List[float] = []
    norm_hl_list: List[float] = []
    norm_cm_list: List[float] = []
    norm_pp_list: List[float] = []

    for hl, cm, hit, max_pid, patent_pass in zip(
        halflife_scores, cmscores, hit_flags, max_pids, patent_pid_flags
    ):
        norm_hl = 0.0 if (hl is None or hl == float('-inf')) else (hl - hl_min) / hl_range
        norm_hl_list.append(norm_hl)

        # Patent PID penalty: if >= 80% to any patent, penalize heavily
        if not patent_pass:
            norm_pp_list.append(0.0)
            fitness.append(patent_pid_penalty)
            continue

        # For passing sequences, normalize patent PID such that lower is better
        # (0 at pp_min, 1 at pp_max); we want to reward sequences with low max_pid
        norm_pp = 1.0 - (max_pid - pp_min) / pp_range
        norm_pp_list.append(norm_pp)

        # No CM hit penalty
        if not hit:
            norm_cm_list.append(0.0)
            fitness.append(no_hit_penalty)
            continue

        norm_cm = (cm - cm_min) / cm_range
        norm_cm_list.append(norm_cm)

        # All three constraints passed: weighted combination
        combined = (
            halflife_weight * norm_hl +
            cm_weight * norm_cm +
            patent_weight * norm_pp
        )
        fitness.append(combined)

    return fitness, norm_hl_list, norm_cm_list, norm_pp_list


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
    """Single-exon transcript: no introns, no internal exons.

    Column set and names must match the real metrics/architecture.py plugin
    output exactly — transcript_id, gene_id, strand, n_exons,
    n_internal_exons, first_exon_length, last_exon_length,
    internal_exon_mean/median/sd, intron_mean/median/sd — since
    01c_predict.py's feature merge requires these specific columns and will
    abort with "feature(s) missing from merged data" if any are absent.
    For a single-exon transcript (our no-GFF fallback case) n_exons=1,
    n_internal_exons=0, and all internal-exon/intron stat columns are NA
    (left blank here so they read back as NaN), matching the real plugin's
    documented single-exon behaviour.
    """
    return {
        'transcript_id':        transcript_id,
        'gene_id':              transcript_id,
        'strand':               '+',
        'n_exons':              1,
        'n_internal_exons':     0,
        'first_exon_length':    '',
        'last_exon_length':     '',
        'internal_exon_mean':   '',
        'internal_exon_median': '',
        'internal_exon_sd':     '',
        'intron_mean':          '',
        'intron_median':        '',
        'intron_sd':            '',
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
    Plot combined fitness (top panel) and its predicted-half-life /
    cmscore components plus the annealed mutation rate (bottom panel) per
    generation, and save as PNG.  Uses matplotlib with a non-interactive
    backend so it works headlessly.  Failure to plot (e.g. matplotlib not
    installed) is logged as a warning, never a fatal error — the GA's
    TSV/FASTA outputs are unaffected.
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
        gens        = [row['generation'] for row in best_per_gen]
        best_fit    = [row['best_fitness'] for row in best_per_gen]
        mean_fit    = [row['mean_fitness'] for row in best_per_gen]
        worst_fit   = [row['worst_fitness'] for row in best_per_gen]
        best_hl     = [row['best_halflife'] for row in best_per_gen]
        best_cm     = [row['best_cmscore'] if row['best_cmscore'] is not None else float('nan')
                       for row in best_per_gen]
        mut_rate    = [row['mutation_rate'] for row in best_per_gen]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8.5), dpi=150, sharex=True)

        # ── Top panel: combined fitness ────────────────────────────────
        ax1.plot(gens, best_fit, label='Best fitness', color='#1b9e3e', linewidth=2, marker='o', markersize=3)
        ax1.plot(gens, mean_fit, label='Mean fitness', color='#2166ac', linewidth=1.5, marker='o', markersize=3)
        ax1.plot(gens, worst_fit, label='Worst fitness', color='#b2182b', linewidth=1.5, linestyle='--', alpha=0.7)
        ax1.fill_between(gens, worst_fit, best_fit, color='#1b9e3e', alpha=0.07)

        best_idx = max(range(len(best_fit)), key=lambda i: best_fit[i])
        ax1.scatter([gens[best_idx]], [best_fit[best_idx]], color='gold', edgecolor='black',
                    zorder=5, s=80, marker='*', label=f"All-time best (gen {gens[best_idx]})")

        ax1.set_ylabel('Combined fitness')
        ax1.set_title("3' UTR GA — Combined Fitness (half-life + cmscore) Over Generations")
        ax1.legend(loc='best', framealpha=0.9)
        ax1.grid(True, alpha=0.3)

        # ── Bottom panel: fitness components + annealed mutation rate ──
        ax2.plot(gens, best_hl, label='Best predicted half-life', color='#762a83', linewidth=1.5, marker='o', markersize=3)
        ax2.set_ylabel('Predicted half-life')
        ax2.tick_params(axis='y', labelcolor='#762a83')

        ax2b = ax2.twinx()
        ax2b.plot(gens, best_cm, label='Best cmscore (bits)', color='#e08214', linewidth=1.5, marker='s', markersize=3)
        ax2b.set_ylabel('cmscore (bits)')
        ax2b.tick_params(axis='y', labelcolor='#e08214')

        ax2c = ax2.twinx()
        ax2c.spines['right'].set_position(('outward', 60))
        ax2c.plot(gens, mut_rate, label='Mutation rate', color='#4d4d4d', linewidth=1.2, linestyle=':')
        ax2c.set_ylabel('Mutation rate')

        lines1, labels1 = ax2.get_legend_handles_labels()
        lines2, labels2 = ax2b.get_legend_handles_labels()
        lines3, labels3 = ax2c.get_legend_handles_labels()
        ax2.legend(lines1 + lines2 + lines3, labels1 + labels2 + labels3, loc='best', framealpha=0.9, fontsize=8)

        ax2.set_xlabel('Generation')
        ax2.grid(True, alpha=0.3)
        if len(gens) <= 30:
            ax2.set_xticks(gens)

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
    mutation_rate_start: float,
    mutation_rate_end: float,
    cooling_schedule: str,
    cooling_rate: float,
    crossover_rate: float,
    elite_frac: float,
    tournament_k: int,
    rng_seed: int,
    cm_model: Optional[Path],
    cm_evalue_threshold: float,
    halflife_weight: float,
    cm_weight: float,
    no_hit_penalty: float,
    patent_seqs: Optional[List[dict]] = None,
    patent_pid_threshold: float = 80.0,
    patent_weight: float = 0.0,
    patent_pid_penalty: float = -1000.0,
    make_plot: bool = True,
) -> None:
    rng = random.Random(rng_seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    results_dir = output_dir / 'results'
    results_dir.mkdir(exist_ok=True)

    n_elite = max(1, int(population_size * elite_frac))
    use_cm = cm_model is not None
    if patent_seqs is None:
        patent_seqs = []
    use_patent = len(patent_seqs) > 0

    log.info("═" * 60)
    log.info("3' UTR Genetic Algorithm")
    log.info(f"  Population        : {population_size}")
    log.info(f"  Generations       : {generations}")
    log.info(f"  Mutation rate     : {mutation_rate_start:.4f} → {mutation_rate_end:.4f} "
             f"({cooling_schedule}, k={cooling_rate})")
    log.info(f"  Crossover         : {crossover_rate:.4f}")
    log.info(f"  Elites            : {n_elite}")
    log.info(f"  Seed 3UTR         : {len(seed_3utr_seq)} nt")
    if use_cm:
        log.info(f"  CM model          : {cm_model}")
        log.info(f"  cmsearch E-value  : <= {cm_evalue_threshold}")
        log.info(f"  Fitness weights   : halflife={halflife_weight}, cmscore={cm_weight}, pid={patent_weight}")
        if use_patent:
            log.info(f", patent={patent_weight}")
        else:
            log.info(f", no patent")
        log.info(f"  No-hit penalty    : {no_hit_penalty}")
    else:
        log.info("  CM model          : (none) — cmscore constraint disabled")
        if use_patent:
            log.info(f"  Fitness weights   : halflife={halflife_weight}, patent={patent_weight}")
        else:
            log.info("  Fitness weights   : halflife only")
    if use_patent:
        log.info(f"  Patent sequences  : {len(patent_seqs)} motifs")
        log.info(f"  Patent PID thresh : < {patent_pid_threshold}%")
        log.info(f"  Patent penalty    : {patent_pid_penalty}")
    else:
        log.info("  Patent sequences  : (none) — patent PID constraint disabled")
    log.info(f"  Output            : {output_dir}")
    log.info("═" * 60)

    init_mutation_rate = min(
        compute_mutation_rate(1, generations, mutation_rate_start, mutation_rate_end,
                               cooling_schedule, cooling_rate) * 5,
        0.75,
    )
    population = [seed_3utr_seq] + [
        point_mutate(seed_3utr_seq, init_mutation_rate, rng)
        for _ in range(population_size - 1)
    ]

    all_rows: List[dict] = []
    best_per_gen: List[dict] = []
    all_time_best_fitness = float('-inf')
    all_time_best_seq     = seed_3utr_seq
    all_time_best_halflife = None
    all_time_best_cmscore  = None
    all_time_best_patent_pid = None

    work_dir = output_dir / '_workspace'
    work_dir.mkdir(exist_ok=True)

    for gen in range(1, generations + 1):
        gen_start = time.time()
        current_mutation_rate = compute_mutation_rate(
            gen, generations, mutation_rate_start, mutation_rate_end,
            cooling_schedule, cooling_rate,
        )
        log.info(f"\n── Generation {gen}/{generations} "
                 f"(mutation rate = {current_mutation_rate:.5f}) ──")

        # Unique sample IDs so the merge in predict_script works
        utr_ids = [f"ind_{gen:04d}_{i:05d}" for i in range(population_size)]

        halflife_scores = evaluate_population(
            population, utr_ids,
            fixed_5utr_id, fixed_5utr_seq,
            fixed_cds_id, fixed_cds_seq,
            species, metrics_script, predict_script, work_dir,
        )

        if halflife_scores is None:
            log.error(f"Pipeline failed on generation {gen}. Aborting.")
            if best_per_gen and make_plot:
                plot_fitness_over_generations(
                    best_per_gen, results_dir / 'fitness_over_generations.png'
                )
            sys.exit(1)

        # ── Structural constraint: cmsearch against the original 3' UTR CM ──
        if use_cm:
            cmscores, hit_flags = evaluate_cmscores(
                cm_model, utr_ids, population,
                work_dir, cm_evalue_threshold, generation=gen,
            )
        else:
            cmscores = [None] * population_size
            hit_flags = [True] * population_size

        # ── Patent PID constraint: avoid re-patenting existing motifs ──
        if use_patent:
            all_pids = []
            patent_pid_flags = []
            for seq in population:
                pids, passes = check_patent_pid(seq, patent_seqs, patent_pid_threshold)
                all_pids.append(pids)
                patent_pid_flags.append(passes)
        else:
            all_pids = [[] for _ in range(population_size)]
            patent_pid_flags = [True] * population_size

        fitness, norm_hl, norm_cm, norm_pp = compute_fitness(
            halflife_scores, cmscores, hit_flags,
            all_pids, patent_pid_flags,
            halflife_weight, cm_weight, patent_weight,
            no_hit_penalty, patent_pid_penalty,
        )

        # ── Record results ────────────────────────────────────────────────
        for sid, seq, hl, cm, hit, pids, p_pass, f, nhl, ncm, npp in zip(
            utr_ids, population, halflife_scores, cmscores, hit_flags,
            all_pids, patent_pid_flags, fitness, norm_hl, norm_cm, norm_pp,
        ):
            max_pid = max(pids) if pids else None
            all_rows.append({
                'generation': gen,
                'rank_in_gen': 0,
                'sample_id': sid,
                'fitness': f,
                'predicted_halflife': hl,
                'cmscore': cm if cm is not None else '',
                'cm_hit': hit,
                'max_patent_pid': max_pid if max_pid is not None else '',
                'patent_pid_pass': p_pass,
                'norm_halflife': nhl,
                'norm_cmscore': ncm,
                'norm_patent_pid': npp,
                'mutation_rate': current_mutation_rate,
                'sequence': seq,
                'length': len(seq),
                'mutations_from_seed': sum(a != b for a, b in zip(seq, seed_3utr_seq)
                                           if len(seq) == len(seed_3utr_seq)),
            })

        # Sort by combined fitness descending
        ranked = sorted(
            zip(fitness, population, utr_ids, halflife_scores, cmscores, hit_flags,
                all_pids, patent_pid_flags),
            key=lambda x: x[0], reverse=True,
        )
        for rank, (f, sq, sid, hl, cm, hit, pids, p_pass) in enumerate(ranked):
            for row in all_rows:
                if row['sample_id'] == sid:
                    row['rank_in_gen'] = rank + 1

        best_fitness, best_seq, best_id, best_halflife, best_cmscore, best_hit, best_pids, best_p_pass = ranked[0]
        best_max_pid = max(best_pids) if best_pids else None
        gen_elapsed = time.time() - gen_start
        n_hits = sum(hit_flags) if use_cm else population_size
        n_patent_pass = sum(patent_pid_flags) if use_patent else population_size
        
        patent_str = ""
        if use_patent and best_max_pid is not None:
            patent_str = f", max_patent_pid={best_max_pid:.1f}%"
        
        log.info(f"  Best fitness this gen : {best_fitness:.4f}  ({best_id}, "
                 f"halflife={best_halflife:.4f}, "
                 f"cmscore={'n/a' if best_cmscore is None else f'{best_cmscore:.2f}'}{patent_str})")
        log.info(f"  Worst fitness this gen: {ranked[-1][0]:.4f}")
        log.info(f"  Mean fitness          : {sum(fitness)/len(fitness):.4f}")
        if use_cm:
            log.info(f"  CM hits               : {n_hits}/{population_size}")
        if use_patent:
            log.info(f"  Patent PID pass       : {n_patent_pass}/{population_size}")
        log.info(f"  Elapsed               : {gen_elapsed:.1f}s")

        best_per_gen_dict = {
            'generation': gen,
            'best_fitness': best_fitness,
            'mean_fitness': sum(fitness) / len(fitness),
            'worst_fitness': ranked[-1][0],
            'best_halflife': best_halflife,
            'best_cmscore': best_cmscore,
            'best_max_patent_pid': best_max_pid,
            'n_cm_hits': n_hits,
            'n_patent_pass': n_patent_pass,
            'mutation_rate': current_mutation_rate,
            'best_id': best_id,
            'best_sequence': best_seq,
        }
        best_per_gen.append(best_per_gen_dict)

        if best_fitness > all_time_best_fitness:
            all_time_best_fitness  = best_fitness
            all_time_best_seq      = best_seq
            all_time_best_halflife = best_halflife
            all_time_best_cmscore  = best_cmscore
            all_time_best_patent_pid = best_max_pid
            log.info(f"  ★ New all-time best fitness: {all_time_best_fitness:.4f}")

        if gen == generations:
            break

        # ── Breed next generation ─────────────────────────────────────────
        elite_seqs = [sq for _, sq, _, _, _, _, _, _ in ranked[:n_elite]]
        pool_fitness = [f for f, _, _, _, _, _, _, _ in ranked]
        pool_seqs    = [sq for _, sq, _, _, _, _, _, _ in ranked]

        next_gen = list(elite_seqs)

        while len(next_gen) < population_size:
            parent_a = tournament_select(pool_seqs, pool_fitness, tournament_k, rng)
            if rng.random() < crossover_rate:
                parent_b = tournament_select(pool_seqs, pool_fitness, tournament_k, rng)
                child = uniform_crossover(parent_a, parent_b, rng)
            else:
                child = parent_a

            child = point_mutate(child, current_mutation_rate, rng)
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
        fh.write(f"Generations         : {generations}\n")
        fh.write(f"Population size     : {population_size}\n")
        fh.write(f"Mutation rate       : {mutation_rate_start} → {mutation_rate_end} "
                 f"({cooling_schedule} cooling, k={cooling_rate})\n")
        fh.write(f"Crossover rate      : {crossover_rate}\n")
        fh.write(f"Elite fraction      : {elite_frac}  ({n_elite} individuals)\n")
        fh.write(f"Tournament k        : {tournament_k}\n")
        fh.write(f"RNG seed            : {rng_seed}\n")
        if use_cm:
            fh.write(f"CM model            : {cm_model}\n")
            fh.write(f"cmsearch E-value    : <= {cm_evalue_threshold}\n")
            fh.write(f"Fitness weights     : halflife={halflife_weight}, cmscore={cm_weight}")
            if use_patent:
                fh.write(f", patent={patent_weight}\n")
            else:
                fh.write("\n")
            fh.write(f"No-hit penalty      : {no_hit_penalty}\n")
        else:
            fh.write("CM model            : (none) — cmscore constraint disabled\n")
            if use_patent:
                fh.write(f"Fitness weights     : halflife={halflife_weight}, patent={patent_weight}\n")
        if use_patent:
            fh.write(f"Patent sequences    : {len(patent_seqs)} motifs\n")
            fh.write(f"Patent PID threshold: < {patent_pid_threshold}%\n")
            fh.write(f"Patent penalty      : {patent_pid_penalty}\n")
        fh.write("\n")
        fh.write(f"Seed 3UTR length : {len(seed_3utr_seq)} nt\n")
        fh.write(f"Best 3UTR length : {len(all_time_best_seq)} nt\n\n")
        fh.write(f"Generation 1 best fitness : {best_per_gen[0]['best_fitness']:.4f}\n")
        fh.write(f"Final best fitness        : {all_time_best_fitness:.4f}\n")
        fh.write(f"Final best predicted half-life : {all_time_best_halflife:.4f}\n")
        if all_time_best_cmscore is not None:
            fh.write(f"Final best cmscore (bits)      : {all_time_best_cmscore:.2f}\n")
        if all_time_best_patent_pid is not None:
            fh.write(f"Final best max patent PID      : {all_time_best_patent_pid:.1f}%\n")
        fh.write("\n")
        mutations = sum(a != b for a, b in zip(seed_3utr_seq, all_time_best_seq)
                        if len(seed_3utr_seq) == len(all_time_best_seq))
        fh.write(f"Point mutations from seed : {mutations} / {len(seed_3utr_seq)}\n\n")
        fh.write("Generation-by-generation best fitness:\n")
        for row in best_per_gen:
            cm_str = 'n/a' if row['best_cmscore'] is None else f"{row['best_cmscore']:.2f}"
            patent_str = f" | patent_pass {row['n_patent_pass']}" if use_patent else ""
            fh.write(
                f"  Gen {row['generation']:>3d} : fitness {row['best_fitness']:.4f} "
                f"(mean {row['mean_fitness']:.4f}, worst {row['worst_fitness']:.4f}) | "
                f"halflife {row['best_halflife']:.4f} | cmscore {cm_str} | "
                f"mut_rate {row['mutation_rate']:.5f} | cm_hits {row['n_cm_hits']}{patent_str}\n"
            )
        fh.write("\nBest evolved 3UTR sequence:\n")
        for i in range(0, len(all_time_best_seq), 60):
            fh.write(all_time_best_seq[i:i + 60] + '\n')

    log.info(f"Wrote {summary_path}")
    log.info("\n" + "═" * 60)
    patent_str_final = ""
    if all_time_best_patent_pid is not None:
        patent_str_final = f", max_patent_pid={all_time_best_patent_pid:.1f}%"
    log.info(f"GA complete.  All-time best fitness: {all_time_best_fitness:.4f} "
             f"(halflife={all_time_best_halflife:.4f}, "
             f"cmscore={'n/a' if all_time_best_cmscore is None else f'{all_time_best_cmscore:.2f}'}{patent_str_final})")
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
    parser.add_argument('--mutation-rate-start', type=float, default=0.05,
                        dest='mutation_rate_start', metavar='RATE',
                        help="Per-nucleotide mutation probability at generation 1 "
                             "(default: 0.05).")
    parser.add_argument('--mutation-rate-end', type=float, default=0.005,
                        dest='mutation_rate_end', metavar='RATE',
                        help="Per-nucleotide mutation probability at the final "
                             "generation (default: 0.005).")
    parser.add_argument('--cooling-schedule', choices=['exponential', 'linear'],
                        default='exponential', dest='cooling_schedule',
                        help="How the mutation rate is annealed from start to end "
                             "(default: exponential, simulated-annealing style).")
    parser.add_argument('--cooling-rate', type=float, default=3.0,
                        dest='cooling_rate', metavar='K',
                        help="Decay constant for exponential cooling; larger values "
                             "cool faster early on (default: 3.0). Ignored for "
                             "--cooling-schedule linear.")
    parser.add_argument('--mutation-rate', type=float, default=None,
                        dest='mutation_rate_legacy', metavar='RATE',
                        help="Deprecated: sets a single constant mutation rate "
                             "(disables annealing) for backward compatibility.")
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

    # ── Structural (cmscore) constraint ──────────────────────────────────
    parser.add_argument('--cm-model', default=None, dest='cm_model', metavar='FILE',
                        help="Calibrated Infernal covariance model (.cm, from "
                             "cmbuild + cmcalibrate on the original 3' UTR structure). "
                             "If omitted, the structural constraint is skipped and "
                             "fitness is normalised predicted half-life alone.")
    parser.add_argument('--cm-evalue-threshold', type=float, default=0.01,
                        dest='cm_evalue_threshold', metavar='EVALUE',
                        help="E-value cutoff for a cmsearch hit against the original "
                             "structure (default: 0.01). cmsearch itself is run "
                             "permissively (-E 1000) so the full tblout/Stockholm "
                             "output is preserved; this threshold is applied "
                             "afterward to decide hit vs. no-hit. Sequences with no "
                             "hit at or below this threshold receive --no-hit-penalty "
                             "as fitness.")
    parser.add_argument('--halflife-weight', type=float, default=0.5,
                        dest='halflife_weight', metavar='W',
                        help="Weight of normalised predicted half-life in the combined "
                             "fitness function (default: 0.5).")
    parser.add_argument('--cm-weight', type=float, default=0.5,
                        dest='cm_weight', metavar='W',
                        help="Weight of normalised cmscore in the combined fitness "
                             "function (default: 0.5).")
    parser.add_argument('--no-hit-penalty', type=float, default=-1000.0,
                        dest='no_hit_penalty', metavar='VALUE',
                        help="Fitness assigned to individuals with no cmsearch hit "
                             "against the original structure (default: -1000.0). "
                             "Should be well below any achievable weighted combination "
                             "of normalised half-life and cmscore (which lies in "
                             "[0, halflife-weight + cm-weight]).")

    # ── Patent sequence constraint ───────────────────────────────────────
    parser.add_argument('--use-default-patents', action='store_true',
                        dest='use_default_patents',
                        help="Use built-in A7 patent sequences (default: False). "
                             "Sequences with >= 80%% PID to any patent receive a penalty.")
    parser.add_argument('--patent-pid-threshold', type=float, default=80.0,
                        dest='patent_pid_threshold', metavar='PCT',
                        help="PID threshold for patent filtering (default: 80.0). "
                             "Sequences with >= this PID to any patent are penalized.")
    parser.add_argument('--patent-weight', type=float, default=0.0,
                        dest='patent_weight', metavar='W',
                        help="Weight of normalised patent PID avoidance in combined "
                             "fitness (default: 0.0, i.e., patent constraint is soft "
                             "penalty-only). Set > 0 to actively reward low PID.")
    parser.add_argument('--patent-pid-penalty', type=float, default=-1000.0,
                        dest='patent_pid_penalty', metavar='VALUE',
                        help="Fitness penalty assigned to sequences with >= patent-pid-threshold "
                             "PID to any patent (default: -1000.0).")

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

    # ── Backward-compat: a single --mutation-rate disables annealing ──────
    mutation_rate_start = args.mutation_rate_start
    mutation_rate_end = args.mutation_rate_end
    if args.mutation_rate_legacy is not None:
        log.warning(
            "--mutation-rate is deprecated; using it as a constant rate "
            "(annealing disabled). Prefer --mutation-rate-start/--mutation-rate-end."
        )
        mutation_rate_start = args.mutation_rate_legacy
        mutation_rate_end = args.mutation_rate_legacy

    # ── Validate paths ────────────────────────────────────────────────────
    required_paths = [
        ('--fasta-5utr', args.fasta_5utr),
        ('--fasta-cds',  args.fasta_cds),
        ('--fasta-3utr', args.fasta_3utr),
        ('--metrics-script', args.metrics_script),
        ('--predict-script', args.predict_script),
    ]
    for label, path_str in required_paths:
        p = Path(path_str)
        if not p.exists():
            log.error(f"{label} not found: {p}")
            sys.exit(1)

    cm_model_path: Optional[Path] = None
    if args.cm_model is not None:
        cm_model_path = Path(args.cm_model)
        if not cm_model_path.exists():
            log.error(f"--cm-model not found: {cm_model_path}")
            sys.exit(1)
        if shutil.which("cmsearch") is None:
            log.warning(
                "--cm-model was given but 'cmsearch' was not found on PATH — "
                "the structural constraint will be skipped (fitness = "
                "normalised half-life only)."
            )
            cm_model_path = None

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
        mutation_rate_start = mutation_rate_start,
        mutation_rate_end   = mutation_rate_end,
        cooling_schedule    = args.cooling_schedule,
        cooling_rate        = args.cooling_rate,
        crossover_rate  = args.crossover_rate,
        elite_frac      = args.elite_frac,
        tournament_k    = args.tournament_k,
        rng_seed        = args.seed,
        cm_model            = cm_model_path,
        cm_evalue_threshold = args.cm_evalue_threshold,
        halflife_weight     = args.halflife_weight,
        cm_weight           = args.cm_weight,
        no_hit_penalty      = args.no_hit_penalty,
        patent_seqs         = PATENT_SEQUENCES if args.use_default_patents else None,
        patent_pid_threshold = args.patent_pid_threshold,
        patent_weight       = args.patent_weight,
        patent_pid_penalty  = args.patent_pid_penalty,
        make_plot       = not args.no_plot,
    )


if __name__ == '__main__':
    main()