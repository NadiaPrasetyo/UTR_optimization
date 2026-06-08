#!/usr/bin/env python3
"""
Upstream HMM Pipeline
=====================
1. Aligns a query sequence to a set of reference genomes (BLAST)
2. Extracts upstream (coding) sequence from each alignment hit
3. Builds a MSA and HMM profile from upstream sequences (MUSCLE + hmmbuild)
4. Searches the HMM profile against a sequence database (hmmsearch)

Dependencies: BLAST+, MAFFT, HMMER3
  conda install -c bioconda blast mafft hmmer
  OR
  apt install ncbi-blast+ hmmer mafft

Usage:
  python upstream_hmm_pipeline.py \\
    --ref-genomes   reference_genomes.fasta \\
    --query         query.fasta \\
    --db            sequence_database.fasta \\
    --outdir        results/ \\
    [--upstream-len 500] \\
    [--min-identity 70] \\
    [--min-coverage 50] \\
    [--threads 4]
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# ── Logging ───────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ── Helpers ───────────────────────────────────────────────────────────────────

def require_tool(name: str) -> str:
    """Return the path to *name* or abort with a helpful message."""
    path = shutil.which(name)
    if path is None:
        sys.exit(
            f"ERROR: '{name}' not found on PATH.\n"
            "Install it with:  conda install -c bioconda blast muscle hmmer"
        )
    return path


def run(cmd: list[str], desc: str = "") -> subprocess.CompletedProcess:
    """Run a shell command, stream stderr, and raise on failure."""
    log.info("%s: %s", desc or "Running", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout.strip():
        log.debug(result.stdout.strip())
    if result.stderr.strip():
        log.debug(result.stderr.strip())
    if result.returncode != 0:
        log.error("STDOUT:\n%s", result.stdout)
        log.error("STDERR:\n%s", result.stderr)
        sys.exit(f"ERROR: command failed (exit {result.returncode}): {' '.join(cmd)}")
    return result


def read_fasta(path: str) -> dict[str, str]:
    """Parse a FASTA file into {header: sequence} (header = full description line)."""
    seqs: dict[str, str] = {}
    current_header = None
    current_seq: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header is not None:
                    seqs[current_header] = "".join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
    if current_header is not None:
        seqs[current_header] = "".join(current_seq)
    return seqs


def write_fasta(seqs: dict[str, str], path: str) -> None:
    """Write {header: sequence} to a FASTA file."""
    with open(path, "w") as fh:
        for header, seq in seqs.items():
            fh.write(f">{header}\n")
            # 80-char line wrap
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def seq_id(header: str) -> str:
    """Extract the first whitespace-delimited token from a FASTA header."""
    return header.split()[0]

# ── Step 1 – BLAST query vs reference genomes ────────────────────────────────

def blast_query(
    query_fasta: str,
    ref_genomes_fasta: str,
    outdir: Path,
    threads: int,
    min_identity: float,
    min_coverage: float,
) -> list[dict]:
    """
    BLASTn query against reference genomes.
    Returns a list of hit dicts with keys:
      sseqid, sstart, send, sstrand, slen, qcovhsp, pident
    """
    db_dir = outdir / "blast_db"
    db_dir.mkdir(exist_ok=True)
    db_path = str(db_dir / "ref_db")

    # Build BLAST database from reference genomes
    run(
        ["makeblastdb", "-in", ref_genomes_fasta, "-dbtype", "nucl", "-out", db_path],
        desc="makeblastdb (reference genomes)",
    )

    blast_out = str(outdir / "blast_hits.tsv")
    fmt = "6 qseqid sseqid pident length qlen slen qstart qend sstart send sstrand qcovhsp"

    run(
        [
            "blastn",
            "-query", query_fasta,
            "-db", db_path,
            "-out", blast_out,
            "-outfmt", fmt,
            "-num_threads", str(threads),
            "-perc_identity", str(min_identity),
        ],
        desc="blastn (query vs references)",
    )

    hits = []
    with open(blast_out) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            (qseqid, sseqid, pident, length, qlen, slen,
             qstart, qend, sstart, send, sstrand, qcovhsp) = parts[:12]

            if float(qcovhsp) < min_coverage:
                continue

            hits.append(
                {
                    "qseqid":  qseqid,
                    "sseqid":  sseqid,
                    "pident":  float(pident),
                    "sstart":  int(sstart),
                    "send":    int(send),
                    "slen":    int(slen),
                    "sstrand": sstrand,   # "plus" or "minus"
                    "qcovhsp": float(qcovhsp),
                }
            )

    log.info("BLAST: %d hit(s) pass identity ≥ %.0f%% and coverage ≥ %.0f%%",
             len(hits), min_identity, min_coverage)
    if not hits:
        sys.exit(
            "No BLAST hits found. Try lowering --min-identity or --min-coverage."
        )
    return hits

# ── Step 2 – Extract upstream sequences ──────────────────────────────────────

def extract_upstream(
    hits: list[dict],
    ref_genomes_fasta: str,
    upstream_len: int,
    outdir: Path,
) -> str:
    """
    For each BLAST hit, extract up to *upstream_len* nt of sequence that lies
    immediately upstream (5′) of the hit on the subject genome.

    'Upstream' means:
      • plus  strand → bases [sstart - upstream_len, sstart - 1]
      • minus strand → bases [send + 1,              send + upstream_len]
                        (then reverse-complemented)

    Returns the path to the extracted FASTA file.
    """
    ref_seqs = read_fasta(ref_genomes_fasta)
    # Build a lookup by sequence ID (first token of header)
    ref_by_id: dict[str, tuple[str, str]] = {}
    for header, seq in ref_seqs.items():
        ref_by_id[seq_id(header)] = (header, seq)

    upstream_seqs: dict[str, str] = {}
    seen: set[str] = set()

    for hit in hits:
        sid = hit["sseqid"]
        if sid not in ref_by_id:
            log.warning("Subject '%s' not found in reference FASTA – skipping", sid)
            continue

        _, genome_seq = ref_by_id[sid]
        genome_len = len(genome_seq)

        if hit["sstrand"] == "plus":
            # hit starts at sstart (1-based) → upstream is to the left
            hit_start_0 = hit["sstart"] - 1          # convert to 0-based
            up_end_0    = hit_start_0                  # exclusive
            up_start_0  = max(0, hit_start_0 - upstream_len)
            subseq = genome_seq[up_start_0:up_end_0]
        else:
            # minus strand: BLAST reports sstart > send; the 'start' of the
            # gene on the minus strand is at the larger coordinate
            hit_end_0   = hit["sstart"]               # already 1-based max coord
            up_start_0  = hit_end_0                    # 0-based start (exclusive of hit)
            up_end_0    = min(genome_len, hit_end_0 + upstream_len)
            subseq = reverse_complement(genome_seq[up_start_0:up_end_0])

        if not subseq:
            log.warning("Empty upstream region for hit on '%s' – skipping", sid)
            continue

        # Deduplicate by genome + approximate position
        key = f"{sid}_{hit['sstart']}_{hit['sstrand']}"
        if key in seen:
            continue
        seen.add(key)

        label = (
            f"{sid}_upstream_{hit['sstart']}_{hit['sstrand']}"
            f"_pident{hit['pident']:.0f}"
        )
        upstream_seqs[label] = subseq
        log.info(
            "  Extracted %d nt upstream of %s (%s strand)",
            len(subseq), sid, hit["sstrand"],
        )

    if not upstream_seqs:
        sys.exit("No upstream sequences could be extracted.")

    out_path = str(outdir / "upstream_sequences.fasta")
    write_fasta(upstream_seqs, out_path)
    log.info("Upstream sequences written to %s (%d total)", out_path, len(upstream_seqs))
    return out_path


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]

# ── Step 3 – MSA with MAFFT then hmmbuild ─────────────────────────────────────

def build_hmm_profile(
    upstream_fasta: str,
    outdir: Path,
    profile_name: str = "upstream_profile",
    mafft_args: list[str] | None = None,
) -> tuple[str, str]:
    """
    Align upstream sequences with MAFFT, then build an HMM profile.
    Returns (msa_path, hmm_path).

    mafft_args: extra MAFFT flags (e.g. ["--localpair", "--maxiterate", "1000"])
                defaults to --auto if not provided.
    """
    msa_path = str(outdir / "upstream_msa.afa")
    hmm_path = str(outdir / f"{profile_name}.hmm")

    mafft_bin = require_tool("mafft")
    extra = mafft_args if mafft_args else ["--auto"]

    # MAFFT writes the alignment to stdout
    cmd = [mafft_bin, "--thread", "-1"] + extra + ["--out", msa_path, upstream_fasta]
    # --thread -1 → auto-detect CPU count
    run(cmd, desc="MAFFT (MSA of upstream sequences)")

    run(
        ["hmmbuild", "--dna", "-n", profile_name, hmm_path, msa_path],
        desc="hmmbuild (HMM profile)",
    )
    log.info("HMM profile written to %s", hmm_path)
    return msa_path, hmm_path

# ── Step 4 – hmmsearch against database ──────────────────────────────────────

def run_hmmsearch(
    hmm_path: str,
    db_fasta: str,
    outdir: Path,
    threads: int,
    profile_name: str = "upstream_profile",
) -> dict[str, str]:
    """
    Search the HMM profile against the sequence database.
    Returns a dict of output file paths.
    """
    outputs = {
        "txt":    str(outdir / f"{profile_name}_hmmsearch.txt"),
        "tblout": str(outdir / f"{profile_name}_hmmsearch_tblout.txt"),
        "msa":    str(outdir / f"{profile_name}_hmmsearch_hits.afa"),
    }

    run(
        [
            "hmmsearch",
            "--cpu",    str(threads),
            "-o",       outputs["txt"],
            "--tblout", outputs["tblout"],
            "-A",       outputs["msa"],
            "--dna",
            hmm_path,
            db_fasta,
        ],
        desc="hmmsearch (HMM vs sequence database)",
    )

    log.info("hmmsearch outputs:")
    for key, path in outputs.items():
        log.info("  %-8s → %s", key, path)

    return outputs

# ── Summary ───────────────────────────────────────────────────────────────────

def print_summary(
    blast_hits: list[dict],
    upstream_fasta: str,
    msa_path: str,
    hmm_path: str,
    hmmsearch_outputs: dict[str, str],
    outdir: Path,
) -> None:
    upstream_count = len(read_fasta(upstream_fasta))

    # Count hmmsearch hits from tblout (non-comment lines)
    tblout_hits = 0
    try:
        with open(hmmsearch_outputs["tblout"]) as fh:
            tblout_hits = sum(
                1 for line in fh
                if line.strip() and not line.startswith("#")
            )
    except FileNotFoundError:
        pass

    print("\n" + "=" * 60)
    print("  Pipeline complete – summary")
    print("=" * 60)
    print(f"  BLAST hits (passing filters) : {len(blast_hits)}")
    print(f"  Upstream sequences extracted : {upstream_count}")
    print(f"  MSA file                     : {msa_path}")
    print(f"  HMM profile                  : {hmm_path}")
    print(f"  hmmsearch hits (tblout)      : {tblout_hits}")
    print()
    print("  Output files")
    print(f"    Upstream FASTA   : {upstream_fasta}")
    print(f"    MSA (upstream)   : {msa_path}")
    print(f"    HMM profile      : {hmm_path}")
    print(f"    hmmsearch TXT    : {hmmsearch_outputs['txt']}")
    print(f"    hmmsearch TBLOUT : {hmmsearch_outputs['tblout']}")
    print(f"    hmmsearch MSA    : {hmmsearch_outputs['msa']}")
    print("=" * 60 + "\n")

# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Upstream HMM pipeline: BLAST → upstream extraction → MSA → hmmbuild → hmmsearch",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--ref-genomes",  required=True, metavar="FASTA",
                        help="FASTA of reference genomes (small set)")
    parser.add_argument("--query",        required=True, metavar="FASTA",
                        help="Query FASTA (single sequence)")
    parser.add_argument("--db",           required=True, metavar="FASTA",
                        help="Sequence database FASTA for hmmsearch")
    parser.add_argument("--outdir",       default="pipeline_output", metavar="DIR",
                        help="Output directory")
    parser.add_argument("--upstream-len", type=int, default=500, metavar="INT",
                        help="Maximum length (nt) of upstream region to extract")
    parser.add_argument("--min-identity", type=float, default=70.0, metavar="FLOAT",
                        help="Minimum BLAST percent identity to keep a hit")
    parser.add_argument("--min-coverage", type=float, default=50.0, metavar="FLOAT",
                        help="Minimum BLAST query coverage (%%) to keep a hit")
    parser.add_argument("--profile-name", default="upstream_profile", metavar="STR",
                        help="Base name for the HMM profile and output files")
    parser.add_argument("--threads",      type=int, default=4, metavar="INT",
                        help="CPU threads for BLAST and hmmsearch")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Validate inputs
    for label, path in [
        ("--ref-genomes", args.ref_genomes),
        ("--query",       args.query),
        ("--db",          args.db),
    ]:
        if not os.path.isfile(path):
            sys.exit(f"ERROR: {label} file not found: {path}")

    # Check required tools
    for tool in ["makeblastdb", "blastn", "mafft", "hmmbuild", "hmmsearch"]:
        require_tool(tool)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log.info("Output directory: %s", outdir.resolve())

    # ── Step 1: BLAST ────────────────────────────────────────────────────────
    log.info("=== Step 1: BLAST query vs reference genomes ===")
    blast_hits = blast_query(
        query_fasta      = args.query,
        ref_genomes_fasta= args.ref_genomes,
        outdir           = outdir,
        threads          = args.threads,
        min_identity     = args.min_identity,
        min_coverage     = args.min_coverage,
    )

    # ── Step 2: Extract upstream sequences ───────────────────────────────────
    log.info("=== Step 2: Extract upstream sequences (len=%d) ===", args.upstream_len)
    upstream_fasta = extract_upstream(
        hits              = blast_hits,
        ref_genomes_fasta = args.ref_genomes,
        upstream_len      = args.upstream_len,
        outdir            = outdir,
    )

    # ── Step 3: MSA + hmmbuild ───────────────────────────────────────────────
    log.info("=== Step 3: MSA (MAFFT) + hmmbuild ===")
    msa_path, hmm_path = build_hmm_profile(
        upstream_fasta = upstream_fasta,
        outdir         = outdir,
        profile_name   = args.profile_name,
    )

    # ── Step 4: hmmsearch ────────────────────────────────────────────────────
    log.info("=== Step 4: hmmsearch vs sequence database ===")
    hmmsearch_outputs = run_hmmsearch(
        hmm_path     = hmm_path,
        db_fasta     = args.db,
        outdir       = outdir,
        threads      = args.threads,
        profile_name = args.profile_name,
    )

    # ── Summary ──────────────────────────────────────────────────────────────
    print_summary(
        blast_hits        = blast_hits,
        upstream_fasta    = upstream_fasta,
        msa_path          = msa_path,
        hmm_path          = hmm_path,
        hmmsearch_outputs = hmmsearch_outputs,
        outdir            = outdir,
    )


if __name__ == "__main__":
    main()