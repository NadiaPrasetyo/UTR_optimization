#!/usr/bin/env python3
"""
Upstream HMM Pipeline
=====================
1. Parse hit coordinates from a Stockholm alignment file
2. Extract upstream (5') sequence from each hit using reference genomes
3. Align upstream sequences with MAFFT
4. Build an HMM profile with hmmbuild
5. Search the HMM profile against a sequence database with hmmsearch

Dependencies: MAFFT, HMMER3
  conda install -c bioconda mafft hmmer
  OR
  apt install hmmer mafft

Usage:
  python upstream_hmm_pipeline.py \\
    --ref-genomes   reference_genomes.fasta \\
    --query-stk     query.stk \\
    --db            sequence_database.fasta \\
    --outdir        results/ \\
    [--upstream-len 500] \\
    [--threads 4]
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
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
            "Install with:  conda install -c bioconda mafft hmmer"
        )
    return path


def run(cmd: list[str], desc: str = "") -> subprocess.CompletedProcess:
    """Run a command, log stderr at DEBUG level, and abort on failure."""
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
    """Parse a FASTA file → {header: sequence} (header = full description line)."""
    seqs: dict[str, str] = {}
    current_header: str | None = None
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
    """Write {header: sequence} to a FASTA file (80-char line wrap)."""
    with open(path, "w") as fh:
        for header, seq in seqs.items():
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def fasta_seq_id(header: str) -> str:
    """Return the first whitespace-delimited token of a FASTA header."""
    return header.split()[0]


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]

# ── Step 1 – Parse Stockholm coordinates ─────────────────────────────────────

def parse_stockholm_hits(stk_file: str) -> list[dict]:
    """
    Extract unique hit coordinates from a Stockholm alignment file.

    Sequence rows are named  ACCESSION/START-END  (1-based, always start < end
    because Rfam/Infernal always reports on the forward strand of the stored
    sequence).  Strand is inferred: if the genome is stored 5'→3' then
    start < end is plus strand; the Stockholm format itself does not encode
    strand, so we assume plus throughout unless you have external evidence.

    Returns a list of dicts:
        accession : str   – sequence accession (e.g. NC_024120.1)
        start     : int   – 1-based start (always ≤ end)
        end       : int   – 1-based end
        strand    : str   – 'plus' (default; extend if needed)
    """
    hits: list[dict] = []
    seen: set[str] = set()
    coord_re = re.compile(r"^([A-Za-z0-9_.]+)/(\d+)-(\d+)")

    with open(stk_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("//"):
                continue

            name = line.split()[0]
            if name in seen:
                continue
            seen.add(name)

            m = coord_re.match(name)
            if not m:
                continue

            accession, start, end = m.groups()
            hits.append(
                {
                    "accession": accession,
                    "start":     int(start),
                    "end":       int(end),
                    "strand":    "plus",   # Stockholm doesn't encode strand
                }
            )

    log.info("Parsed %d unique hits from %s", len(hits), stk_file)
    return hits

# ── Step 2 – Extract upstream sequences ──────────────────────────────────────

def extract_upstream(
    hits: list[dict],
    ref_genomes_fasta: str,
    upstream_len: int,
    outdir: Path,
) -> str:
    """
    For each hit, extract up to *upstream_len* nt immediately 5' of the hit.

    Plus strand  → bases [start - upstream_len, start - 1]  (1-based)
    Minus strand → bases [end + 1, end + upstream_len]  then rev-comp

    Returns the path to the extracted FASTA file.
    """
    ref_seqs = read_fasta(ref_genomes_fasta)
    # Build lookup by sequence ID (first token of header)
    ref_by_id: dict[str, str] = {
        fasta_seq_id(h): seq for h, seq in ref_seqs.items()
    }

    upstream_seqs: dict[str, str] = {}
    seen: set[str] = set()

    for hit in hits:
        acc = hit["accession"]
        if acc not in ref_by_id:
            log.warning("Accession '%s' not found in reference FASTA – skipping", acc)
            continue

        genome = ref_by_id[acc]
        genome_len = len(genome)

        if hit["strand"] == "plus":
            # hit starts at 'start' (1-based) → upstream is to the left
            hit_start_0 = hit["start"] - 1           # 0-based inclusive start of hit
            up_end_0    = hit_start_0                 # exclusive end of upstream window
            up_start_0  = max(0, hit_start_0 - upstream_len)
            subseq = genome[up_start_0:up_end_0]
        else:
            # minus strand: 'end' is the rightmost coordinate (1-based)
            hit_end_0  = hit["end"]                   # 0-based exclusive end of hit
            up_start_0 = hit_end_0                    # 0-based start of upstream window
            up_end_0   = min(genome_len, hit_end_0 + upstream_len)
            subseq = reverse_complement(genome[up_start_0:up_end_0])

        if not subseq:
            log.warning("Empty upstream region for '%s' – skipping", acc)
            continue

        # Deduplicate by accession + hit position
        key = f"{acc}_{hit['start']}_{hit['strand']}"
        if key in seen:
            continue
        seen.add(key)

        label = f"{acc}/{hit['start']}-{hit['end']}_upstream_{hit['strand']}"
        upstream_seqs[label] = subseq
        log.info(
            "  Extracted %d nt upstream of %s (%s strand)",
            len(subseq), acc, hit["strand"],
        )

    if not upstream_seqs:
        sys.exit("ERROR: no upstream sequences could be extracted. "
                 "Check that accessions in the Stockholm file match the reference FASTA.")

    out_path = str(outdir / "upstream_sequences.fasta")
    write_fasta(upstream_seqs, out_path)
    log.info("Upstream sequences written to %s (%d total)", out_path, len(upstream_seqs))
    return out_path

# ── Step 3 – MSA with MAFFT ──────────────────────────────────────────────────

def build_msa(
    upstream_fasta: str,
    outdir: Path,
    mafft_args: list[str] | None = None,
) -> str:
    """
    Align upstream sequences with MAFFT.
    Returns the path to the aligned FASTA (.afa).
    """
    msa_path = str(outdir / "upstream_msa.afa")
    mafft_bin = require_tool("mafft")
    extra = mafft_args if mafft_args else ["--auto"]

    cmd = [mafft_bin, "--thread", "-1"] + extra + ["--out", msa_path, upstream_fasta]
    run(cmd, desc="MAFFT (MSA of upstream sequences)")
    log.info("MSA written to %s", msa_path)
    return msa_path

# ── Step 4 – hmmbuild ────────────────────────────────────────────────────────

def build_hmm_profile(
    msa_path: str,
    outdir: Path,
    profile_name: str = "upstream_profile",
) -> str:
    """
    Build an HMM profile from the upstream MSA.
    Returns the path to the .hmm file.
    """
    hmm_path = str(outdir / f"{profile_name}.hmm")
    run(
        ["hmmbuild", "--dna", "-n", profile_name, hmm_path, msa_path],
        desc="hmmbuild (HMM profile from upstream MSA)",
    )
    log.info("HMM profile written to %s", hmm_path)
    return hmm_path

# ── Step 5 – hmmsearch ───────────────────────────────────────────────────────

def run_hmmsearch(
    hmm_path: str,
    db_fasta: str,
    outdir: Path,
    threads: int,
    profile_name: str = "upstream_profile",
) -> dict[str, str]:
    """
    Search the HMM profile against the sequence database with hmmsearch.
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
            hmm_path,
            db_fasta,
        ],
        desc="hmmsearch (HMM profile vs sequence database)",
    )

    log.info("hmmsearch outputs:")
    for key, path in outputs.items():
        log.info("  %-8s → %s", key, path)

    return outputs

# ── Summary ───────────────────────────────────────────────────────────────────

def print_summary(
    hits: list[dict],
    upstream_fasta: str,
    msa_path: str,
    hmm_path: str,
    hmmsearch_outputs: dict[str, str],
) -> None:
    upstream_count = len(read_fasta(upstream_fasta))

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
    print(f"  Stockholm hits parsed        : {len(hits)}")
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
        description=(
            "Upstream HMM pipeline: Stockholm coords → upstream extraction "
            "→ MAFFT MSA → hmmbuild → hmmsearch"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--ref-genomes", required=True, metavar="FASTA",
        help="FASTA of reference genomes (accessions must match Stockholm names)",
    )
    parser.add_argument(
        "--query-stk", required=True, metavar="STK",
        help="Stockholm alignment file containing hit coordinates",
    )
    parser.add_argument(
        "--db", required=True, metavar="FASTA",
        help="Sequence database FASTA for hmmsearch",
    )
    parser.add_argument(
        "--outdir", default="pipeline_output", metavar="DIR",
        help="Output directory (created if absent)",
    )
    parser.add_argument(
        "--upstream-len", type=int, default=500, metavar="INT",
        help="Maximum upstream region length (nt) to extract",
    )
    parser.add_argument(
        "--profile-name", default="upstream_profile", metavar="STR",
        help="Base name for the HMM profile and output files",
    )
    parser.add_argument(
        "--threads", type=int, default=4, metavar="INT",
        help="CPU threads for hmmsearch",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Validate input files
    for label, path in [
        ("--ref-genomes", args.ref_genomes),
        ("--query-stk",   args.query_stk),
        ("--db",          args.db),
    ]:
        if not os.path.isfile(path):
            sys.exit(f"ERROR: {label} file not found: {path}")

    # Check required external tools
    for tool in ["mafft", "hmmbuild", "hmmsearch"]:
        require_tool(tool)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log.info("Output directory: %s", outdir.resolve())

    # ── Step 1: Parse Stockholm ───────────────────────────────────────────────
    log.info("=== Step 1: Parse hit coordinates from Stockholm file ===")
    hits = parse_stockholm_hits(args.query_stk)
    if not hits:
        sys.exit("ERROR: no hits parsed from Stockholm file.")

    # ── Step 2: Extract upstream sequences ───────────────────────────────────
    log.info("=== Step 2: Extract upstream sequences (max len=%d nt) ===", args.upstream_len)
    upstream_fasta = extract_upstream(
        hits              = hits,
        ref_genomes_fasta = args.ref_genomes,
        upstream_len      = args.upstream_len,
        outdir            = outdir,
    )

    # ── Step 3: MSA with MAFFT ───────────────────────────────────────────────
    log.info("=== Step 3: Align upstream sequences with MAFFT ===")
    msa_path = build_msa(
        upstream_fasta = upstream_fasta,
        outdir         = outdir,
    )

    # ── Step 4: Build HMM profile ────────────────────────────────────────────
    log.info("=== Step 4: Build HMM profile with hmmbuild ===")
    hmm_path = build_hmm_profile(
        msa_path     = msa_path,
        outdir       = outdir,
        profile_name = args.profile_name,
    )

    # ── Step 5: Search with hmmsearch ────────────────────────────────────────
    log.info("=== Step 5: Search HMM profile against sequence database ===")
    hmmsearch_outputs = run_hmmsearch(
        hmm_path     = hmm_path,
        db_fasta     = args.db,
        outdir       = outdir,
        threads      = args.threads,
        profile_name = args.profile_name,
    )

    # ── Summary ──────────────────────────────────────────────────────────────
    print_summary(
        hits              = hits,
        upstream_fasta    = upstream_fasta,
        msa_path          = msa_path,
        hmm_path          = hmm_path,
        hmmsearch_outputs = hmmsearch_outputs,
    )


if __name__ == "__main__":
    main()