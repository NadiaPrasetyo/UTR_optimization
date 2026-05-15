#!/usr/bin/env python3
"""
Fetch MANE transcript sequences (5'UTR, CDS, 3'UTR) from Ensembl TARK API.
Uses Ensembl Gene ID as a sanity check.
"""

import requests
import sys
import os
import json
import pandas as pd
import argparse
import time
import logging
from pathlib import Path

ENDPOINT_URL = "http://tark.ensembl.org/api/transcript/"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger(__name__)


def fetch_transcript(ensembl_nuc: str, ensembl_gene: str) -> dict | None:
    """
    Query the TARK API for a transcript by stable_id, expanding sequence and genes.
    Returns the best-matching result after validating the gene ID, or None on failure.

    Parameters
    ----------
    ensembl_nuc  : Ensembl transcript stable ID (with or without version suffix).
    ensembl_gene : Expected Ensembl gene stable ID (used as sanity check).
    """
    # Strip version suffix for the query (API accepts bare stable_id)
    stable_id = ensembl_nuc.split(".")[0]
    version = ensembl_nuc.split(".")[-1]

    params = {
        "stable_id": stable_id,
        "stable_id_version": version,
        "expand": "sequence,genes",
        "assembly_name": "GRCh38",
    }

    try:
        response = requests.get(ENDPOINT_URL, params=params, timeout=30)
        response.raise_for_status()
    except requests.RequestException as exc:
        log.error("HTTP error fetching %s: %s", stable_id, exc)
        return None

    data = response.json()
    results = data.get("results", [])

    if not results:
        log.warning("No results returned for %s", stable_id)
        return None

    # Take the first result returned
    record = results[0]

    # --- Sanity check: warn if the gene ID doesn't match, but proceed anyway ---
    gene_bare = ensembl_gene.split(".")[0]
    found_gene_ids = [g.get("stable_id", "") for g in record.get("genes", [])]
    if not any(gid.split(".")[0] == gene_bare for gid in found_gene_ids):
        log.warning(
            "Gene ID mismatch for transcript %s: expected %s, found %s — using first result anyway",
            stable_id,
            ensembl_gene,
            found_gene_ids,
        )

    return record


def extract_cds(record: dict) -> str | None:
    """
    Derive the CDS sequence by trimming the 5'UTR and 3'UTR from the full transcript sequence.
    Falls back to None if any required field is missing.
    """
    seq_block = record.get("sequence", {})
    full_seq = seq_block.get("sequence") if seq_block else None

    five_utr = record.get("five_prime_utr_seq", "")
    three_utr = record.get("three_prime_utr_seq", "")

    if not full_seq:
        return None

    five_len = len(five_utr) if five_utr else 0
    three_len = len(three_utr) if three_utr else 0

    cds_end = len(full_seq) - three_len if three_len else len(full_seq)
    cds = full_seq[five_len:cds_end]

    return cds if cds else None


def write_fasta(filepath: str, header: str, sequence: str) -> None:
    """Write a single-record FASTA file with 60-char line wrapping."""
    with open(filepath, "w") as fh:
        fh.write(f">{header}\n")
        for i in range(0, len(sequence), 60):
            fh.write(sequence[i : i + 60] + "\n")


def main(input_file: str, output_dir: str) -> None:
    df = pd.read_csv(input_file, sep="\t", header=0)

    # Normalise column names (strip whitespace)
    df.columns = df.columns.str.strip()

    required_cols = {"Ensembl_Gene", "Ensembl_nuc", "symbol"}
    missing = required_cols - set(df.columns)
    if missing:
        log.error("Input file is missing required columns: %s", missing)
        sys.exit(1)

    out_path = Path(output_dir)
    summary_records = []

    for _, row in df.iterrows():
        ensembl_gene = str(row["Ensembl_Gene"]).strip()
        ensembl_nuc = str(row["Ensembl_nuc"]).strip()
        symbol = str(row["symbol"]).strip()

        log.info("Processing %s (%s / %s)", symbol, ensembl_gene, ensembl_nuc)

        record = fetch_transcript(ensembl_nuc, ensembl_gene)

        status = {
            "symbol": symbol,
            "Ensembl_Gene": ensembl_gene,
            "Ensembl_nuc": ensembl_nuc,
            "five_utr_written": False,
            "cds_written": False,
            "three_utr_written": False,
            "error": "",
        }

        if record is None:
            status["error"] = "fetch_failed"
            summary_records.append(status)
            continue

        five_utr = record.get("five_prime_utr_seq", "")
        three_utr = record.get("three_prime_utr_seq", "")
        cds = extract_cds(record)

        transcript_id = record.get("stable_id", ensembl_nuc.split(".")[0])
        version = record.get("stable_id_version", "")
        versioned_id = f"{transcript_id}.{version}" if version else transcript_id
        base_header = f"{versioned_id}|{ensembl_gene}|{symbol}"

        if five_utr:
            fpath = out_path / f"{symbol}_{transcript_id}_5UTR.fa"
            write_fasta(str(fpath), f"{base_header}|5UTR", five_utr)
            log.info("  5'UTR  -> %s (%d nt)", fpath.name, len(five_utr))
            status["five_utr_written"] = True
        else:
            log.warning("  No 5'UTR sequence for %s", symbol)

        if cds:
            fpath = out_path / f"{symbol}_{transcript_id}_CDS.fa"
            write_fasta(str(fpath), f"{base_header}|CDS", cds)
            log.info("  CDS    -> %s (%d nt)", fpath.name, len(cds))
            status["cds_written"] = True
        else:
            log.warning("  Could not derive CDS for %s", symbol)

        if three_utr:
            fpath = out_path / f"{symbol}_{transcript_id}_3UTR.fa"
            write_fasta(str(fpath), f"{base_header}|3UTR", three_utr)
            log.info("  3'UTR  -> %s (%d nt)", fpath.name, len(three_utr))
            status["three_utr_written"] = True
        else:
            log.warning("  No 3'UTR sequence for %s", symbol)

        summary_records.append(status)

        # Polite delay to avoid hammering the API
        time.sleep(0.5)

    # Write summary TSV
    summary_df = pd.DataFrame(summary_records)
    summary_path = out_path / "fetch_summary.tsv"
    summary_df.to_csv(str(summary_path), sep="\t", index=False)
    log.info("Summary written to %s", summary_path)

    ok = summary_df[["five_utr_written", "cds_written", "three_utr_written"]].all(axis=1).sum()
    log.info("Done. %d/%d transcripts fully retrieved.", ok, len(summary_df))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch MANE transcript sequences (5'UTR, CDS, 3'UTR) from Ensembl TARK API"
    )
    parser.add_argument(
        "-i", "--input-file",
        required=True,
        type=str,
        help="Input TSV file containing MANE summary fields (must include Ensembl_Gene, Ensembl_nuc, symbol)",
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        type=str,
        help="Output directory for per-region FASTA files and summary TSV",
    )
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    main(args.input_file, args.output_dir)