#!/usr/bin/env python3
"""
fetch_phylop_scores.py
----------------------
For each Ensembl transcript ID in an input CSV, this script:
  1. Fetches genomic coordinates (chrom, start, end, strand) via the
     Ensembl REST API (GRCh38 / hg38).
  2. Queries the Zoonomia 241-mammalian PhyloP bigWig for scores over
     those coordinates.
  3. Records max_phylop and median_phylop per transcript and writes
     results to an output CSV.

Usage
-----
    python fetch_phylop_scores.py --input transcripts.csv \
                                  --id_col transcript_id \
                                  --output results.csv

Input CSV
---------
Must contain at least one column with Ensembl transcript IDs
(e.g. ENST00000310256 or ENST00000310256.7 — version suffixes are stripped).
All other columns are preserved in the output.

Example input (the format present in the provided data):
    gene,gene_id,transcript_id,...
    ARV1,ENSG00000173409.16,ENST00000310256.7,...

Dependencies
------------
    pip install pyBigWig requests pandas numpy tqdm

PhyloP bigWig (Zoonomia 241-mammalian, hg38)
--------------------------------------------
    https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/
        Homo_sapiens/241-mammalian-2020v2.bigWig
"""

import argparse
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig
import requests
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ENSEMBL_REST_BASE = "https://rest.ensembl.org"
ENSEMBL_LOOKUP_URL = ENSEMBL_REST_BASE + "/lookup/id/{transcript_id}"

# Ensembl REST imposes a rate limit; pause between requests to be polite
ENSEMBL_SLEEP_S = 0.1          # seconds between requests
ENSEMBL_RETRY_MAX = 3           # retries on transient HTTP errors
ENSEMBL_RETRY_SLEEP_S = 5       # seconds to wait before retrying

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


# ---------------------------------------------------------------------------
# Ensembl coordinate lookup
# ---------------------------------------------------------------------------

def strip_version(transcript_id: str) -> str:
    """Remove Ensembl version suffix (e.g. ENST00000310256.7 → ENST00000310256)."""
    return transcript_id.split(".")[0]


def fetch_transcript_coords(transcript_id: str, session: requests.Session) -> dict | None:
    """
    Query the Ensembl REST API and return a dict with keys:
        chrom, start, end, strand, display_name
    Returns None if the lookup fails.

    Ensembl coordinates are 1-based inclusive; pyBigWig uses 0-based half-open,
    so we convert: bw_start = ensembl_start - 1, bw_end = ensembl_end.
    """
    tid = strip_version(transcript_id)
    url = ENSEMBL_LOOKUP_URL.format(transcript_id=tid)
    headers = {"Content-Type": "application/json"}

    for attempt in range(1, ENSEMBL_RETRY_MAX + 1):
        try:
            resp = session.get(url, headers=headers, timeout=20)

            if resp.status_code == 200:
                data = resp.json()
                chrom = data.get("seq_region_name", "")
                # Normalise chromosome name to UCSC style (chr1, chrX …)
                if not chrom.startswith("chr"):
                    chrom = f"chr{chrom}"
                return {
                    "chrom": chrom,
                    # Convert 1-based Ensembl → 0-based for pyBigWig
                    "start": int(data["start"]) - 1,
                    "end": int(data["end"]),
                    "strand": data.get("strand", 1),
                    "display_name": data.get("display_name", tid),
                }

            elif resp.status_code == 429:
                # Rate-limited — honour the Retry-After header if present
                retry_after = int(resp.headers.get("Retry-After", ENSEMBL_RETRY_SLEEP_S))
                logging.warning("Rate-limited by Ensembl; sleeping %ds …", retry_after)
                time.sleep(retry_after)

            elif resp.status_code == 404:
                logging.warning("Transcript not found: %s", tid)
                return None

            else:
                logging.warning(
                    "Ensembl HTTP %d for %s (attempt %d/%d)",
                    resp.status_code, tid, attempt, ENSEMBL_RETRY_MAX,
                )
                time.sleep(ENSEMBL_RETRY_SLEEP_S)

        except requests.RequestException as exc:
            logging.warning("Request error for %s (attempt %d): %s", tid, attempt, exc)
            time.sleep(ENSEMBL_RETRY_SLEEP_S)

    logging.error("All %d attempts failed for %s", ENSEMBL_RETRY_MAX, tid)
    return None


# ---------------------------------------------------------------------------
# PhyloP score extraction
# ---------------------------------------------------------------------------

def open_bigwig(file: str) -> pyBigWig.pyBigWig:
    """Open a local file or remote URL as a pyBigWig object."""
    try:
        bw = pyBigWig.open(file)
        if bw is None:
            raise RuntimeError("pyBigWig returned None — check the URL/path.")
        return bw
    except Exception as exc:
        logging.error("Cannot open bigWig at %s: %s", file, exc)
        raise


def get_phylop_stats(
    bw: pyBigWig.pyBigWig,
    chrom: str,
    start: int,
    end: int,
) -> tuple[float | None, float | None]:
    """
    Fetch per-base PhyloP values for chrom:start-end (0-based, half-open)
    and return (max_phylop, median_phylop).

    Returns (None, None) when:
    - the chromosome is absent from the bigWig
    - no data exists in the interval
    - the interval is invalid (start >= end)
    """
    if start >= end:
        logging.warning("Invalid interval %s:%d-%d (start >= end)", chrom, start, end)
        return None, None

    available_chroms = bw.chroms()
    if chrom not in available_chroms:
        logging.warning("Chromosome %s not found in bigWig", chrom)
        return None, None

    chrom_size = available_chroms[chrom]
    # Clamp to chromosome boundaries
    start = max(0, start)
    end = min(end, chrom_size)

    try:
        values = bw.values(chrom, start, end, numpy=True)   # returns np.ndarray
    except Exception as exc:
        logging.warning("pyBigWig error for %s:%d-%d: %s", chrom, start, end, exc)
        return None, None

    if values is None or len(values) == 0:
        return None, None

    # Filter out NaN positions (no data)
    valid = values[~np.isnan(values)]
    if len(valid) == 0:
        logging.info("No valid PhyloP data for %s:%d-%d", chrom, start, end)
        return None, None

    return float(np.max(valid)), float(np.median(valid))


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def main(args):

    # --- Load input ---
    input_path = Path(args.input)
    if not input_path.exists():
        logging.error("Input file not found: %s", input_path)
        sys.exit(1)

    logging.info("Reading input: %s", input_path)
    df = pd.read_csv(input_path, sep=args.sep, dtype=str)

    if args.id_col not in df.columns:
        logging.error(
            "Column '%s' not found. Available columns: %s",
            args.id_col, list(df.columns),
        )
        sys.exit(1)

    transcript_ids = df[args.id_col].tolist()
    logging.info("Found %d transcript IDs to process", len(transcript_ids))

    # --- Open bigWig ---
    logging.info("Opening bigWig: %s", args.bigwig)
    bw = open_bigwig(args.bigwig)
    logging.info("bigWig opened successfully.")

    # --- Process each transcript ---
    results = []
    http_session = requests.Session()

    for tid in tqdm(transcript_ids, desc="Transcripts", unit="tx"):
        row: dict = {"transcript_id_query": tid}

        # 1. Get coordinates
        coords = fetch_transcript_coords(tid, http_session)
        time.sleep(ENSEMBL_SLEEP_S)   # be polite to Ensembl

        if coords is None:
            row.update({
                "chrom": None, "start": None, "end": None, "strand": None,
                "max_phylop": None, "median_phylop": None,
                "status": "coord_lookup_failed",
            })
            results.append(row)
            continue

        row.update(coords)

        # 2. Query bigWig
        max_p, med_p = get_phylop_stats(bw, coords["chrom"], coords["start"], coords["end"])

        row["max_phylop"] = max_p
        row["median_phylop"] = med_p
        row["status"] = "ok" if max_p is not None else "no_bigwig_data"
        results.append(row)

        logging.debug(
            "%s  %s:%d-%d  max=%.4f  median=%.4f",
            tid, coords["chrom"], coords["start"], coords["end"],
            max_p or float("nan"), med_p or float("nan"),
        )

    bw.close()

    # --- Merge results back onto original dataframe ---
    results_df = pd.DataFrame(results)

    # Align on the original id column
    df = df.copy()
    df["_merge_key"] = df[args.id_col]
    results_df["_merge_key"] = results_df["transcript_id_query"]

    out_df = df.merge(
        results_df[["_merge_key", "chrom", "start", "end", "strand",
                    "max_phylop", "median_phylop", "status"]],
        on="_merge_key",
        how="left",
    ).drop(columns=["_merge_key"])

    # --- Write output ---
    output_path = Path(args.output)
    out_df.to_csv(output_path, index=False)
    logging.info("Results written to %s  (%d rows)", output_path, len(out_df))

    # Summary
    ok = (out_df["status"] == "ok").sum()
    fail = len(out_df) - ok
    logging.info("Done. Success: %d | Failed/no-data: %d", ok, fail)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch Zoonomia PhyloP scores for Ensembl transcript IDs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Input CSV file path. Must contain a column with Ensembl transcript IDs.",
    )
    parser.add_argument(
        "--id_col", "-c", default="transcript_id",
        help="Name of the column containing Ensembl transcript IDs.",
    )
    parser.add_argument(
        "--output", "-o", default="phylop_results.csv",
        help="Output CSV file path.",
    )
    parser.add_argument(
        "--bigwig", "-b", required=True,
        help="local path to the PhyloP bigWig file (download from https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/Homo_sapiens/241-mammalian-2020v2.bigWig)",
    )
    parser.add_argument(
        "--sep", default=",",
        help="CSV delimiter for the input file.",
    )

    setup_logging()
    main(args=parser.parse_args())