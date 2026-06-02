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

Ensembl REST API
----------------
Coordinates are fetched from:
    https://rest.ensembl.org/lookup/id/{ENST_ID}?content-type=application/json

Key fields used from the response:
    seq_region_name  → chromosome (e.g. "1" → normalised to "chr1")
    start            → 1-based inclusive start (converted to 0-based for pyBigWig)
    end              → 1-based inclusive end
    strand           → 1 or -1
    display_name     → e.g. "ARV1-201"

PhyloP bigWig (Zoonomia 241-mammalian, hg38)
--------------------------------------------
    https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/
        Homo_sapiens/241-mammalian-2020v2.bigWig

    Chromosomes present: chr1–chr22, chrX, chrY
    Scaffolds / patches (chrUn_*, chr*_random, etc.) are absent — transcripts
    on those regions will be recorded with status "chrom_not_in_bigwig".

Dependencies
------------
    pip install pyBigWig requests pandas numpy tqdm
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

BIGWIG_URL = (
    "https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/"
    "Homo_sapiens/241-mammalian-2020v2.bigWig"
)

ENSEMBL_REST_BASE = "https://rest.ensembl.org"
ENSEMBL_LOOKUP_URL = ENSEMBL_REST_BASE + "/lookup/id/{transcript_id}"
ENSEMBL_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}

# Rate-limiting / retry settings
ENSEMBL_SLEEP_S = 0.1        # polite pause between every request
ENSEMBL_RETRY_MAX = 3        # retries on transient errors
ENSEMBL_RETRY_SLEEP_S = 5    # wait before each retry

# Canonical autosomes + sex chromosomes present in this bigWig
BIGWIG_VALID_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

# ---------------------------------------------------------------------------
# Logging — module-level so all functions share the same logger
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def strip_version(transcript_id: str) -> str:
    """Remove Ensembl version suffix  (ENST00000310256.7 → ENST00000310256)."""
    return transcript_id.split(".")[0]


def ensembl_chrom_to_ucsc(seq_region_name: str) -> str:
    """
    Convert Ensembl seq_region_name to UCSC-style chromosome name.

    Ensembl uses bare numbers and letters ("1", "X", "MT").
    UCSC / this bigWig uses "chr" prefixes ("chr1", "chrX").
    Mitochondrial DNA is "MT" in Ensembl and "chrM" in UCSC — mapped explicitly.
    Unplaced scaffolds are passed through with a "chr" prefix; they will be
    absent from the bigWig and handled gracefully downstream.
    """
    if seq_region_name == "MT":
        return "chrM"
    if seq_region_name.startswith("chr"):
        return seq_region_name          # already UCSC-style
    return f"chr{seq_region_name}"


# ---------------------------------------------------------------------------
# Ensembl coordinate lookup
# ---------------------------------------------------------------------------

def fetch_transcript_coords(transcript_id: str, session: requests.Session) -> dict | None:
    """
    Query the Ensembl REST API and return a dict:
        {chrom, start, end, strand, display_name, ensembl_biotype}

    Coordinates are converted to 0-based half-open (pyBigWig convention):
        bw_start = ensembl_start - 1
        bw_end   = ensembl_end        (unchanged)

    Returns None on unrecoverable failure (404, max retries exceeded, etc.).
    """
    tid = strip_version(transcript_id)
    url = ENSEMBL_LOOKUP_URL.format(transcript_id=tid)

    for attempt in range(1, ENSEMBL_RETRY_MAX + 1):
        try:
            resp = session.get(url, headers=ENSEMBL_HEADERS, timeout=30)

            if resp.status_code == 200:
                data = resp.json()

                # Validate required fields are present
                for field in ("seq_region_name", "start", "end"):
                    if field not in data:
                        log.error("Ensembl response for %s missing field '%s'", tid, field)
                        return None

                chrom = ensembl_chrom_to_ucsc(data["seq_region_name"])
                ensembl_start = int(data["start"])
                ensembl_end   = int(data["end"])

                # Sanity check
                if ensembl_start > ensembl_end:
                    log.warning(
                        "%s: Ensembl returned start > end (%d > %d) — skipping",
                        tid, ensembl_start, ensembl_end,
                    )
                    return None

                return {
                    "chrom":          chrom,
                    "start":          ensembl_start - 1,   # 0-based
                    "end":            ensembl_end,          # half-open
                    "strand":         int(data.get("strand", 1)),
                    "display_name":   data.get("display_name", tid),
                    "ensembl_biotype": data.get("biotype", ""),
                }

            elif resp.status_code == 429:
                retry_after = int(resp.headers.get("Retry-After", ENSEMBL_RETRY_SLEEP_S))
                log.warning("Ensembl rate limit hit — waiting %ds …", retry_after)
                time.sleep(retry_after)
                # Don't count rate-limit waits as an attempt

            elif resp.status_code == 404:
                log.warning("Transcript not found in Ensembl: %s", tid)
                return None

            else:
                log.warning(
                    "Ensembl HTTP %d for %s (attempt %d/%d)",
                    resp.status_code, tid, attempt, ENSEMBL_RETRY_MAX,
                )
                time.sleep(ENSEMBL_RETRY_SLEEP_S)

        except requests.RequestException as exc:
            log.warning("Network error for %s (attempt %d/%d): %s", tid, attempt, ENSEMBL_RETRY_MAX, exc)
            time.sleep(ENSEMBL_RETRY_SLEEP_S)

    log.error("All %d attempts failed for %s", ENSEMBL_RETRY_MAX, tid)
    return None


# ---------------------------------------------------------------------------
# PhyloP score extraction
# ---------------------------------------------------------------------------

def open_bigwig(path_or_url: str) -> pyBigWig.pyBigWig:
    """
    Open a bigWig file from a local path or remote URL.
    Logs the available chromosomes for diagnostic purposes.
    Raises RuntimeError if the file cannot be opened.
    """
    log.info("Opening bigWig: %s", path_or_url)
    bw = pyBigWig.open(path_or_url)
    if bw is None:
        raise RuntimeError(
            f"pyBigWig could not open '{path_or_url}'. "
            "Check the path/URL and that the file is a valid bigWig."
        )
    chroms_in_file = sorted(bw.chroms().keys())
    log.info("Chromosomes: %s", ", ".join(chroms_in_file))
    log.info("bigWig opened successfully.")
    return bw


def get_phylop_stats(
    bw: pyBigWig.pyBigWig,
    chrom: str,
    start: int,
    end: int,
) -> tuple[float | None, float | None, str]:
    """
    Extract per-base PhyloP scores for chrom:start-end (0-based half-open)
    and return (max_phylop, median_phylop, status_detail).

    status_detail is one of:
        "ok"                  — scores computed successfully
        "chrom_not_in_bigwig" — chromosome absent from the file
        "no_data"             — interval exists but all positions are NaN
        "invalid_interval"    — start >= end after clamping

    NaN positions (bases with no conservation score) are excluded from
    both max and median calculations.
    """
    available_chroms = bw.chroms()

    if chrom not in available_chroms:
        log.debug("Chromosome %s not found in bigWig", chrom)
        return None, None, "chrom_not_in_bigwig"

    # Clamp to chromosome boundaries (guards against off-by-one at chr ends)
    chrom_size = available_chroms[chrom]
    start = max(0, start)
    end   = min(end, chrom_size)

    if start >= end:
        log.warning("Interval %s:%d-%d invalid after clamping (chrom size=%d)", chrom, start, end, chrom_size)
        return None, None, "invalid_interval"

    try:
        values = bw.values(chrom, start, end, numpy=True)
    except Exception as exc:
        log.warning("pyBigWig error at %s:%d-%d: %s", chrom, start, end, exc)
        return None, None, "bigwig_error"

    if values is None or len(values) == 0:
        return None, None, "no_data"

    valid = values[~np.isnan(values)]
    if len(valid) == 0:
        log.debug("All NaN scores at %s:%d-%d", chrom, start, end)
        return None, None, "no_data"

    return float(np.max(valid)), float(np.median(valid)), "ok"


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch Zoonomia 241-mammalian PhyloP scores for Ensembl transcript IDs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Input CSV file containing Ensembl transcript IDs.",
    )
    parser.add_argument(
        "--id_col", "-c", default="transcript_id",
        help="Column name that contains Ensembl transcript IDs.",
    )
    parser.add_argument(
        "--output", "-o", default="phylop_results.csv",
        help="Output CSV file path.",
    )
    parser.add_argument(
        "--bigwig", "-b", default=BIGWIG_URL,
        help="Path or URL to the PhyloP bigWig file.",
    )
    parser.add_argument(
        "--sep", default=",",
        help="Delimiter used in the input CSV.",
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Enable DEBUG-level logging.",
    )
    return parser.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # --- Load input --------------------------------------------------------
    input_path = Path(args.input)
    if not input_path.exists():
        log.error("Input file not found: %s", input_path)
        sys.exit(1)

    log.info("Reading input: %s", input_path)
    df = pd.read_csv(input_path, sep=args.sep, dtype=str)

    if args.id_col not in df.columns:
        log.error(
            "Column '%s' not found. Available columns: %s",
            args.id_col, list(df.columns),
        )
        sys.exit(1)

    transcript_ids = df[args.id_col].tolist()
    log.info("Found %d transcript IDs to process", len(transcript_ids))

    # --- Open bigWig -------------------------------------------------------
    bw = open_bigwig(args.bigwig)

    # --- Process each transcript ------------------------------------------
    results = []
    http_session = requests.Session()

    for tid in tqdm(transcript_ids, desc="Transcripts", unit="tx"):
        row: dict = {"transcript_id_query": tid}

        # Step 1: coordinates from Ensembl REST
        coords = fetch_transcript_coords(tid, http_session)
        time.sleep(ENSEMBL_SLEEP_S)

        if coords is None:
            row.update({
                "chrom": None, "start": None, "end": None,
                "strand": None, "display_name": None, "ensembl_biotype": None,
                "max_phylop": None, "median_phylop": None,
                "status": "coord_lookup_failed",
            })
            results.append(row)
            continue

        row.update(coords)

        # Warn if chrom is not in the canonical set (won't be in bigWig)
        if coords["chrom"] not in BIGWIG_VALID_CHROMS:
            log.warning(
                "%s maps to %s — not a canonical chromosome; "
                "PhyloP scores will be unavailable.",
                tid, coords["chrom"],
            )

        # Step 2: PhyloP scores from bigWig
        max_p, med_p, bw_status = get_phylop_stats(
            bw, coords["chrom"], coords["start"], coords["end"]
        )

        row["max_phylop"]    = max_p
        row["median_phylop"] = med_p
        row["status"]        = bw_status

        log.debug(
            "%s  %s:%d-%d  max=%s  median=%s  [%s]",
            tid, coords["chrom"], coords["start"], coords["end"],
            f"{max_p:.4f}" if max_p is not None else "N/A",
            f"{med_p:.4f}" if med_p is not None else "N/A",
            bw_status,
        )

        results.append(row)

    bw.close()

    # --- Merge results onto original dataframe ----------------------------
    results_df = pd.DataFrame(results)
    results_df["_merge_key"] = results_df["transcript_id_query"]

    out_df = df.copy()
    out_df["_merge_key"] = out_df[args.id_col]
    out_df = out_df.merge(
        results_df[[
            "_merge_key", "chrom", "start", "end", "strand",
            "display_name", "ensembl_biotype",
            "max_phylop", "median_phylop", "status",
        ]],
        on="_merge_key",
        how="left",
    ).drop(columns=["_merge_key"])

    # --- Write output -----------------------------------------------------
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_path, index=False)
    log.info("Results written to %s  (%d rows)", output_path, len(out_df))

    # --- Summary ----------------------------------------------------------
    status_counts = out_df["status"].value_counts()
    log.info("Status summary:\n%s", status_counts.to_string())


if __name__ == "__main__":
    main()