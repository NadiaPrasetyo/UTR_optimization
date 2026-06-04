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
ENSEMBL_LOOKUP_BATCH_URL = (
    ENSEMBL_REST_BASE
    + "/lookup/id"
)
ENSEMBL_BATCH_SIZE = 100
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
def fetch_transcript_coords_batch(
    transcript_ids: list[str],
    session: requests.Session,
) -> dict:
    """
    Returns:

        {
            ENST... : coords_dict,
            ENST... : coords_dict,
        }

    Missing transcripts are omitted.
    """

    payload = {
        "ids": [
            strip_version(tid)
            for tid in transcript_ids
        ]
    }

    params = {
        "expand": 1,
    }

    for attempt in range(1, ENSEMBL_RETRY_MAX + 1):

        try:

            resp = session.post(
                ENSEMBL_LOOKUP_BATCH_URL,
                params=params,
                json=payload,
                headers=ENSEMBL_HEADERS,
                timeout=120,
            )

            if resp.status_code == 200:

                raw = resp.json()

                results = {}

                for tid, data in raw.items():

                    if data is None:
                        continue

                    try:

                        chrom = ensembl_chrom_to_ucsc(
                            data["seq_region_name"]
                        )

                        ensembl_start = int(data["start"])
                        ensembl_end = int(data["end"])

                        translation = data.get("Translation")

                        cds_start = None
                        cds_end = None

                        if translation:

                            cds_start = int(
                                translation["start"]
                            )

                            cds_end = int(
                                translation["end"]
                            )

                        exons = []

                        for exon in data.get("Exon", []):

                            exons.append({
                                "start": int(exon["start"]),
                                "end": int(exon["end"]),
                            })

                        results[tid] = {
                            "chrom": chrom,
                            "start": ensembl_start - 1,
                            "end": ensembl_end,
                            "strand": int(
                                data.get("strand", 1)
                            ),
                            "display_name": data.get(
                                "display_name",
                                tid,
                            ),
                            "ensembl_biotype": data.get(
                                "biotype",
                                "",
                            ),
                            "cds_start": cds_start,
                            "cds_end": cds_end,
                            "exons": exons,
                        }

                    except Exception as exc:

                        log.warning(
                            "Failed parsing %s: %s",
                            tid,
                            exc,
                        )

                return results

            elif resp.status_code == 429:

                retry_after = int(
                    resp.headers.get(
                        "Retry-After",
                        ENSEMBL_RETRY_SLEEP_S,
                    )
                )

                log.warning(
                    "Rate limited, waiting %ds",
                    retry_after,
                )

                time.sleep(retry_after)

            else:

                log.warning(
                    "Ensembl batch lookup returned HTTP %d",
                    resp.status_code,
                )

                time.sleep(
                    ENSEMBL_RETRY_SLEEP_S
                )

        except requests.RequestException as exc:

            log.warning(
                "Batch lookup failed (attempt %d/%d): %s",
                attempt,
                ENSEMBL_RETRY_MAX,
                exc,
            )

            time.sleep(
                ENSEMBL_RETRY_SLEEP_S
            )

    return {}


# ---------------------------------------------------------------------------
# PhyloP score extraction
# ---------------------------------------------------------------------------
def calculate_utr_ranges(
    tx_start,
    tx_end,
    cds_start,
    cds_end,
    strand,
):
    if cds_start is None or cds_end is None:
        return None, None

    if strand == 1:

        utr5 = (
            tx_start,
            cds_start - 1,
        )

        utr3 = (
            cds_end + 1,
            tx_end,
        )

    else:

        utr3 = (
            tx_start,
            cds_start - 1,
        )

        utr5 = (
            cds_end + 1,
            tx_end,
        )

    if utr5[0] > utr5[1]:
        utr5 = None

    if utr3[0] > utr3[1]:
        utr3 = None

    return utr5, utr3


def intersect_exons_with_range(
    exons,
    region_start,
    region_end,
):
    fragments = []

    for exon in exons:

        start = max(exon["start"], region_start)
        end = min(exon["end"], region_end)

        if start <= end:
            fragments.append((start, end))

    return fragments


def get_utr_fragments(
    exons,
    tx_start,
    tx_end,
    cds_start,
    cds_end,
    strand,
):
    utr5_range, utr3_range = calculate_utr_ranges(
        tx_start,
        tx_end,
        cds_start,
        cds_end,
        strand,
    )

    utr5_fragments = []
    utr3_fragments = []

    if utr5_range is not None:
        utr5_fragments = intersect_exons_with_range(
            exons,
            utr5_range[0],
            utr5_range[1],
        )

    if utr3_range is not None:
        utr3_fragments = intersect_exons_with_range(
            exons,
            utr3_range[0],
            utr3_range[1],
        )

    return utr5_fragments, utr3_fragments


def get_fragmented_phylop_stats(
    bw,
    chrom,
    fragments,
):
    if not fragments:
        return None, None, "no_utr"

    chroms = bw.chroms()

    if chrom not in chroms:
        return None, None, "chrom_not_in_bigwig"

    chrom_size = chroms[chrom]

    all_scores = []

    for start_1based, end_1based in fragments:

        start = max(0, start_1based - 1)
        end = min(end_1based, chrom_size)

        if start >= end:
            continue

        try:
            vals = bw.values(
                chrom,
                start,
                end,
                numpy=True,
            )
        except Exception:
            continue

        if vals is None:
            continue

        vals = vals[~np.isnan(vals)]

        if len(vals):
            all_scores.append(vals)

    if not all_scores:
        return None, None, "no_data"

    all_scores = np.concatenate(all_scores)

    return (
        float(np.max(all_scores)),
        float(np.median(all_scores)),
        "ok",
    )

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
def chunks(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i + size]

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

    # ------------------------------------------------------------------
    # Load input
    # ------------------------------------------------------------------

    input_path = Path(args.input)

    if not input_path.exists():
        log.error("Input file not found: %s", input_path)
        sys.exit(1)

    log.info("Reading input: %s", input_path)

    df = pd.read_csv(
        input_path,
        sep=args.sep,
        dtype=str,
    )

    if args.id_col not in df.columns:
        log.error(
            "Column '%s' not found. Available columns: %s",
            args.id_col,
            list(df.columns),
        )
        sys.exit(1)

    transcript_ids = df[args.id_col].tolist()

    log.info(
        "Found %d transcript IDs",
        len(transcript_ids),
    )

    # ------------------------------------------------------------------
    # Open PhyloP bigWig
    # ------------------------------------------------------------------

    bw = open_bigwig(args.bigwig)

    # ------------------------------------------------------------------
    # Create HTTP session
    # ------------------------------------------------------------------

    http_session = requests.Session()

    # ------------------------------------------------------------------
    # Batch Ensembl lookup
    # ------------------------------------------------------------------

    coord_lookup = {}

    batches = list(
        chunks(
            transcript_ids,
            ENSEMBL_BATCH_SIZE,
        )
    )

    log.info(
        "Fetching transcript annotations from Ensembl "
        "(%d batches)",
        len(batches),
    )

    for batch in tqdm(
        batches,
        desc="Ensembl lookup",
        unit="batch",
    ):
        log.info(
            "Submitting Ensembl batch of %d transcripts",
            len(batch),
        )

        batch_results = fetch_transcript_coords_batch(
            batch,
            http_session,
        )

        coord_lookup.update(
            batch_results
        )

        time.sleep(
            ENSEMBL_SLEEP_S
        )

    log.info(
        "Retrieved annotations for %d transcripts",
        len(coord_lookup),
    )

    # ------------------------------------------------------------------
    # Process transcripts
    # ------------------------------------------------------------------

    results = []

    for tid in tqdm(
        transcript_ids,
        desc="Transcripts",
        unit="tx",
    ):

        lookup_id = strip_version(tid)

        row = {
            "transcript_id_query": tid
        }

        coords = coord_lookup.get(
            lookup_id
        )

        # --------------------------------------------------------------
        # Annotation lookup failed
        # --------------------------------------------------------------

        if coords is None:

            row.update({
                "chrom": None,
                "start": None,
                "end": None,
                "strand": None,

                "display_name": None,
                "ensembl_biotype": None,

                "cds_start": None,
                "cds_end": None,

                "utr5_exons": None,
                "utr3_exons": None,

                "max_phylop": None,
                "median_phylop": None,

                "utr5_max_phylop": None,
                "utr5_med_phylop": None,

                "utr3_max_phylop": None,
                "utr3_med_phylop": None,

                "utr5_status": None,
                "utr3_status": None,

                "status": "coord_lookup_failed",
            })

            results.append(row)
            continue

        row.update(coords)

        # --------------------------------------------------------------
        # Determine exon-only UTR fragments
        # --------------------------------------------------------------

        tx_start = coords["start"] + 1
        tx_end = coords["end"]

        utr5_fragments, utr3_fragments = (
            get_utr_fragments(
                exons=coords["exons"],
                tx_start=tx_start,
                tx_end=tx_end,
                cds_start=coords["cds_start"],
                cds_end=coords["cds_end"],
                strand=coords["strand"],
            )
        )

        row["utr5_exons"] = ";".join(
            f"{start}-{end}"
            for start, end in utr5_fragments
        )

        row["utr3_exons"] = ";".join(
            f"{start}-{end}"
            for start, end in utr3_fragments
        )

        # --------------------------------------------------------------
        # Warn if chromosome absent from bigWig
        # --------------------------------------------------------------

        if coords["chrom"] not in BIGWIG_VALID_CHROMS:

            log.warning(
                "%s maps to %s — chromosome absent from PhyloP bigWig",
                tid,
                coords["chrom"],
            )

        # --------------------------------------------------------------
        # Whole transcript PhyloP
        # --------------------------------------------------------------

        max_p, med_p, bw_status = (
            get_phylop_stats(
                bw,
                coords["chrom"],
                coords["start"],
                coords["end"],
            )
        )

        row["max_phylop"] = max_p
        row["median_phylop"] = med_p
        row["status"] = bw_status

        # --------------------------------------------------------------
        # 5' UTR PhyloP
        # --------------------------------------------------------------

        utr5_max, utr5_med, utr5_status = (
            get_fragmented_phylop_stats(
                bw,
                coords["chrom"],
                utr5_fragments,
            )
        )

        row["utr5_max_phylop"] = utr5_max
        row["utr5_med_phylop"] = utr5_med
        row["utr5_status"] = utr5_status

        # --------------------------------------------------------------
        # 3' UTR PhyloP
        # --------------------------------------------------------------

        utr3_max, utr3_med, utr3_status = (
            get_fragmented_phylop_stats(
                bw,
                coords["chrom"],
                utr3_fragments,
            )
        )

        row["utr3_max_phylop"] = utr3_max
        row["utr3_med_phylop"] = utr3_med
        row["utr3_status"] = utr3_status

        results.append(row)

    bw.close()

    # ------------------------------------------------------------------
    # Merge results onto original dataframe
    # ------------------------------------------------------------------

    results_df = pd.DataFrame(results)

    results_df["_merge_key"] = (
        results_df["transcript_id_query"]
    )

    out_df = df.copy()

    out_df["_merge_key"] = (
        out_df[args.id_col]
    )

    out_df = out_df.merge(
        results_df[[
            "_merge_key",

            "chrom",
            "start",
            "end",
            "strand",

            "display_name",
            "ensembl_biotype",

            "cds_start",
            "cds_end",

            "utr5_exons",
            "utr3_exons",

            "max_phylop",
            "median_phylop",

            "utr5_max_phylop",
            "utr5_med_phylop",

            "utr3_max_phylop",
            "utr3_med_phylop",

            "utr5_status",
            "utr3_status",

            "status",
        ]],
        on="_merge_key",
        how="left",
    ).drop(
        columns=["_merge_key"]
    )

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