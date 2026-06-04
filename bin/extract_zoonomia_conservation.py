#!/usr/bin/env python3
"""
fetch_phylop_scores.py
----------------------
For each Ensembl transcript ID in an input CSV, this script:
  1. Fetches genomic coordinates (chrom, start, end, strand) via the
     Ensembl REST API (GRCh38 / hg38).
  2. BLATs the UTR sequences (utr5, utr3 columns) against hg38 using a
     LOCAL blat binary + hg38.2bit file.
  3. Queries the Zoonomia 241-mammalian PhyloP bigWig for scores over
     the full transcript AND each UTR region.
  4. Writes results to an output CSV, adding:
       utr5_coords, utr3_coords          — chr:start-end (0-based half-open)
       utr5_blat_status, utr3_blat_status
       utr5_max_phylop, utr5_med_phylop
       utr3_max_phylop, utr3_med_phylop
       max_phylop, median_phylop         — whole-transcript scores (unchanged)

Usage
-----
    python fetch_phylop_scores.py --input transcripts.csv \\
                                  --id_col transcript_id \\
                                  --output results.csv \\
                                  --bigwig 241-mammalian-2020v2.bigWig \\
                                  --blat_bin ./blat \\
                                  --twobit hg38.2bit

Input CSV
---------
Must contain:
  • A column with Ensembl transcript IDs (default: transcript_id).
  • A column named  utr5  with the 5′ UTR nucleotide sequence string.
  • A column named  utr3  with the 3′ UTR nucleotide sequence string.
All other columns are preserved in the output.

Local BLAT
----------
Requires:
  • A blat binary (download from https://hgdownload.gi.ucsc.edu/admin/exe/).
  • The hg38.2bit genome file
    (https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit).

For each UTR sequence the script:
  1. Writes the sequence to a temporary FASTA file.
  2. Runs: blat <hg38.2bit> <query.fa> <output.psl> -noHead
  3. Parses the PSL output to pick the best-scoring hit, preferring hits
     on the expected chromosome from the Ensembl lookup.

PSL column indices used (0-based):
  0  matches          8  qNumInsert
  9  tName (chrom)   15  tStart        16  tEnd

Hit selection strategy:
  1. Parse all PSL rows.
  2. Prefer hits on the expected chromosome (from Ensembl lookup).
  3. Among remaining hits, pick the one with the highest score
     (matches − misMatches − qNumInsert).
  4. If no on-chrom hit exists, fall back to the global top hit and flag
     with blat_status="blat_off_chrom".
  5. Very short sequences (<20 nt) are flagged as "blat_seq_too_short"
     and skipped (BLAT is unreliable below this length).

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
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Optional, Tuple

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
ENSEMBL_HEADERS = {
    "Content-Type": "application/json",
    "Accept": "application/json",
}

# Rate-limiting / retry settings for Ensembl
ENSEMBL_SLEEP_S = 0.1
ENSEMBL_RETRY_MAX = 3
ENSEMBL_RETRY_SLEEP_S = 5

BLAT_MIN_SEQ_LEN = 20       # BLAT unreliable below this length

BIGWIG_VALID_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

# PSL column indices (0-based) used during hit parsing
PSL_MATCHES     = 0
PSL_MISMATCHES  = 1
PSL_Q_NUM_INS   = 8
PSL_T_NAME      = 13
PSL_T_START     = 15
PSL_T_END       = 16

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------

def strip_version(transcript_id: str) -> str:
    """Remove Ensembl version suffix (ENST00000310256.7 → ENST00000310256)."""
    return transcript_id.split(".")[0]


def ensembl_chrom_to_ucsc(seq_region_name: str) -> str:
    """Convert Ensembl seq_region_name to UCSC-style chromosome name."""
    if seq_region_name == "MT":
        return "chrM"
    if seq_region_name.startswith("chr"):
        return seq_region_name
    return f"chr{seq_region_name}"


def fmt_coords(chrom: Optional[str], start: Optional[int], end: Optional[int]) -> Optional[str]:
    """Format a coordinate triple as 'chr:start-end' (0-based half-open), or None."""
    if chrom is None or start is None or end is None:
        return None
    return f"{chrom}:{start}-{end}"


# ---------------------------------------------------------------------------
# Ensembl coordinate lookup
# ---------------------------------------------------------------------------

def fetch_transcript_coords(transcript_id: str, session: requests.Session) -> Optional[dict]:
    """
    Query the Ensembl REST API and return a dict with keys:
        chrom, start (0-based), end (half-open), strand, display_name, ensembl_biotype

    Returns None on unrecoverable failure.
    """
    tid = strip_version(transcript_id)
    url = ENSEMBL_LOOKUP_URL.format(transcript_id=tid)

    for attempt in range(1, ENSEMBL_RETRY_MAX + 1):
        try:
            resp = session.get(url, headers=ENSEMBL_HEADERS, timeout=30)

            if resp.status_code == 200:
                data = resp.json()

                for field in ("seq_region_name", "start", "end"):
                    if field not in data:
                        log.error("Ensembl response for %s missing field '%s'", tid, field)
                        return None

                chrom = ensembl_chrom_to_ucsc(data["seq_region_name"])
                ensembl_start = int(data["start"])
                ensembl_end   = int(data["end"])

                if ensembl_start > ensembl_end:
                    log.warning(
                        "%s: Ensembl returned start > end (%d > %d) — skipping",
                        tid, ensembl_start, ensembl_end,
                    )
                    return None

                return {
                    "chrom":           chrom,
                    "start":           ensembl_start - 1,   # convert to 0-based
                    "end":             ensembl_end,          # half-open
                    "strand":          int(data.get("strand", 1)),
                    "display_name":    data.get("display_name", tid),
                    "ensembl_biotype": data.get("biotype", ""),
                }

            elif resp.status_code == 429:
                retry_after = int(resp.headers.get("Retry-After", ENSEMBL_RETRY_SLEEP_S))
                log.warning("Ensembl rate limit — waiting %ds …", retry_after)
                time.sleep(retry_after)

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
# Local BLAT
# ---------------------------------------------------------------------------

def _psl_score(fields: list) -> int:
    """
    Compute a BLAT score from a parsed PSL row (list of strings).

    Formula: matches - misMatches - qNumInsert
    """
    try:
        return (
            int(fields[PSL_MATCHES])
            - int(fields[PSL_MISMATCHES])
            - int(fields[PSL_Q_NUM_INS])
        )
    except (IndexError, ValueError):
        return 0


def blat_sequence(
    sequence: str,
    expected_chrom: Optional[str],
    blat_bin: str,
    twobit: str,
    label: str = "",
) -> Tuple[Optional[str], Optional[int], Optional[int], str]:
    """
    BLAT a nucleotide sequence against hg38 using a local blat binary.

    Parameters
    ----------
    sequence       : nucleotide string (should be ≥ BLAT_MIN_SEQ_LEN nt)
    expected_chrom : chromosome expected from Ensembl lookup (e.g. "chr7");
                     used to prefer on-chromosome hits
    blat_bin       : path to the local blat executable
    twobit         : path to the hg38.2bit genome file
    label          : human-readable label for log messages (e.g. "utr5")

    Returns
    -------
    (chrom, start, end, status)
      chrom/start/end  — 0-based half-open coordinates of the best hit,
                         or (None, None, None) on failure
      status           — one of:
                           "ok"
                           "blat_no_sequence"
                           "blat_seq_too_short"
                           "blat_no_hit"
                           "blat_off_chrom"    (hit found but on wrong chrom)
                           "blat_error"
    """
    seq = str(sequence).strip()

    if not seq or seq.lower() in ("nan", "none", ""):
        return None, None, None, "blat_no_sequence"

    if len(seq) < BLAT_MIN_SEQ_LEN:
        log.debug("%s sequence too short for BLAT (%d nt)", label, len(seq))
        return None, None, None, "blat_seq_too_short"

    with tempfile.TemporaryDirectory() as tmpdir:
        query_fa  = Path(tmpdir) / "query.fa"
        output_psl = Path(tmpdir) / "output.psl"

        # Write query FASTA
        query_fa.write_text(f">query\n{seq}\n")

        # Run local blat
        cmd = [
            blat_bin,
            twobit,
            str(query_fa),
            str(output_psl),
            "-noHead",          # suppress PSL header lines for easy parsing
            "-stepSize=5",      # sensitive setting; matches UCSC web defaults
            "-repMatch=2253",   # hg38 default repeat threshold
            "-minScore=20",
            "-minIdentity=0",
        ]

        log.debug("%s: running blat: %s", label, " ".join(cmd))

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120,
            )
        except subprocess.TimeoutExpired:
            log.warning("%s: blat timed out", label)
            return None, None, None, "blat_error"
        except FileNotFoundError:
            log.error(
                "blat binary not found at '%s'. "
                "Download from https://hgdownload.gi.ucsc.edu/admin/exe/ "
                "and pass the path via --blat_bin.",
                blat_bin,
            )
            raise   # fatal — no point continuing if blat is missing

        if proc.returncode != 0:
            log.warning(
                "%s: blat exited with code %d\nstderr: %s",
                label, proc.returncode, proc.stderr.strip(),
            )
            return None, None, None, "blat_error"

        # Parse PSL rows
        psl_lines = output_psl.read_text().splitlines()

    hits = []
    for line in psl_lines:
        line = line.strip()
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 17:
            continue  # malformed row
        hits.append(fields)

    if not hits:
        log.debug("%s: blat returned no hits", label)
        return None, None, None, "blat_no_hit"

    # Prefer hits on the expected chromosome
    on_chrom  = [h for h in hits if h[PSL_T_NAME] == expected_chrom]
    pool      = on_chrom if on_chrom else hits
    best      = max(pool, key=_psl_score)

    chrom = best[PSL_T_NAME]

    # PSL tStart / tEnd are 0-based half-open (same convention as BED / bigWig)
    try:
        start = int(best[PSL_T_START])
        end   = int(best[PSL_T_END])
    except (IndexError, ValueError) as exc:
        log.warning("%s: could not parse PSL coordinates: %s", label, exc)
        return None, None, None, "blat_error"

    if start >= end:
        log.warning("%s: BLAT hit has start >= end (%d >= %d)", label, start, end)
        return None, None, None, "blat_error"

    status = "ok" if on_chrom else "blat_off_chrom"
    if not on_chrom:
        log.warning(
            "%s: best BLAT hit is on %s, expected %s",
            label, chrom, expected_chrom,
        )

    log.debug(
        "%s BLAT → %s:%d-%d  score=%d  status=%s",
        label, chrom, start, end, _psl_score(best), status,
    )
    return chrom, start, end, status


# ---------------------------------------------------------------------------
# PhyloP score extraction
# ---------------------------------------------------------------------------

def open_bigwig(path: str) -> pyBigWig.pyBigWig:
    """Open a bigWig file from a local path or remote URL."""
    log.info("Opening bigWig: %s", path)
    bw = pyBigWig.open(path)
    if bw is None:
        raise RuntimeError(
            f"pyBigWig could not open '{path}'. "
            "Check the path and that the file is a valid bigWig."
        )
    chroms_in_file = sorted(bw.chroms().keys())
    log.info("Chromosomes in bigWig: %s", ", ".join(chroms_in_file))
    log.info("bigWig opened successfully.")
    return bw


def get_phylop_stats(
    bw: pyBigWig.pyBigWig,
    chrom: Optional[str],
    start: Optional[int],
    end: Optional[int],
) -> Tuple[Optional[float], Optional[float], str]:
    """
    Extract per-base PhyloP scores for chrom:start-end (0-based half-open)
    and return (max_phylop, median_phylop, status_detail).

    Gracefully handles None inputs (returns (None, None, "no_coords")).
    """
    if chrom is None or start is None or end is None:
        return None, None, "no_coords"

    available_chroms = bw.chroms()

    if chrom not in available_chroms:
        log.debug("Chromosome %s not found in bigWig", chrom)
        return None, None, "chrom_not_in_bigwig"

    # Clamp to chromosome boundaries (guards against off-by-one at chr ends)
    chrom_size = available_chroms[chrom]
    start = max(0, start)
    end   = min(end, chrom_size)

    if start >= end:
        log.warning(
            "Interval %s:%d-%d invalid after clamping (chrom size=%d)",
            chrom, start, end, chrom_size,
        )
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
        description=(
            "Fetch Zoonomia 241-mammalian PhyloP scores for Ensembl transcript IDs, "
            "including per-UTR scores derived by BLATting UTR sequences against hg38 "
            "using a local blat binary."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Input CSV file containing Ensembl transcript IDs and UTR sequences.",
    )
    parser.add_argument(
        "--id_col", "-c", default="transcript_id",
        help="Column name that contains Ensembl transcript IDs.",
    )
    parser.add_argument(
        "--utr5_col", default="utr5",
        help="Column name containing 5′ UTR nucleotide sequences.",
    )
    parser.add_argument(
        "--utr3_col", default="utr3",
        help="Column name containing 3′ UTR nucleotide sequences.",
    )
    parser.add_argument(
        "--output", "-o", default="phylop_results.csv",
        help="Output CSV file path.",
    )
    parser.add_argument(
        "--bigwig", "-b", required=True,
        help="Path to the PhyloP bigWig file.",
    )
    parser.add_argument(
        "--blat_bin", default="./blat",
        help=(
            "Path to the local blat executable. "
            "Download from https://hgdownload.gi.ucsc.edu/admin/exe/ "
            "and ensure it is executable (chmod +x blat)."
        ),
    )
    parser.add_argument(
        "--twobit", default="./hg38.2bit",
        help=(
            "Path to the hg38.2bit genome file. "
            "Download from https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit"
        ),
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

    # --- Validate local BLAT resources ------------------------------------
    blat_bin = args.blat_bin
    twobit   = args.twobit

    if not Path(blat_bin).exists():
        log.error(
            "blat binary not found: %s\n"
            "Download from https://hgdownload.gi.ucsc.edu/admin/exe/ "
            "and pass the path via --blat_bin.",
            blat_bin,
        )
        sys.exit(1)

    if not Path(twobit).exists():
        log.error(
            "hg38.2bit file not found: %s\n"
            "Download from https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit "
            "and pass the path via --twobit.",
            twobit,
        )
        sys.exit(1)

    # Quick smoke-test: make sure blat is executable
    try:
        test = subprocess.run(
            [blat_bin],
            capture_output=True, text=True, timeout=10,
        )
        # blat exits non-zero when called with no args but prints usage — that's fine
    except (FileNotFoundError, PermissionError) as exc:
        log.error("Cannot execute blat at '%s': %s", blat_bin, exc)
        sys.exit(1)

    log.info("Using local blat: %s", blat_bin)
    log.info("Using 2bit genome: %s", twobit)

    # --- Load input --------------------------------------------------------
    input_path = Path(args.input)
    if not input_path.exists():
        log.error("Input file not found: %s", input_path)
        sys.exit(1)

    log.info("Reading input: %s", input_path)
    df = pd.read_csv(input_path, sep=args.sep, dtype=str)

    missing = [
        col for col in (args.id_col, args.utr5_col, args.utr3_col)
        if col not in df.columns
    ]
    if missing:
        log.error(
            "Required column(s) not found: %s. Available: %s",
            missing, list(df.columns),
        )
        sys.exit(1)

    transcript_ids = df[args.id_col].tolist()
    utr5_seqs      = df[args.utr5_col].tolist()
    utr3_seqs      = df[args.utr3_col].tolist()
    log.info("Found %d transcript IDs to process", len(transcript_ids))

    # --- Open bigWig -------------------------------------------------------
    bw = open_bigwig(args.bigwig)

    # --- Process each transcript ------------------------------------------
    results = []
    http_session = requests.Session()

    for tid, utr5_seq, utr3_seq in tqdm(
        zip(transcript_ids, utr5_seqs, utr3_seqs),
        total=len(transcript_ids),
        desc="Transcripts",
        unit="tx",
    ):
        row: dict = {"transcript_id_query": tid}

        # ------------------------------------------------------------------
        # Step 1: full-transcript coordinates from Ensembl REST
        # ------------------------------------------------------------------
        coords = fetch_transcript_coords(tid, http_session)
        time.sleep(ENSEMBL_SLEEP_S)

        if coords is None:
            row.update({
                "chrom": None, "start": None, "end": None,
                "strand": None, "display_name": None, "ensembl_biotype": None,
                "max_phylop": None, "median_phylop": None,
                "status": "coord_lookup_failed",
                "utr5_coords": None, "utr5_blat_status": "skipped",
                "utr5_max_phylop": None, "utr5_med_phylop": None,
                "utr3_coords": None, "utr3_blat_status": "skipped",
                "utr3_max_phylop": None, "utr3_med_phylop": None,
            })
            results.append(row)
            continue

        row.update(coords)
        expected_chrom = coords["chrom"]

        if expected_chrom not in BIGWIG_VALID_CHROMS:
            log.warning(
                "%s maps to %s — not a canonical chromosome; "
                "PhyloP scores will be unavailable.",
                tid, expected_chrom,
            )

        # ------------------------------------------------------------------
        # Step 2: whole-transcript PhyloP
        # ------------------------------------------------------------------
        max_p, med_p, bw_status = get_phylop_stats(
            bw, coords["chrom"], coords["start"], coords["end"]
        )

        row["max_phylop"]    = max_p
        row["median_phylop"] = med_p
        row["status"]        = bw_status

        # ------------------------------------------------------------------
        # Step 3: BLAT utr5 → coordinates → PhyloP
        # ------------------------------------------------------------------
        u5_chrom, u5_start, u5_end, u5_blat_status = blat_sequence(
            utr5_seq, expected_chrom,
            blat_bin=blat_bin, twobit=twobit,
            label=f"{tid}/utr5",
        )
        u5_max, u5_med, _ = get_phylop_stats(bw, u5_chrom, u5_start, u5_end)

        row["utr5_coords"]      = fmt_coords(u5_chrom, u5_start, u5_end)
        row["utr5_blat_status"] = u5_blat_status
        row["utr5_max_phylop"]  = u5_max
        row["utr5_med_phylop"]  = u5_med

        # ------------------------------------------------------------------
        # Step 4: BLAT utr3 → coordinates → PhyloP
        # ------------------------------------------------------------------
        u3_chrom, u3_start, u3_end, u3_blat_status = blat_sequence(
            utr3_seq, expected_chrom,
            blat_bin=blat_bin, twobit=twobit,
            label=f"{tid}/utr3",
        )
        u3_max, u3_med, _ = get_phylop_stats(bw, u3_chrom, u3_start, u3_end)

        row["utr3_coords"]      = fmt_coords(u3_chrom, u3_start, u3_end)
        row["utr3_blat_status"] = u3_blat_status
        row["utr3_max_phylop"]  = u3_max
        row["utr3_med_phylop"]  = u3_med

        log.debug(
            "%s  tx=%s:%d-%d [%s]  "
            "utr5=%s [%s] max=%.4g med=%.4g  "
            "utr3=%s [%s] max=%.4g med=%.4g",
            tid,
            coords["chrom"], coords["start"], coords["end"], bw_status,
            row["utr5_coords"], u5_blat_status,
            u5_max if u5_max is not None else float("nan"),
            u5_med if u5_med is not None else float("nan"),
            row["utr3_coords"], u3_blat_status,
            u3_max if u3_max is not None else float("nan"),
            u3_med if u3_med is not None else float("nan"),
        )

        results.append(row)

    bw.close()

    # --- Merge results onto original dataframe ----------------------------
    results_df = pd.DataFrame(results)
    results_df["_merge_key"] = results_df["transcript_id_query"]

    out_df = df.copy()
    out_df["_merge_key"] = out_df[args.id_col]

    merge_cols = [
        "_merge_key",
        "chrom", "start", "end", "strand",
        "display_name", "ensembl_biotype",
        "max_phylop", "median_phylop", "status",
        "utr5_coords", "utr5_blat_status", "utr5_max_phylop", "utr5_med_phylop",
        "utr3_coords", "utr3_blat_status", "utr3_max_phylop", "utr3_med_phylop",
    ]

    out_df = out_df.merge(
        results_df[merge_cols],
        on="_merge_key",
        how="left",
    ).drop(columns=["_merge_key"])

    # --- Write output -----------------------------------------------------
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(output_path, index=False)
    log.info("Results written to %s  (%d rows)", output_path, len(out_df))

    # --- Summary ----------------------------------------------------------
    log.info("=== Transcript status summary ===")
    log.info("\n%s", out_df["status"].value_counts().to_string())

    for utr_label in ("utr5_blat_status", "utr3_blat_status"):
        if utr_label in out_df.columns:
            log.info("=== %s summary ===", utr_label)
            log.info("\n%s", out_df[utr_label].value_counts().to_string())


if __name__ == "__main__":
    main()