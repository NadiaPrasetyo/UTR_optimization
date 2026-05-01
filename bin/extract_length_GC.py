"""
fetch_sequences_Uniprot.py
Command-line tool to fetch protein sequence and metadata from UniProt based on protein data.

Overview:
    - Loads protein records (protein names, gene names, UniProt IDs) from a compiled CSV file.
    - Queries the UniProt API to retrieve full protein information for each protein.
    - Fetches nucleotide (CDS) sequences using the following priority:
        1. Ensembl transcript IDs from UniProt xref cross-references (xref_ensembl).
        2. UniProt ID mapping API (UniProtKB_AC-ID → Ensembl_Transcript) to obtain a transcript ID.
        3. RefSeq nucleotide accession from UniProt cross-references → NCBI efetch.
    - For accessions that had no Ensembl xref after the main pass, attempts a symbol-based
      lookup via the Ensembl REST xrefs/symbol endpoint and retries CDS fetching + CSV update.
    - Parses and standardizes protein metadata including sequence, organism, and nucleotide data.
    - For Ensembl-sourced sequences, separates the transcript into 5' UTR, CDS, and 3' UTR by:
        • Fetching the full cdna with mask_feature=true (introns lowercase, exons uppercase).
        • Fetching the CDS sequence separately.
        • Locating the CDS within the uppercase-only (exonic) portion of the cdna to determine
          UTR boundaries, then mapping those boundaries back to the original masked cdna string.
    - Writes each result immediately to the output CSV so progress is preserved if interrupted.
    - Resumes from where it left off by skipping accessions already present in the output file.

Arguments:
    --org (str): Full name of the organism (used in UniProt queries). Default: "null".
    --ensembl-species (str): Ensembl species name for symbol lookups (e.g. "bos_taurus").
                             Default: "bos_taurus".
    --output (str): Output CSV file path. Default: "data/GC_length_seq.csv".
    --input (str): Required. Input TSV file path with protein data.
    --fasta (bool): If specified, outputs protein data in FASTA format instead of CSV.

Requirements:
    - Input TSV file with 'bovine_gene' and 'bovine_id' columns.
    - Python packages: argparse, csv, requests, os, re, time.

Usage Example:
    python fetch_sequences_Uniprot.py --org "Bos taurus" --ensembl-species bos_taurus \
        --output proteins.csv --input proteins.tsv

Outputs:
    - A CSV file (or FASTA file if --fasta is specified) containing compiled protein metadata,
      including UniProt accession, protein name, sequence, organism, Ensembl/RefSeq nucleotide
      sequence, sequence length, GC content, and (for Ensembl sources) the separated 5' UTR,
      CDS, and 3' UTR sequences.

Author: Nadia
"""

import csv
import os
import re
import time
import requests
import argparse
from requests.adapters import HTTPAdapter, Retry


# ---------------------------------------------------------------------------
# Field names (single source of truth)
# ---------------------------------------------------------------------------

FIELDNAMES = [
    "uniprot_accession",
    "ensembl_id",
    "gene_name",
    "protein_name",
    "sequence",
    "organism_name",
    "nucleotide_sequence",
    "nucleotide_source",
    "length",
    "gc_content",
    "utr5",
    "cds",
    "utr3",
]

# Maps accession → gene_name for accessions that had no Ensembl xref during
# the main pass.  Populated inside parse_uniprot_response / main().
MISSING_ENSEMBL_XREF: dict[str, str] = {}   # accession → gene_name

# ---------------------------------------------------------------------------
# Shared HTTP session with retry logic
# ---------------------------------------------------------------------------

_retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=_retries))

UNIPROT_API      = "https://rest.uniprot.org"
ENSEMBL_REST     = "https://rest.ensembl.org"
POLLING_INTERVAL = 3  # seconds between ID-mapping job status polls


# ---------------------------------------------------------------------------
# Ensembl ID classification helpers
# ---------------------------------------------------------------------------

# Ensembl stable ID patterns by type:
#   Transcript: ENST000…  or species-specific like ENSBTAT000…
#   Protein:    ENSP000…  or ENSBTAP000…
#   Gene:       ENSG000…  or ENSBTAG000…
_TRANSCRIPT_RE = re.compile(r"ENS[A-Z]*T\d+", re.IGNORECASE)
_PROTEIN_RE    = re.compile(r"ENS[A-Z]*P\d+", re.IGNORECASE)
_GENE_RE       = re.compile(r"ENS[A-Z]*G\d+", re.IGNORECASE)


def ensembl_id_type(eid: str) -> str:
    """Return 'transcript', 'protein', 'gene', or 'unknown' for a stable Ensembl ID."""
    base = eid.split(".")[0]
    if _TRANSCRIPT_RE.fullmatch(base):
        return "transcript"
    if _PROTEIN_RE.fullmatch(base):
        return "protein"
    if _GENE_RE.fullmatch(base):
        return "gene"
    return "unknown"


# ---------------------------------------------------------------------------
# UniProt search helpers
# ---------------------------------------------------------------------------

def fetch_uniprot_data(query, retries=3, delay=5):
    """Fetch data from UniProt API with retries.

    Args:
        query (str): Query string for UniProt API.
        retries (int): Number of retry attempts on failure.
        delay (int): Delay in seconds between retries.

    Returns:
        dict | None: JSON response from UniProt API, or None if all retries fail.
    """
    base_url = f"{UNIPROT_API}/uniprotkb/search"
    params = {
        "query": query,
        "fields": "accession,gene_names,protein_name,sequence,xref_refseq,xref_ensembl,organism_name",
        "format": "json",
        "size": 1,
    }

    for attempt in range(retries):
        try:
            response = session.get(base_url, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"  [WARN] UniProt attempt {attempt + 1} failed: {e}")
            if attempt < retries - 1:
                print(f"  Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                print("  [ERROR] All UniProt retries failed.")
                return None


# ---------------------------------------------------------------------------
# NCBI / RefSeq helpers
# ---------------------------------------------------------------------------

def fetch_refseq_nucleotide(refseq_id):
    """Fetch RefSeq nucleotide sequence from NCBI, skipping complete genomes.

    Args:
        refseq_id (str): RefSeq nucleotide accession (e.g. NM_001234).

    Returns:
        str: Nucleotide sequence string, or "" if unavailable / skipped.
    """
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    fetch_url   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    try:
        summary_resp = session.get(
            summary_url,
            params={"db": "nuccore", "id": refseq_id, "retmode": "json"},
            timeout=30,
        )
        summary_resp.raise_for_status()
        summary_data = summary_resp.json()

        uid   = next(iter(summary_data["result"].keys() - {"uids"}), None)
        title = summary_data["result"].get(uid, {}).get("title", "").lower()

        if any(kw in title for kw in ("complete genome", "complete sequence", "whole genome")):
            print(f"  [SKIP] {refseq_id}: complete genome/sequence — skipping.")
            return ""

    except requests.exceptions.RequestException as e:
        print(f"  [WARN] Could not fetch summary for {refseq_id}: {e}")
        return ""

    try:
        response = session.get(
            fetch_url,
            params={"db": "nuccore", "id": refseq_id, "rettype": "fasta", "retmode": "text"},
            timeout=60,
        )
        response.raise_for_status()
        lines = response.text.strip().splitlines()
        return "".join(line for line in lines if not line.startswith(">"))

    except requests.exceptions.RequestException as e:
        print(f"  [WARN] Failed to fetch nucleotide RefSeq {refseq_id}: {e}")
        return ""


# ---------------------------------------------------------------------------
# Ensembl REST helpers
# ---------------------------------------------------------------------------

def extract_transcript_ids_from_xrefs(entry) -> list:
    """Extract Ensembl *transcript* IDs from UniProt xref_ensembl cross-references.

    Protein IDs (ENSBTAP…) and gene IDs (ENSBTAG…) are intentionally excluded:
    requesting type=cds for a protein-stable-ID returns an amino-acid sequence,
    not a coding DNA sequence.

    Args:
        entry (dict): A single UniProt result entry.

    Returns:
        list[str]: Ensembl transcript IDs (may be versioned, e.g. "ENSBTAT00000063874.2").
    """
    transcript_ids = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "Ensembl":
            eid = xref.get("id", "").strip()
            if eid and ensembl_id_type(eid) == "transcript":
                transcript_ids.append(eid)
    return transcript_ids


def fetch_ensembl_sequence(ensembl_id, seq_type="cdna"):
    """Fetch a nucleotide sequence from the Ensembl REST API.

    The version suffix (e.g. ".2") is stripped before the request because the
    Ensembl REST endpoint does not accept versioned IDs.

    mask_feature=true is always sent so that introns are returned in lowercase
    and exons (UTR + CDS) in uppercase — this is required by split_utr_cds().

    Args:
        ensembl_id (str): Ensembl stable ID, optionally versioned.
        seq_type (str): One of "cds", "cdna", "genomic", "protein".

    Returns:
        str: Sequence string (may contain mixed case when seq_type="cdna"),
             or "" if unavailable.
    """
    base_id = ensembl_id.split(".")[0]
    url     = f"{ENSEMBL_REST}/sequence/id/{base_id}"
    params  = {"type": seq_type, "content-type": "text/plain", "mask_feature": "true"}

    try:
        response = session.get(url, params=params, timeout=30)

        if response.status_code == 400:
            if seq_type != "cdna":
                print(f"  [WARN] Ensembl 400 for {base_id} (type={seq_type}); retrying with type=cdna.")
                return fetch_ensembl_sequence(ensembl_id, seq_type="cdna")
            print(f"  [WARN] Ensembl 400 for {base_id}. Skipping.")
            return ""

        response.raise_for_status()
        sequence = response.text.strip()

        if sequence.startswith(">"):
            lines    = sequence.splitlines()
            sequence = "".join(l for l in lines if not l.startswith(">"))

        return sequence

    except requests.exceptions.RequestException as e:
        print(f"  [WARN] Failed to fetch Ensembl sequence for {base_id}: {e}")
        return ""


# ---------------------------------------------------------------------------
# UTR / CDS splitter
# ---------------------------------------------------------------------------

def split_utr_cds(cdna_masked: str, ensembl_id: str) -> tuple[str, str, str]:
    """Separate 5' UTR, CDS, and 3' UTR from a masked cdna sequence.

    The Ensembl REST API returns cdna with ``mask_feature=true`` such that:
      - Intronic / non-coding bases are **lowercase**
      - Exonic bases (UTR + CDS) are **UPPERCASE**

    The CDS is fetched separately (type=cds, also uppercase) and located
    within the exon-only string derived from the cdna.  The position of the
    CDS within the exon string maps back to positions in the full masked cdna,
    allowing precise 5' UTR and 3' UTR extraction.

    Strategy
    --------
    1. Build ``exon_seq`` by concatenating only the uppercase characters of
       ``cdna_masked`` while recording their original indices.
    2. Fetch the CDS (type=cds) for the same transcript — it will be a
       contiguous uppercase string that must appear somewhere in exon_seq.
    3. Find the CDS within exon_seq (case-insensitive, both are all-caps).
    4. Use the exon index map to recover the start/end positions in the
       original masked cdna string.
    5. Slice the original cdna to extract:
         • 5' UTR  = cdna_masked[0 : cds_start_in_cdna]          (may include intronic lowercase)
         • CDS     = exon_seq[cds_start_in_exon : cds_end_in_exon] (uppercase only, no introns)
         • 3' UTR  = cdna_masked[cds_end_in_cdna :]               (may include intronic lowercase)

    Args:
        cdna_masked (str): Full masked cdna string (mixed case) returned by
                           fetch_ensembl_sequence(id, seq_type="cdna").
        ensembl_id (str): The transcript stable ID (used for the CDS fetch).

    Returns:
        tuple[str, str, str]: (utr5, cds, utr3).
            Each is an empty string if the split could not be performed.
    """
    if not cdna_masked:
        return "", "", ""

    # ------------------------------------------------------------------
    # 1. Build exon-only sequence and index map
    # ------------------------------------------------------------------
    exon_indices: list[int] = []   # exon_indices[i] = position in cdna_masked
    exon_seq_chars: list[str] = []

    for i, ch in enumerate(cdna_masked):
        if ch.isupper():
            exon_indices.append(i)
            exon_seq_chars.append(ch)

    exon_seq = "".join(exon_seq_chars)

    if not exon_seq:
        print(f"  [WARN] split_utr_cds: no uppercase (exonic) characters found in cdna for {ensembl_id}.")
        return "", "", ""

    # ------------------------------------------------------------------
    # 2. Fetch the CDS sequence for this transcript
    # ------------------------------------------------------------------
    cds_seq = fetch_ensembl_sequence(ensembl_id, seq_type="cds")

    if not cds_seq:
        print(f"  [WARN] split_utr_cds: could not fetch CDS for {ensembl_id}; UTR split skipped.")
        return "", "", ""

    # Normalise to uppercase for matching (cds is already uppercase, but be safe)
    cds_upper  = cds_seq.upper()
    exon_upper = exon_seq.upper()

    # ------------------------------------------------------------------
    # 3. Locate CDS within the exon sequence
    # ------------------------------------------------------------------
    cds_start_in_exon = exon_upper.find(cds_upper)

    if cds_start_in_exon == -1:
        print(
            f"  [WARN] split_utr_cds: CDS not found in exon sequence for {ensembl_id}. "
            "UTR split skipped."
        )
        return "", "", ""

    cds_end_in_exon = cds_start_in_exon + len(cds_upper)

    # ------------------------------------------------------------------
    # 4. Map exon positions back to original cdna positions
    # ------------------------------------------------------------------
    cds_start_in_cdna = exon_indices[cds_start_in_exon]

    # cds_end_in_exon points one past the last CDS character in exon_seq.
    # The last CDS base is at exon_indices[cds_end_in_exon - 1] in cdna_masked;
    # the cdna slice end is one position after that.
    cds_end_in_cdna = exon_indices[cds_end_in_exon - 1] + 1

    # ------------------------------------------------------------------
    # 5. Slice cdna into UTR / CDS / UTR regions
    # ------------------------------------------------------------------
    utr5 = cdna_masked[:cds_start_in_cdna]
    cds  = cds_seq          # use the clean CDS (no intronic lowercase)
    utr3 = cdna_masked[cds_end_in_cdna:]

    print(
        f"  [INFO] UTR split for {ensembl_id}: "
        f"5'UTR={len(utr5)} nt, CDS={len(cds)} nt, 3'UTR={len(utr3)} nt "
        f"(cdna total={len(cdna_masked)} nt)"
    )

    return utr5, cds, utr3


# ---------------------------------------------------------------------------
# Ensembl symbol → transcript ID lookup
# ---------------------------------------------------------------------------

def fetch_ensembl_id_by_symbol(symbol: str, species: str) -> tuple[list[str], list[str]]:
    """Look up Ensembl IDs for a gene symbol via the xrefs/symbol endpoint.

    Uses:
        GET /xrefs/symbol/:species/:symbol

    Returns transcript IDs and gene IDs separately.  For bovine genes the API
    typically returns only a gene ID (e.g. ENSBTAG); transcript IDs are expanded
    from the gene in a follow-up call if needed.

    Args:
        symbol (str): Gene symbol / display name (e.g. "BRCA2", "SNX11").
        species (str): Ensembl species string (e.g. "bos_taurus", "homo_sapiens").

    Returns:
        tuple[list[str], list[str]]: (transcript_ids, gene_ids)
    """
    url    = f"{ENSEMBL_REST}/xrefs/symbol/{species}/{symbol}"
    params = {"content-type": "application/json"}

    try:
        response = session.get(url, params=params, timeout=30)

        if response.status_code == 400:
            print(f"  [WARN] Ensembl xrefs/symbol 400 for symbol={symbol!r}, species={species!r}.")
            return [], []

        response.raise_for_status()
        hits = response.json()

        transcript_ids = [h["id"] for h in hits if h.get("id") and ensembl_id_type(h["id"]) == "transcript"]
        gene_ids       = [h["id"] for h in hits if h.get("id") and ensembl_id_type(h["id"]) == "gene"]

        if transcript_ids:
            print(f"  [INFO] xrefs/symbol found {len(transcript_ids)} transcript ID(s) for {symbol!r}: {transcript_ids}")
        if gene_ids:
            print(f"  [INFO] xrefs/symbol found {len(gene_ids)} gene ID(s) for {symbol!r}: {gene_ids}")
        if not transcript_ids and not gene_ids:
            print(f"  [INFO] xrefs/symbol returned no usable IDs for {symbol!r}.")

        return transcript_ids, gene_ids

    except requests.exceptions.RequestException as e:
        print(f"  [WARN] xrefs/symbol request failed for {symbol!r}: {e}")
        return [], []


def fetch_transcripts_for_gene(gene_id: str) -> list[str]:
    """Return all transcript stable IDs for an Ensembl gene ID.

    Uses:
        GET /lookup/id/:gene_id?expand=1

    Args:
        gene_id (str): Ensembl gene stable ID (e.g. "ENSBTAG00000020321").

    Returns:
        list[str]: Transcript stable IDs belonging to this gene, or [] on failure.
    """
    base_id = gene_id.split(".")[0]
    url     = f"{ENSEMBL_REST}/lookup/id/{base_id}"
    params  = {"content-type": "application/json", "expand": 1}

    try:
        response = session.get(url, params=params, timeout=30)

        if response.status_code == 400:
            print(f"  [WARN] Ensembl lookup/id 400 for gene {base_id}.")
            return []

        response.raise_for_status()
        data = response.json()

        transcript_ids = [t["id"] for t in data.get("Transcript", []) if t.get("id")]
        print(f"  [INFO] Gene {base_id} has {len(transcript_ids)} transcript(s): {transcript_ids}")
        return transcript_ids

    except requests.exceptions.RequestException as e:
        print(f"  [WARN] Failed to expand transcripts for gene {base_id}: {e}")
        return []


# ---------------------------------------------------------------------------
# Post-pass resolver for accessions that had no Ensembl xref
# ---------------------------------------------------------------------------

def resolve_missing_ensembl_xrefs(
    missing: dict,          # accession → gene_name
    species: str,
    output_path: str,
) -> None:
    """Attempt Ensembl symbol lookups for accessions that had no xref_ensembl hit.

    For each accession in *missing*:
      1. Query ``GET /xrefs/symbol/:species/:gene_name`` for transcript IDs.
      2. Try ``fetch_ensembl_sequence`` on each transcript ID until one succeeds.
      3. Run split_utr_cds() to separate 5' UTR, CDS, and 3' UTR.
      4. If a CDS is found, patch the matching row in the output CSV in-place
         (ensembl_id, nucleotide_sequence, nucleotide_source, length, gc_content,
          utr5, cds, utr3).

    Rows that still cannot be resolved are reported in a final summary.

    Args:
        missing (dict): Mapping of uniprot_accession → gene_name collected during
                        the main processing pass.
        species (str): Ensembl species string for the symbol lookup endpoint
                       (e.g. "bos_taurus").
        output_path (str): Path to the CSV output file to patch.
    """
    if not missing:
        print("\n[ENSEMBL SYMBOL LOOKUP] No accessions to resolve — all had xref hits.")
        return

    print(f"\n[ENSEMBL SYMBOL LOOKUP] Attempting symbol-based lookup for "
          f"{len(missing)} accession(s) with no Ensembl xref …")

    # ------------------------------------------------------------------ #
    # 1. Load the current CSV into memory so we can patch rows            #
    # ------------------------------------------------------------------ #
    try:
        with open(output_path, newline="") as fh:
            rows = list(csv.DictReader(fh))
    except FileNotFoundError:
        print(f"  [ERROR] Output file {output_path!r} not found — cannot patch.")
        return

    # Index rows by accession for O(1) lookup
    row_index: dict[str, dict] = {r["uniprot_accession"]: r for r in rows}

    still_missing: list[str] = []

    for accession, gene_name in missing.items():
        print(f"\n  [ACCESSION] {accession}  gene_name={gene_name!r}")

        if not gene_name:
            print("  [SKIP] No gene name available for symbol lookup.")
            still_missing.append(accession)
            continue

        # ---------------------------------------------------------------- #
        # 2. Symbol lookup → transcript IDs (direct) or gene IDs          #
        # ---------------------------------------------------------------- #
        transcript_ids, gene_ids = fetch_ensembl_id_by_symbol(gene_name, species)

        if not transcript_ids and gene_ids:
            for gid in gene_ids:
                transcript_ids.extend(fetch_transcripts_for_gene(gid))

        if not transcript_ids:
            print(f"  [WARN] No transcript IDs found for {gene_name!r} — skipping.")
            still_missing.append(accession)
            continue

        # ---------------------------------------------------------------- #
        # 3. Try cdna fetch + UTR split for each transcript ID             #
        # ---------------------------------------------------------------- #
        nucleotide_sequence = ""
        ensembl_id_used     = ""
        utr5 = cds = utr3  = ""

        for tid in transcript_ids:
            print(f"  [INFO] Trying cdna for transcript {tid} …")
            cdna_masked = fetch_ensembl_sequence(tid, seq_type="cdna")
            if cdna_masked:
                u5, cd, u3 = split_utr_cds(cdna_masked, tid)
                if cd:
                    nucleotide_sequence = cdna_masked
                    ensembl_id_used     = tid
                    utr5, cds, utr3     = u5, cd, u3
                    break
            # Fallback: try CDS-only if cdna gave nothing useful
            if not nucleotide_sequence:
                print(f"  [INFO] cdna/split failed for {tid}; trying cds …")
                seq = fetch_ensembl_sequence(tid, seq_type="cds")
                if seq:
                    nucleotide_sequence = seq
                    ensembl_id_used     = tid
                    cds                 = seq
                    break

        if not nucleotide_sequence:
            print(f"  [WARN] Could not fetch sequence for any transcript of {gene_name!r}.")
            still_missing.append(accession)
            continue

        # ---------------------------------------------------------------- #
        # 4. Patch the in-memory row                                       #
        # ---------------------------------------------------------------- #
        length, gc_content = extract_length_gc(nucleotide_sequence)
        nucleotide_source  = f"Ensembl:{ensembl_id_used}"

        if accession in row_index:
            row_index[accession]["ensembl_id"]          = ensembl_id_used
            row_index[accession]["nucleotide_sequence"] = nucleotide_sequence
            row_index[accession]["nucleotide_source"]   = nucleotide_source
            row_index[accession]["length"]              = length
            row_index[accession]["gc_content"]          = gc_content
            row_index[accession]["utr5"]                = utr5
            row_index[accession]["cds"]                 = cds
            row_index[accession]["utr3"]                = utr3
            print(f"  [✓] Patched {accession}: source={nucleotide_source}, length={length}")
        else:
            print(f"  [WARN] Accession {accession} not found in CSV rows — skipping patch.")
            still_missing.append(accession)

    # ------------------------------------------------------------------ #
    # 5. Write patched CSV back to disk                                   #
    # ------------------------------------------------------------------ #
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n[ENSEMBL SYMBOL LOOKUP] Patched CSV written to {output_path!r}.")

    if still_missing:
        print(
            f"[ENSEMBL SYMBOL LOOKUP] {len(still_missing)} accession(s) still have no "
            f"nucleotide sequence after symbol lookup:\n  " + ", ".join(still_missing)
        )
    else:
        print("[ENSEMBL SYMBOL LOOKUP] All missing accessions resolved successfully.")


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

def extract_length_gc(sequence):
    """Compute length and GC content of a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.

    Returns:
        tuple[int, float]: (length, gc_content_percent).
    """
    length = len(sequence)
    if length == 0:
        return 0, 0.0
    gc = sequence.upper().count("G") + sequence.upper().count("C")
    return length, round(gc / length * 100, 2)


# ---------------------------------------------------------------------------
# UniProt response parser
# ---------------------------------------------------------------------------

def parse_uniprot_response(data):
    """Parse a UniProt API response, fetching nucleotide sequences.

    Nucleotide sequence priority:
        1. Ensembl transcript IDs from xref_ensembl cross-references → Ensembl REST
           (type=cdna with mask_feature=true), then split into 5'UTR / CDS / 3'UTR.
        2. RefSeq nucleotide accession from xref_refseq → NCBI efetch (no UTR split).

    Args:
        data (dict): JSON response from UniProt API.

    Returns:
        list[dict]: Parsed protein records.
    """
    results = []
    if not data or "results" not in data:
        return results

    for entry in data.get("results", []):
        accession    = entry.get("primaryAccession", "")
        protein_name = (
            entry.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "")
        )
        sequence = entry.get("sequence", {}).get("value", "")
        organism = entry.get("organism", {}).get("scientificName", "")

        nucleotide_sequence = ""
        nucleotide_source   = ""
        ensembl_id_used     = ""
        utr5 = cds = utr3  = ""

        # ------------------------------------------------------------------
        # Strategy 1: xref_ensembl — transcript IDs only
        # ------------------------------------------------------------------
        xref_transcript_ids = extract_transcript_ids_from_xrefs(entry)
        if xref_transcript_ids:
            print(f"  [INFO] Found {len(xref_transcript_ids)} xref transcript ID(s): {xref_transcript_ids}")
        else:
            print("  [INFO] No transcript IDs found in xref_ensembl.")
            if accession:
                MISSING_ENSEMBL_XREF[accession] = ""

        for eid in xref_transcript_ids:
            print(f"  [INFO] Trying Ensembl xref transcript (cdna): {eid}")
            cdna_masked = fetch_ensembl_sequence(eid, seq_type="cdna")
            if cdna_masked:
                u5, cd, u3 = split_utr_cds(cdna_masked, eid)
                if cd:
                    nucleotide_sequence = cdna_masked
                    ensembl_id_used     = eid
                    nucleotide_source   = f"Ensembl:{eid}"
                    utr5, cds, utr3     = u5, cd, u3
                    break
            # If cdna returned nothing or split failed, try cds-only as fallback
            if not nucleotide_sequence:
                print(f"  [INFO] cdna/split failed for {eid}; trying cds …")
                seq = fetch_ensembl_sequence(eid, seq_type="cds")
                if seq:
                    nucleotide_sequence = seq
                    ensembl_id_used     = eid
                    nucleotide_source   = f"Ensembl:{eid}"
                    cds                 = seq
                    break

        # ------------------------------------------------------------------
        # Strategy 2: RefSeq via NCBI (no UTR split available)
        # ------------------------------------------------------------------
        if not nucleotide_sequence:
            refseq_nuc_id = ""
            for xref in entry.get("uniProtKBCrossReferences", []):
                if xref.get("database") == "RefSeq":
                    for prop in xref.get("properties", []):
                        if prop.get("key") == "NucleotideSequenceId":
                            refseq_nuc_id = prop.get("value", "")
                            break
                if refseq_nuc_id:
                    break

            if refseq_nuc_id:
                print(f"  [INFO] Trying RefSeq fallback: {refseq_nuc_id}")
                nucleotide_sequence = fetch_refseq_nucleotide(refseq_nuc_id)
                if nucleotide_sequence:
                    nucleotide_source = f"RefSeq:{refseq_nuc_id}"
                    # UTR separation not possible from RefSeq without annotation;
                    # leave utr5/cds/utr3 empty.

        if not nucleotide_sequence:
            print(f"  [WARN] No nucleotide sequence found for {accession}.")

        length, gc_content = extract_length_gc(nucleotide_sequence)

        results.append({
            "uniprot_accession":   accession,
            "ensembl_id":          ensembl_id_used,
            "gene_name":           "",   # filled in by main()
            "protein_name":        protein_name,
            "sequence":            sequence,
            "organism_name":       organism,
            "nucleotide_sequence": nucleotide_sequence,
            "nucleotide_source":   nucleotide_source,
            "length":              length,
            "gc_content":          gc_content,
            "utr5":                utr5,
            "cds":                 cds,
            "utr3":                utr3,
        })

    return results


# ---------------------------------------------------------------------------
# CSV loader
# ---------------------------------------------------------------------------

def load_protein_records(file_path):
    """Load protein records from a TSV file.

    Expected columns: ``bovine_gene``, ``bovine_id``.

    Args:
        file_path (str): Path to the TSV file.

    Returns:
        list[dict]: Records with keys ``gene_name`` and ``uniprot_id``.
    """
    records = []
    with open(file_path, newline="") as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            records.append({
                "gene_name":  row.get("bovine_gene", "").strip(),
                "uniprot_id": row.get("bovine_id", "").strip(),
            })
    return records


# ---------------------------------------------------------------------------
# Resume helpers
# ---------------------------------------------------------------------------

def load_completed_accessions(output_path):
    """Read the output CSV (if it exists) and return already-processed accessions.

    Also populates ``MISSING_ENSEMBL_XREF`` with any rows that are already written
    but still lack an ``ensembl_id``, so the post-pass resolver can fill them in
    even when every record was processed on a previous run.

    Args:
        output_path (str): Path to the output CSV file.

    Returns:
        tuple[set[str], bool]: (seen_accessions, file_exists)
    """
    if not os.path.exists(output_path):
        return set(), False

    seen = set()
    try:
        with open(output_path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                acc = row.get("uniprot_accession", "").strip()
                if not acc:
                    continue
                seen.add(acc)

                if not row.get("ensembl_id", "").strip():
                    MISSING_ENSEMBL_XREF[acc] = row.get("gene_name", "").strip()

    except Exception as e:
        print(f"[WARN] Could not read existing output file for resume: {e}")
        return set(), False

    n_missing = len(MISSING_ENSEMBL_XREF)
    print(
        f"[RESUME] Found {len(seen)} already-processed accession(s) in {output_path} — will skip these."
    )
    if n_missing:
        print(f"[RESUME] {n_missing} of those row(s) have no ensembl_id and will be resolved after the main pass.")

    return seen, True


def open_output_csv(output_path, file_exists):
    """Open the output CSV for appending, writing the header only when the file is new.

    Args:
        output_path (str): Destination CSV file path.
        file_exists (bool): True if the file already exists with a valid header row.

    Returns:
        tuple[TextIO, csv.DictWriter]: Open file handle and configured DictWriter.
    """
    mode    = "a" if file_exists else "w"
    outfile = open(output_path, mode, newline="")
    writer  = csv.DictWriter(outfile, fieldnames=FIELDNAMES)
    if not file_exists:
        writer.writeheader()
        outfile.flush()
    return outfile, writer


# ---------------------------------------------------------------------------
# FASTA writer
# ---------------------------------------------------------------------------

def write_fasta(protein_data, output_path):
    """Write protein (amino-acid) sequences to a FASTA file.

    Args:
        protein_data (list[dict]): Parsed protein records.
        output_path (str): Destination file path.
    """
    with open(output_path, "w") as f:
        for p in protein_data:
            header = f">{p['uniprot_accession']} {p['protein_name']} [{p['organism_name']}]"
            f.write(header + "\n")
            seq = p["sequence"]
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
    print(f"[DONE] Wrote {len(protein_data)} sequences to FASTA: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(organism="null", ensembl_species="bos_taurus", output=None, input_file=None, fasta=False):
    """Fetch and compile protein data from UniProt.

    Nucleotide sequences are fetched via Ensembl (cdna with mask_feature=true),
    then split into 5' UTR, CDS, and 3' UTR using split_utr_cds().  If Ensembl
    is unavailable, RefSeq is used as a fallback (no UTR split).

    After the main pass, any accessions that had no Ensembl xref are resolved via
    the Ensembl xrefs/symbol endpoint and the output CSV is patched in-place.

    Results are written to the output CSV incrementally so a killed run can be
    resumed by re-running the same command.

    Args:
        organism (str): Full organism name for UniProt queries.
        ensembl_species (str): Ensembl species name for symbol lookups (e.g. "bos_taurus").
        output (str): Output file path (CSV or FASTA).
        input_file (str): Input TSV file path.
        fasta (bool): If True, write FASTA instead of CSV.
    """
    protein_records = load_protein_records(input_file)

    # ------------------------------------------------------------------
    # FASTA mode: collect everything then write (not resumable)
    # ------------------------------------------------------------------
    if fasta:
        protein_data    = []
        seen_accessions = set()

        for record in protein_records:
            gene_name  = record.get("gene_name")
            uniprot_id = record.get("uniprot_id")
            print(f"\n[INFO] Processing gene: {gene_name!r}, UniProt ID: {uniprot_id!r}")

            if uniprot_id:
                query = f"accession:{uniprot_id}"
            elif gene_name:
                query = f"gene:{gene_name}"
            else:
                print("  [SKIP] No gene name or UniProt ID — skipping record.")
                continue

            if organism.lower() != "null":
                query += f' AND organism_name:"{organism}"'

            entry  = fetch_uniprot_data(query)
            parsed = parse_uniprot_response(entry)

            for p in parsed:
                if p["uniprot_accession"] not in seen_accessions:
                    p["gene_name"] = gene_name
                    protein_data.append(p)
                    seen_accessions.add(p["uniprot_accession"])
                    print(f"  [✓] Stored: {p['uniprot_accession']} | source: {p['nucleotide_source'] or 'none'}")
                else:
                    print(f"  [WARN] Duplicate accession {p['uniprot_accession']} — skipped.")

        if not protein_data:
            print("\n[WARN] No matching proteins found. Output file not written.")
            return
        write_fasta(protein_data, output)
        return

    # ------------------------------------------------------------------
    # CSV mode: incremental write with resume support
    # ------------------------------------------------------------------
    seen_accessions, file_exists = load_completed_accessions(output)

    input_ids = {
        r.get("uniprot_id") or r.get("gene_name")
        for r in protein_records
        if r.get("uniprot_id") or r.get("gene_name")
    }
    all_done = file_exists and input_ids and input_ids.issubset(seen_accessions)

    output_dir = os.path.dirname(output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    total_written        = 0
    total_skipped_resume = 0

    if all_done:
        print(
            f"\n[RESUME] All {len(seen_accessions)} accession(s) already processed — "
            "skipping fetch loop."
        )
    else:
        outfile, writer = open_output_csv(output, file_exists)
        try:
            for record in protein_records:
                gene_name  = record.get("gene_name")
                uniprot_id = record.get("uniprot_id")
                print(f"\n[INFO] Processing gene: {gene_name!r}, UniProt ID: {uniprot_id!r}")

                if uniprot_id and uniprot_id in seen_accessions:
                    print(f"  [RESUME] Accession {uniprot_id!r} already in output — skipping.")
                    total_skipped_resume += 1
                    continue

                if uniprot_id:
                    query = f"accession:{uniprot_id}"
                elif gene_name:
                    query = f"gene:{gene_name}"
                else:
                    print("  [SKIP] No gene name or UniProt ID — skipping record.")
                    continue

                if organism.lower() != "null":
                    query += f' AND organism_name:"{organism}"'

                entry  = fetch_uniprot_data(query)
                parsed = parse_uniprot_response(entry)

                for p in parsed:
                    accession = p["uniprot_accession"]
                    if accession in seen_accessions:
                        print(f"  [WARN] Duplicate accession {accession} — skipped.")
                        continue

                    p["gene_name"] = gene_name

                    if accession in MISSING_ENSEMBL_XREF:
                        MISSING_ENSEMBL_XREF[accession] = gene_name

                    writer.writerow(p)
                    outfile.flush()
                    os.fsync(outfile.fileno())

                    seen_accessions.add(accession)
                    total_written += 1
                    print(f"  [✓] Written: {accession} | source: {p['nucleotide_source'] or 'none'}")

        finally:
            outfile.close()

        if total_written == 0 and total_skipped_resume == 0:
            print("\n[WARN] No matching proteins found.")
        else:
            print(
                f"\n[DONE] Finished. "
                f"{total_written} new record(s) written, "
                f"{total_skipped_resume} record(s) skipped (already done). "
                f"Output: {output}"
            )

    # ------------------------------------------------------------------
    # Post-pass: resolve rows that are missing an ensembl_id via symbol
    # ------------------------------------------------------------------
    unresolved_missing: dict[str, str] = {}
    if MISSING_ENSEMBL_XREF:
        try:
            with open(output, newline="") as fh:
                for row in csv.DictReader(fh):
                    acc = row.get("uniprot_accession", "").strip()
                    if acc in MISSING_ENSEMBL_XREF and not row.get("ensembl_id", "").strip():
                        unresolved_missing[acc] = MISSING_ENSEMBL_XREF[acc]
        except FileNotFoundError:
            pass

    resolve_missing_ensembl_xrefs(unresolved_missing, ensembl_species, output)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch protein sequences and metadata from UniProt."
    )
    parser.add_argument(
        "--org",
        default="null",
        help='Organism name for UniProt queries (e.g. "Bos taurus"). Default: "null".',
    )
    parser.add_argument(
        "--ensembl-species",
        default="bos_taurus",
        dest="ensembl_species",
        help=(
            'Ensembl species name used for xrefs/symbol lookups '
            '(e.g. "bos_taurus", "homo_sapiens"). Default: "bos_taurus".'
        ),
    )
    parser.add_argument(
        "--output",
        default="data/GC_length_seq.csv",
        help="Output CSV (or FASTA) file path.",
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input TSV file with 'bovine_gene' and 'bovine_id' columns.",
    )
    parser.add_argument(
        "--fasta",
        action="store_true",
        help="Write output as FASTA instead of CSV.",
    )
    args = parser.parse_args()

    main(
        organism=args.org,
        ensembl_species=args.ensembl_species,
        output=args.output,
        input_file=args.input,
        fasta=args.fasta,
    )