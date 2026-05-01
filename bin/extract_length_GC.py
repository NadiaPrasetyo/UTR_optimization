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
    - Parses and standardizes protein metadata including sequence, organism, and nucleotide data.
    - Writes each result immediately to the output CSV so progress is preserved if interrupted.
    - Resumes from where it left off by skipping accessions already present in the output file.

Arguments:
    --org (str): Full name of the organism (used in UniProt queries). Default: "null".
    --output (str): Output CSV file path. Default: "data/GC_length_seq.csv".
    --input (str): Required. Input TSV file path with protein data.
    --fasta (bool): If specified, outputs protein data in FASTA format instead of CSV.

Requirements:
    - Input TSV file with 'bovine_gene' and 'bovine_id' columns.
    - Python packages: argparse, csv, requests, os, re, time.

Usage Example:
    python fetch_sequences_Uniprot.py --org "Bos taurus" --output proteins.csv --input proteins.tsv

Outputs:
    - A CSV file (or FASTA file if --fasta is specified) containing compiled protein metadata,
      including UniProt accession, protein name, sequence, organism, Ensembl/RefSeq nucleotide
      sequence, sequence length, and GC content.

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
]

MISSING_ENSEMBL_XREF = set()

# ---------------------------------------------------------------------------
# Shared HTTP session with retry logic
# ---------------------------------------------------------------------------

_retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=_retries))

UNIPROT_API      = "https://rest.uniprot.org"
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


def fetch_ensembl_sequence(ensembl_id, seq_type="cds"):
    """Fetch a nucleotide sequence from the Ensembl REST API.

    The version suffix (e.g. ".2") is stripped before the request because the
    Ensembl REST endpoint does not accept versioned IDs.

    Args:
        ensembl_id (str): Ensembl stable ID, optionally versioned.
        seq_type (str): One of "cds", "cdna", "genomic", "protein".

    Returns:
        str: Sequence string, or "" if unavailable.
    """
    base_id = ensembl_id.split(".")[0]
    url     = f"https://rest.ensembl.org/sequence/id/{base_id}"
    params  = {"type": seq_type, "content-type": "text/plain"}

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
    """Parse a UniProt API response, fetching nucleotide CDS sequences.

    Nucleotide sequence priority:
        1. Ensembl transcript IDs from xref_ensembl cross-references → Ensembl REST (type=cds).
        2. UniProt ID mapping API (UniProtKB_AC-ID → Ensembl_Transcript) → Ensembl REST (type=cds).
        3. RefSeq nucleotide accession from xref_refseq → NCBI efetch.

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

        # ------------------------------------------------------------------
        # Strategy 1: xref_ensembl — transcript IDs only
        # ------------------------------------------------------------------
        xref_transcript_ids = extract_transcript_ids_from_xrefs(entry)
        if xref_transcript_ids:
            print(f"  [INFO] Found {len(xref_transcript_ids)} xref transcript ID(s): {xref_transcript_ids}")
        else:
            print("  [INFO] No transcript IDs found in xref_ensembl.")
            if accession:
                MISSING_ENSEMBL_XREF.add(accession)

        for eid in xref_transcript_ids:
            print(f"  [INFO] Trying Ensembl xref transcript: {eid}")
            nucleotide_sequence = fetch_ensembl_sequence(eid, seq_type="cds")
            if nucleotide_sequence:
                ensembl_id_used   = eid
                nucleotide_source = f"Ensembl:{eid}"
                break

        # ------------------------------------------------------------------
        # Strategy 2: RefSeq via NCBI
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
    """Read the output CSV (if it exists) and return the set of already-processed accessions.

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
                if acc:
                    seen.add(acc)
    except Exception as e:
        print(f"[WARN] Could not read existing output file for resume: {e}")
        return set(), False

    print(f"[RESUME] Found {len(seen)} already-processed accession(s) in {output_path} — will skip these.")
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

def main(organism="null", output=None, input_file=None, fasta=False):
    """Fetch and compile protein data from UniProt.

    Nucleotide CDS sequences are fetched via Ensembl transcript IDs (xref first,
    then UniProt ID mapping API), falling back to RefSeq if both Ensembl strategies fail.

    Results are written to the output CSV incrementally so a killed run can be
    resumed by re-running the same command.

    Args:
        organism (str): Full organism name for UniProt queries.
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

    output_dir = os.path.dirname(output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    total_written        = 0
    total_skipped_resume = 0

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

    if MISSING_ENSEMBL_XREF:
        print("\n[SUMMARY] Accessions with NO Ensembl xref:")
        print(",".join(sorted(MISSING_ENSEMBL_XREF)))


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
        output=args.output,
        input_file=args.input,
        fasta=args.fasta,
    )