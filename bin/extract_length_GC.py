"""
fetch_sequences_Uniprot.py
Command-line tool to fetch protein sequence and metadata from UniProt based on protein data.

Overview:
    - Loads protein records (protein names, gene names, UniProt IDs) from a compiled CSV file.
    - Queries the UniProt API to retrieve full protein information for each protein.
    - Fetches RefSeq nucleotide sequences for matching entries via NCBI, if available.
    - Falls back to Ensembl REST API (using Ensembl IDs from gene evidence) if RefSeq is unavailable.
    - Parses and standardizes protein metadata including sequence, organism, and nucleotide data.
    - Compiles and saves the protein data into a new CSV file for downstream analysis.

Arguments:
    --org (str): Full name of the organism (used in UniProt queries). Default: "null".
    --output (str): Output CSV file path. Default: "data/GC_length_seq.csv".
    --input (str): Required. Input TSV file path with protein data.
    --fasta (bool): If specified, outputs protein data in FASTA format instead of CSV.

Requirements:
    - Input TSV file with 'bovine_gene' and 'bovine_id' columns.
    - Python packages: argparse, csv, requests, os, re, unicodedata, time.

Usage Example:
    python fetch_sequences_Uniprot.py --org "Bos taurus" --output proteins.csv --input proteins.tsv

Outputs:
    - A CSV file (or FASTA file if --fasta is specified) containing compiled protein metadata,
      including UniProt accession, protein name, sequence, organism, RefSeq/Ensembl nucleotide
      sequence, sequence length, and GC content.

Author: Nadia
"""

import csv
import requests
import argparse
import time


# ---------------------------------------------------------------------------
# UniProt helpers
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
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "fields": "accession,gene_names,protein_name,sequence,xref_refseq,xref_ensembl,organism_name",
        "format": "json",
        "size": 1,  # only fetch the top hit
    }

    for attempt in range(retries):
        try:
            response = requests.get(base_url, params=params, timeout=30)
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
    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    # Step 1: check whether this is a whole-genome record
    try:
        summary_resp = requests.get(
            summary_url,
            params={"db": "nuccore", "id": refseq_id, "retmode": "json"},
            timeout=30,
        )
        summary_resp.raise_for_status()
        summary_data = summary_resp.json()

        uid = next(iter(summary_data["result"].keys() - {"uids"}), None)
        title = summary_data["result"].get(uid, {}).get("title", "").lower()

        if any(kw in title for kw in ("complete genome", "complete sequence", "whole genome")):
            print(f"  [SKIP] {refseq_id}: complete genome/sequence — skipping.")
            return ""

    except requests.exceptions.RequestException as e:
        print(f"  [WARN] Could not fetch summary for {refseq_id}: {e}")
        return ""

    # Step 2: fetch FASTA sequence
    try:
        response = requests.get(
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
# Ensembl helpers  (backup)
# ---------------------------------------------------------------------------

def extract_ensembl_id_from_xrefs(entry):
    """Extract the first Ensembl transcript ID from UniProt cross-references (xref_ensembl).

    Args:
        entry (dict): A single UniProt result entry.

    Returns:
        str: First Ensembl transcript ID found, or "" if none.
    """
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "Ensembl":
            eid = xref.get("id", "").strip()
            if eid:
                return eid
    return ""

def extract_ensembl_ids_from_entry(entry):
    """Extract Ensembl transcript/protein IDs from a UniProt entry's gene evidence.

    Looks in ``entry["genes"][*]["geneName"]["evidences"]`` for any evidence
    sourced from Ensembl, returning the associated stable IDs.

    Args:
        entry (dict): A single UniProt result entry.

    Returns:
        list[str]: Ensembl stable IDs (e.g. ["ENSBTAP00000044041.5"]).
    """
    ensembl_ids = []
    for gene in entry.get("genes", []):
        for evidence in gene.get("geneName", {}).get("evidences", []):
            if evidence.get("source", "").lower() == "ensembl":
                eid = evidence.get("id", "").strip()
                if eid:
                    ensembl_ids.append(eid)
    return ensembl_ids


def fetch_ensembl_sequence(ensembl_id, seq_type="cds"):
    """Fetch a nucleotide (or protein) sequence from the Ensembl REST API.

    The version suffix (e.g. ".5") is stripped before the request because the
    Ensembl REST endpoint does not accept versioned IDs.

    Args:
        ensembl_id (str): Ensembl stable ID, optionally versioned (e.g. "ENSBTAP00000044041.5").
        seq_type (str): Sequence type to retrieve. One of:
            - "cds"     – spliced coding sequence without UTR  (default)
            - "cdna"    – spliced transcript sequence with UTR
            - "genomic" – genomic sequence
            - "protein" – translated amino-acid sequence

    Returns:
        str: Sequence string, or "" if unavailable.
    """
    # Strip version suffix if present (e.g. "ENSBTAP00000044041.5" → "ENSBTAP00000044041")
    base_id = ensembl_id.split(".")[0]

    url = f"https://rest.ensembl.org/sequence/id/{base_id}"
    params = {"type": seq_type, "content-type": "text/plain"}

    try:
        response = requests.get(url, params=params, timeout=30)

        if response.status_code == 400:
            # Ensembl may reject the type for this ID class; try "cdna" as fallback
            if seq_type != "cdna":
                print(f"  [WARN] Ensembl returned 400 for {base_id} (type={seq_type}); retrying with type=cdna.")
                return fetch_ensembl_sequence(ensembl_id, seq_type="cdna")
            print(f"  [WARN] Ensembl returned 400 for {base_id}. Skipping.")
            return ""

        response.raise_for_status()
        sequence = response.text.strip()

        # If Ensembl returns FASTA, strip headers
        if sequence.startswith(">"):
            lines = sequence.splitlines()
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
            Returns (0, 0.0) for an empty sequence.
    """
    length = len(sequence)
    if length == 0:
        return 0, 0.0
    gc = sequence.upper().count("G") + sequence.upper().count("C")
    gc_content = gc / length * 100
    return length, round(gc_content, 2)


# ---------------------------------------------------------------------------
# UniProt response parser
# ---------------------------------------------------------------------------

def parse_uniprot_response(data):
    """Parse a UniProt API response, fetching nucleotide sequences where possible.

    Primary strategy: Ensembl ID in gene evidence → Ensembl REST API.
    Backup strategy:  RefSeq cross-reference → NCBI efetch.

    Args:
        data (dict): JSON response from UniProt API.

    Returns:
        list[dict]: Parsed protein records.
    """
    results = []
    if not data or "results" not in data:
        return results

    for entry in data.get("results", []):
        accession = entry.get("primaryAccession", "")
        protein_name = (
            entry.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", "")
        )
        sequence = entry.get("sequence", {}).get("value", "")
        organism = entry.get("organism", {}).get("scientificName", "")

        nucleotide_sequence = ""
        nucleotide_source = ""
        ensembl_id_used = ""

        # Primary: Ensembl — try gene evidence IDs first, then xref cross-references
        ensembl_ids = extract_ensembl_ids_from_entry(entry)
        if not ensembl_ids:
            print("  [INFO] No Ensembl IDs in gene evidence — trying xref_ensembl cross-references.")
            xref_eid = extract_ensembl_id_from_xrefs(entry)
            if xref_eid:
                ensembl_ids = [xref_eid]  # wrap in list so the loop below works unchanged

        for eid in ensembl_ids:
            print(f"  [INFO] Trying Ensembl (primary): {eid}")
            nucleotide_sequence = fetch_ensembl_sequence(eid, seq_type="cds")
            if nucleotide_sequence:
                ensembl_id_used = eid
                nucleotide_source = f"Ensembl:{eid}"
                break

        # Backup: RefSeq
        if not nucleotide_sequence:
            nucleotide_id = ""
            for xref in entry.get("uniProtKBCrossReferences", []):
                if xref.get("database") == "RefSeq":
                    for prop in xref.get("properties", []):
                        if prop.get("key") == "NucleotideSequenceId":
                            nucleotide_id = prop.get("value", "")
                            break
                if nucleotide_id:
                    break
            if nucleotide_id:
                print(f"  [INFO] Ensembl unavailable — trying RefSeq backup: {nucleotide_id}")
                nucleotide_sequence = fetch_refseq_nucleotide(nucleotide_id)
                if nucleotide_sequence:
                    nucleotide_source = f"RefSeq:{nucleotide_id}"

        if not nucleotide_sequence:
            print(f"  [WARN] No nucleotide sequence found for {accession}.")

        length, gc_content = extract_length_gc(nucleotide_sequence)

        results.append({
            "uniprot_accession": accession,
            "ensembl_id": ensembl_id_used,
            "gene_name": "",          # filled in by main()
            "protein_name": protein_name,
            "sequence": sequence,
            "organism_name": organism,
            "nucleotide_sequence": nucleotide_sequence,
            "nucleotide_source": nucleotide_source,
            "length": length,
            "gc_content": gc_content,
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
                "gene_name": row.get("bovine_gene", "").strip(),
                "uniprot_id": row.get("bovine_id", "").strip(),
            })
    return records


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
            header = (
                f">{p['uniprot_accession']} "
                f"{p['protein_name']} "
                f"[{p['organism_name']}]"
            )
            f.write(header + "\n")
            # Wrap sequence at 60 characters
            seq = p["sequence"]
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
    print(f"[DONE] Wrote {len(protein_data)} sequences to FASTA: {output_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(organism="null", output=None, input_file=None, fasta=False):
    """Fetch and compile protein data from UniProt (with Ensembl backup).

    Args:
        organism (str): Full organism name for UniProt queries.
        output (str): Output file path (CSV or FASTA).
        input_file (str): Input TSV file path.
        fasta (bool): If True, write FASTA instead of CSV.
    """
    protein_records = load_protein_records(input_file)
    protein_data = []
    seen_accessions = set()

    for record in protein_records:
        gene_name = record.get("gene_name")
        uniprot_id = record.get("uniprot_id")

        print(f"\n[INFO] Processing gene: {gene_name!r}, UniProt ID: {uniprot_id!r}")

        # Build query
        if uniprot_id:
            query = f"accession:{uniprot_id}"
        elif gene_name:
            query = f"gene:{gene_name}"
        else:
            print("  [SKIP] No gene name or UniProt ID — skipping record.")
            continue

        if organism.lower() != "null":
            query += f" AND organism_name:\"{organism}\""

        entry = fetch_uniprot_data(query)
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

    if fasta:
        write_fasta(protein_data, output)
    else:
        fieldnames = [
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
        with open(output, "w", newline="") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(protein_data)
        print(f"\n[DONE] Wrote {len(protein_data)} proteins to: {output}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch protein sequences and metadata from UniProt (with Ensembl backup)."
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