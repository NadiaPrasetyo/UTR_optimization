#!/usr/bin/env python3
"""
validate_est.py
───────────────
Validates CDS/UTR sequences from a gene CSV against bovine EST data.

Pipeline:
  1. Load gene CSV  (ensembl_gene_id, cds, utr5, utr3, bovine_id, …)
  2. Load refGene PSL/table  → map accession → chrom, strand, coords
  3. Load all_EST PSL table  → index EST alignments by chrom region
  4. For each gene:
       a. Look up its refGene entry to get genomic coordinates
       b. Find overlapping ESTs (from all_EST)
       c. Fetch the RefSeq mRNA sequence from NCBI (eutils REST API)
          – skips complete genomes / whole-genome sequences
       d. Concatenate utr5 + cds + utr3 and align the full mRNA
          against the NCBI sequence
       e. Report % identity and PASS / FAIL (threshold ≥ 96 %)

Usage:
  python validate_est.py \
      --csv       genes.csv \
      --est       all_est.txt \
      --refgene   refGene.txt \
      --out       results.csv \
      [--identity 96.0]

Dependencies:
  pip install biopython pandas requests
"""

import argparse
import sys
import time
from typing import Optional

import pandas as pd
import requests
from Bio import pairwise2

# Persistent HTTP session with sensible headers
session = requests.Session()
session.headers.update({"User-Agent": "validate_est/1.0 (research; contact@example.com)"})


# ── helpers ──────────────────────────────────────────────────────────────────

def parse_psl_file(path: str, table_type: str) -> pd.DataFrame:
    """
    Parse a PSL-format txt file (all_EST or refGene).
    Skips comment / header lines (those starting with # or 'bin').
    """
    rows = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("bin"):
                continue
            parts = line.split("\t")
            rows.append(parts)

    if table_type == "est":
        # PSL columns from schema (0-based)
        cols = [
            "bin", "matches", "misMatches", "repMatches", "nCount",
            "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert",
            "strand", "qName", "qSize", "qStart", "qEnd",
            "tName", "tSize", "tStart", "tEnd",
            "blockCount", "blockSizes", "qStarts", "tStarts",
        ]
        df = pd.DataFrame(rows, columns=cols[:len(rows[0])] if rows else cols)
        for c in ["matches", "misMatches", "repMatches", "nCount",
                  "qStart", "qEnd", "tStart", "tEnd", "qSize"]:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    else:  # refGene
        cols = [
            "bin", "name", "chrom", "strand",
            "txStart", "txEnd", "cdsStart", "cdsEnd",
            "exonCount", "exonStarts", "exonEnds",
            "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames",
        ]
        df = pd.DataFrame(rows, columns=cols[:len(rows[0])] if rows else cols)
        for c in ["txStart", "txEnd", "cdsStart", "cdsEnd"]:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


def pct_identity(matches: int, misMatches: int, repMatches: int) -> float:
    """Classic BLAT % identity formula."""
    total = matches + misMatches + repMatches
    return (matches + repMatches) / total * 100 if total > 0 else 0.0


def seq_pct_identity(seq_a: str, seq_b: str) -> float:
    """
    Pairwise % identity between two sequences using a simple
    global alignment (Biopython pairwise2).
    Returns identity % over the alignment length.
    """
    if not seq_a or not seq_b:
        return 0.0
    alignments = pairwise2.align.globalms(
        seq_a.upper(), seq_b.upper(),
        2, -1, -2, -0.5,          # match, mismatch, gap-open, gap-extend
        one_alignment_only=True,
    )
    if not alignments:
        return 0.0
    aln = alignments[0]
    aln_len = len(aln.seqA)
    matches = sum(a == b for a, b in zip(aln.seqA, aln.seqB) if a != "-" and b != "-")
    return matches / aln_len * 100 if aln_len else 0.0


def fetch_refseq_nucleotide(refseq_id: str) -> str:
    """
    Fetch RefSeq nucleotide sequence from NCBI eutils, skipping complete genomes.

    Args:
        refseq_id: RefSeq accession (e.g. NM_001234).
    Returns:
        Nucleotide sequence string, or "" if unavailable / skipped.
    """
    summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    fetch_url   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    # 1. Check summary — skip whole-genome / complete-sequence records
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

    # 2. Fetch the FASTA sequence
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


def find_overlapping_ests(est_df: pd.DataFrame, chrom: str,
                           start: int, end: int,
                           min_identity: float = 96.0) -> pd.DataFrame:
    """
    Return EST rows that overlap [start, end) on the given chromosome
    and meet the minimum % identity threshold.
    """
    mask = (
        (est_df["tName"] == chrom) &
        (est_df["tStart"] < end) &
        (est_df["tEnd"] > start)
    )
    hits = est_df[mask].copy()
    if hits.empty:
        return hits
    hits["pct_identity"] = hits.apply(
        lambda r: pct_identity(r["matches"], r["misMatches"], r["repMatches"]),
        axis=1,
    )
    return hits[hits["pct_identity"] >= min_identity]


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--csv",      required=True, help="Gene CSV file")
    parser.add_argument("--est",      required=True, help="all_EST PSL txt file")
    parser.add_argument("--refgene",  required=True, help="refGene PSL txt file")
    parser.add_argument("--out",      default="est_validation_results.csv",
                        help="Output CSV path (default: est_validation_results.csv)")
    parser.add_argument("--identity", type=float, default=96.0,
                        help="Minimum %% identity threshold (default: 96.0)")
    args = parser.parse_args()

    # ── 1. Load gene CSV ──────────────────────────────────────────────────
    print("Loading gene CSV …")
    genes_df = pd.read_csv(args.csv)
    required_cols = {"ensembl_gene_id", "cds", "utr5", "utr3"}
    missing = required_cols - set(genes_df.columns)
    if missing:
        sys.exit(f"ERROR: gene CSV is missing columns: {missing}")

    # Use bovine_id as the refGene lookup key if present, else ensembl_gene_id
    id_col = "bovine_id" if "bovine_id" in genes_df.columns else "ensembl_gene_id"
    print(f"  {len(genes_df)} genes loaded. Using '{id_col}' for refGene lookup.")

    # ── 2. Load refGene ───────────────────────────────────────────────────
    print("Loading refGene table …")
    refgene_df = parse_psl_file(args.refgene, "refgene")
    print(f"  {len(refgene_df)} refGene entries loaded.")

    # ── 3. Load all_EST ───────────────────────────────────────────────────
    print("Loading all_EST table …")
    est_df = parse_psl_file(args.est, "est")
    print(f"  {len(est_df)} EST alignments loaded.")

    # ── 4. Process each gene ──────────────────────────────────────────────
    results = []
    ncbi_cache: dict[str, Optional[str]] = {}   # accession → sequence

    for idx, row in genes_df.iterrows():
        gene_id   = row["ensembl_gene_id"]
        accession = str(row.get(id_col, "")).strip()
        cds_seq   = str(row.get("cds",  "") or "").strip()
        utr5_seq  = str(row.get("utr5", "") or "").strip()
        utr3_seq  = str(row.get("utr3", "") or "").strip()

        print(f"\n[{idx+1}/{len(genes_df)}] {gene_id}  ({accession})")

        result = {
            "ensembl_gene_id":    gene_id,
            "accession":          accession,
            "chrom":              "",
            "strand":             "",
            "tx_start":           "",
            "tx_end":             "",
            # EST overlap info
            "est_overlap_count":  0,
            "est_best_identity":  "",
            "est_names":          "",
            # Sequence identity vs NCBI mRNA  (utr5 + cds + utr3 concatenated)
            "ncbi_seq_fetched":   False,
            "query_mrna_length":  "",
            "ncbi_mrna_length":   "",
            "mrna_identity":      "",
            "mrna_pass":          "",
            "overall_pass":       "",
            "notes":              "",
        }

        # ── 4a. Look up refGene ───────────────────────────────────────────
        rg_hits = refgene_df[refgene_df["name"] == accession]
        if rg_hits.empty:
            # Try matching on name2 (gene symbol) as fallback
            rg_hits = refgene_df[refgene_df["name2"] == accession]
        if rg_hits.empty:
            print(f"  [WARN] No refGene entry found for '{accession}' — skipping coordinate lookup.")
            result["notes"] = "No refGene entry"
            results.append(result)
            continue

        rg = rg_hits.iloc[0]
        chrom    = rg["chrom"]
        strand   = rg["strand"]
        tx_start = int(rg["txStart"])
        tx_end   = int(rg["txEnd"])
        result.update({"chrom": chrom, "strand": strand,
                        "tx_start": tx_start, "tx_end": tx_end})
        print(f"  refGene: {chrom}:{tx_start}-{tx_end} ({strand})")

        # ── 4b. Find overlapping ESTs ─────────────────────────────────────
        overlapping = find_overlapping_ests(
            est_df, chrom, tx_start, tx_end, args.identity
        )
        result["est_overlap_count"] = len(overlapping)
        if not overlapping.empty:
            result["est_best_identity"] = round(overlapping["pct_identity"].max(), 2)
            result["est_names"] = ";".join(overlapping["qName"].tolist()[:10])  # cap at 10
            print(f"  ESTs overlapping (≥{args.identity}%): {len(overlapping)}  "
                  f"best={result['est_best_identity']:.2f}%")
        else:
            print(f"  No ESTs overlapping at ≥{args.identity}% identity.")
            result["notes"] = "No supporting ESTs"

        # ── 4c. Fetch NCBI mRNA sequence ──────────────────────────────────
        if accession not in ncbi_cache:
            print(f"  Fetching NCBI sequence for {accession} …")
            ncbi_cache[accession] = fetch_refseq_nucleotide(accession)
            time.sleep(0.4)   # stay within NCBI rate limits (~3 req/s)
        ncbi_seq = ncbi_cache[accession]

        if not ncbi_seq:
            result["notes"] = (result["notes"] + "; NCBI fetch failed").lstrip("; ")
            results.append(result)
            continue

        result["ncbi_seq_fetched"] = True
        result["ncbi_mrna_length"] = len(ncbi_seq)
        print(f"  NCBI sequence length: {len(ncbi_seq)} bp")

        # ── 4d. Align utr5 + cds + utr3 (full mRNA) against NCBI seq ──────
        # Concatenate in transcript order: 5'UTR → CDS → 3'UTR
        query_mrna = (utr5_seq + cds_seq + utr3_seq).upper()
        result["query_mrna_length"] = len(query_mrna)

        if not query_mrna:
            result["notes"] = (result["notes"] + "; empty query mRNA").lstrip("; ")
            results.append(result)
            continue

        print(f"  Query mRNA (5'UTR+CDS+3'UTR): {len(query_mrna)} bp")
        mrna_id = seq_pct_identity(query_mrna, ncbi_seq)
        mrna_pass = mrna_id >= args.identity
        print(f"  mRNA identity: {mrna_id:.2f}%  {'PASS' if mrna_pass else 'FAIL'}")

        result["mrna_identity"] = round(mrna_id, 2)
        result["mrna_pass"]     = mrna_pass
        result["overall_pass"]  = mrna_pass and result["est_overlap_count"] > 0

        results.append(result)

    # ── 5. Write output CSV ───────────────────────────────────────────────
    out_df = pd.DataFrame(results)
    out_df.to_csv(args.out, index=False)
    print(f"\n✓ Results written to: {args.out}")

    # ── 6. Summary ────────────────────────────────────────────────────────
    total   = len(out_df)
    fetched = out_df["ncbi_seq_fetched"].sum()
    if fetched > 0:
        passed = out_df[out_df["ncbi_seq_fetched"]]["mrna_pass"].sum()
        print(f"\nSummary")
        print(f"  Total genes:                  {total}")
        print(f"  NCBI seq fetched:             {fetched}")
        print(f"  mRNA identity PASS (≥{args.identity}%): {passed}/{fetched}")
        with_ests = out_df["est_overlap_count"].gt(0).sum()
        print(f"  Genes with EST support:       {with_ests}/{total}")
    else:
        print(f"\nNo sequences fetched from NCBI — check accession IDs and network access.")


if __name__ == "__main__":
    main()