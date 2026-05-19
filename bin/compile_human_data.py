import argparse
import pandas as pd
import os
from Bio import SeqIO


def extract_length_gc(sequence):
    """Compute length and GC content of a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.

    Returns:
        tuple[int, float]: (length, gc_content_percent).
    """
    if not sequence:
        return 0, 0.0
    seq = str(sequence).upper()
    length = len(seq)
    if length == 0:
        return 0, 0.0
    gc = seq.count("G") + seq.count("C")
    return length, round(gc / length * 100, 2)


def read_fasta_dir(directory, suffix, key_index=0):
    """Read all .fa files in a directory and return a dict of {accession: sequence}.

    Args:
        directory (str): Path to directory containing .fa files.
        suffix (str): Column name suffix for logging.
        key_index (int): Index of the accession field in the filename (split by '_').

    Returns:
        dict[str, str]: Mapping of accession -> sequence string.
    """
    records = {}
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Directory not found: {directory}")
    for fname in os.listdir(directory):
        if not fname.endswith(".fa"):
            continue
        parts = fname.split("_")
        accession = parts[key_index]
        fpath = os.path.join(directory, fname)
        record = SeqIO.read(fpath, "fasta")
        records[accession] = str(record.seq)
    print(f"  Loaded {len(records)} {suffix} sequences from {directory}")
    return records


def extract_CDS_UTRs(input_dir):
    """Parse CDS, 5'UTR, and 3'UTR FASTA files from subdirectories.

    Expected structure:
        input_dir/
            3UTR/  <accession>_<ensemblID>_3UTR.fa
            5UTR/  <accession>_<ensemblID>_5UTR.fa
            CDS/   <accession>_<ensemblID>_CDS.fa

    Args:
        input_dir (str): Root directory containing CDS, 3UTR, and 5UTR subdirs.

    Returns:
        pd.DataFrame: One row per accession with columns:
            accession, ensembl_gene_id, ensembl_transcript_id, utr3, utr5, cds.
    """
    utr3_dir = os.path.join(input_dir, "3UTR")
    utr5_dir = os.path.join(input_dir, "5UTR")
    cds_dir  = os.path.join(input_dir, "CDS")

    # Build base table from 3UTR files (carries ensembl_gene_id)
    rows = []
    for fname in sorted(os.listdir(utr3_dir)):
        if not fname.endswith(".fa"):
            continue
        parts = fname.split("_")
        accession  = parts[0]
        fpath = os.path.join(utr3_dir, fname)
        ensembl_gene_id = SeqIO.read(fpath, "fasta").id.split("|")[1]
        ensembl_transcript_id = SeqIO.read(fpath, "fasta").id.split("|")[0]
        seq = str(SeqIO.read(fpath, "fasta").seq)
        rows.append({"accession": accession, "ensembl_gene_id": ensembl_gene_id, "ensembl_transcript_id": ensembl_transcript_id, "utr3": seq})

    df = pd.DataFrame(rows)
    print(f"  Loaded {len(df)} 3UTR sequences from {utr3_dir}")

    # Merge 5UTR sequences
    utr5_map = read_fasta_dir(utr5_dir, suffix="5UTR", key_index=0)
    df["utr5"] = df["accession"].map(utr5_map)

    # Merge CDS sequences
    cds_map = read_fasta_dir(cds_dir, suffix="CDS", key_index=0)
    df["cds"] = df["accession"].map(cds_map)

    missing_utr5 = df["utr5"].isna().sum()
    missing_cds  = df["cds"].isna().sum()
    if missing_utr5:
        print(f"  WARNING: {missing_utr5} accessions have no matching 5UTR sequence.")
    if missing_cds:
        print(f"  WARNING: {missing_cds} accessions have no matching CDS sequence.")

    print(f"Total accessions extracted: {len(df)}")
    return df


def strip_version(ensembl_id):
    """Strip the version suffix from an Ensembl ID (e.g. ENSG00000121410.14 -> ENSG00000121410).

    Args:
        ensembl_id (str): Ensembl ID with or without a version suffix.

    Returns:
        str: Ensembl ID without version suffix.
    """
    if isinstance(ensembl_id, str) and "." in ensembl_id:
        return ensembl_id.split(".")[0]
    return ensembl_id


def compile_te_data(te_csv):
    """Load translation efficiency data and return a DataFrame keyed on
    (gene_id_base, transcript_id_base), with Ensembl version suffixes stripped.

    Expected CSV fields: gene_id, transcript_id, mean_te

    Args:
        te_csv (str): Path to the TE CSV file.

    Returns:
        pd.DataFrame: Columns [gene_id_base, transcript_id_base, mean_te].
    """
    te = pd.read_csv(te_csv)
    required = {"gene_id", "transcript_id", "mean_te"}
    missing = required - set(te.columns)
    if missing:
        raise ValueError(f"TE file is missing columns: {missing}")

    te["gene_id_base"]       = te["gene_id"].map(strip_version)
    te["transcript_id_base"] = te["transcript_id"].map(strip_version)
    print(f"  Loaded {len(te)} TE entries from {te_csv}")
    return te[["gene_id_base", "transcript_id_base", "mean_te"]]


def compile_hl_data(hl_csv):
    """Load mRNA half-life data and return a DataFrame keyed on gene_id_base,
    with Ensembl version suffixes stripped.

    Expected CSV fields: gene_id, half_life

    Args:
        hl_csv (str): Path to the half-life CSV file.

    Returns:
        pd.DataFrame: Columns [gene_id_base, half_life].
    """
    hl = pd.read_csv(hl_csv)
    required = {"gene_id", "half_life"}
    missing = required - set(hl.columns)
    if missing:
        raise ValueError(f"Half-life file is missing columns: {missing}")

    hl["gene_id_base"] = hl["gene_id"].map(strip_version)
    print(f"  Loaded {len(hl)} half-life entries from {hl_csv}")
    return hl[["gene_id_base", "half_life"]]


def main(input_dir, output_file, te_csv=None, hl_csv=None):
    """Extract sequences, compute length/GC metrics, optionally merge TE and HL
    data, then write to CSV.

    Args:
        input_dir (str): Root directory with CDS/, 5UTR/, 3UTR/ subdirectories.
        output_file (str): Path for the output CSV.
        te_csv (str | None): Path to the TE CSV file (optional).
        hl_csv (str | None): Path to the half-life CSV file (optional).
    """
    df = extract_CDS_UTRs(input_dir)

    for region in ("utr5", "utr3", "cds"):
        df[f"{region}_length"], df[f"{region}_gc"] = zip(
            *df[region].map(extract_length_gc)
        )

    # Strip versions from the IDs extracted from FASTA headers so we can join
    # against the TE/HL tables regardless of version mismatches.
    df["gene_id_base"]       = df["ensembl_gene_id"].map(strip_version)
    df["transcript_id_base"] = df["ensembl_transcript_id"].map(strip_version)

    if te_csv:
        te_df = compile_te_data(te_csv)
        before = len(df)
        df = df.merge(te_df, on=["gene_id_base", "transcript_id_base"], how="left")
        matched = df["mean_te"].notna().sum()
        print(f"  TE merge: {matched}/{before} rows matched.")

    if hl_csv:
        hl_df = compile_hl_data(hl_csv)
        before = len(df)
        df = df.merge(hl_df, on="gene_id_base", how="left")
        matched = df["half_life"].notna().sum()
        print(f"  Half-life merge: {matched}/{before} rows matched.")

    # Drop the temporary base-ID columns used for merging
    df.drop(columns=["gene_id_base", "transcript_id_base"], inplace=True)

    df.to_csv(output_file, index=False)
    print(f"Output written to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Extract length and GC content for CDS, 5UTR, and 3UTR sequences, "
            "and optionally merge translation efficiency and mRNA half-life data."
        ),
    )
    parser.add_argument(
        "-i", "--input-dir",
        required=True,
        type=str,
        help="Input directory containing CDS/, 5UTR/, and 3UTR/ subdirectories of .fa files.",
    )
    parser.add_argument(
        "-o", "--output-file",
        default="data/seq_length_gc.csv",
        type=str,
        help="Output CSV file path (default: data/seq_length_gc.csv).",
    )
    parser.add_argument(
        "--te-csv",
        default=None,
        type=str,
        help="CSV with translation efficiency data (fields: gene_id, transcript_id, mean_te).",
    )
    parser.add_argument(
        "--hl-csv",
        default=None,
        type=str,
        help="CSV with mRNA half-life data (fields: gene_id, half_life).",
    )
    args = parser.parse_args()

    out_dir = os.path.dirname(args.output_file)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    main(args.input_dir, args.output_file, te_csv=args.te_csv, hl_csv=args.hl_csv)