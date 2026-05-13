import argparse
import pandas as pd
import csv
from collections import defaultdict
import os
from Bio import SeqIO


def extract_length_gc(sequence):
    """Compute length and GC content of a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence string.

    Returns:
        tuple[int, float]: (length, gc_content_percent).
    """
    if not isinstance(sequence, str):
        return 0, 0.0
    length = len(sequence)
    if length == 0:
        return 0, 0.0
    gc = sequence.upper().count("G") + sequence.upper().count("C")
    return length, round(gc / length * 100, 2)


def extract_CDS_UTRs(input_dir):
    utr3_file = os.path.join(input_dir, 'extracted_3UTR.fa')
    utr5_file = os.path.join(input_dir, 'extracted_5UTR.fa')
    cds_file = os.path.join(input_dir, 'extracted_CDS.fa')

    utr3_dict = {record.id.split('_')[0]: str(record.seq) for record in SeqIO.parse(utr3_file, "fasta")}
    utr5_dict = {record.id.split('_')[0]: str(record.seq) for record in SeqIO.parse(utr5_file, "fasta")}
    cds_dict = {record.id.split('_')[0]: str(record.seq) for record in SeqIO.parse(cds_file, "fasta")}

    df = pd.DataFrame(list(cds_dict.items()), columns=['ensembl_gene_id', 'cds'])
    df['utr5'] = df['ensembl_gene_id'].map(utr5_dict)
    df['utr3'] = df['ensembl_gene_id'].map(utr3_dict)

    print (f"Number of extracted sequences: {len(df)}")
    return df


def main(input_dir, tsv_file, output_file):
    # combine the input file, half life and transcript efficiency, and get the length and GC content for cds, utr5, and utr3
    df = extract_CDS_UTRs(input_dir)

    df['utr5_length'], df['utr5_gc'] = zip(*df['utr5'].map(extract_length_gc))
    df['utr3_length'], df['utr3_gc'] = zip(*df['utr3'].map(extract_length_gc))

    if tsv_file:
        df = pd.merge(df, pd.read_csv(tsv_file, sep='\t'), how='left', on='ensembl_gene_id')


    df.to_csv(output_file, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract length and GC content from a CSV summary file')
    parser.add_argument('-i', '--input-dir', required=True, type=str, help='Input Fasta Directory containing CDS, UTR5 and UTR3 sequences')
    parser.add_argument('-f', '--tsv-file', required=True, type=str, help='Input tsv file containing ensembl ids, half life, and transcript efficiency')
    parser.add_argument('-o', '--output-file', default='data/seq_length_gc.csv', type=str, help='Output CSV file (default: length_gc.csv)')

    args = parser.parse_args()

    # check that the output directory exists
    if not os.path.exists(os.path.dirname(args.output_file)):
        os.makedirs(os.path.dirname(args.output_file))

    main(args.input_dir, args.tsv_file, args.output_file)