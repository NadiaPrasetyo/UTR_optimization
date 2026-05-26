import csv
import pandas as pd
import numpy as np
import os
import sys
import argparse

def main(input_file, output_file):
    df = pd.read_csv(input_file)
    df = df.replace(np.nan, '', regex=True)

    # get the accession, accession,ensembl_gene_id,ensembl_transcript_id, sequence = utr5
    df = df[['accession', 'ensembl_gene_id', 'ensembl_transcript_id', 'utr5', 'half_life', 'mean_te']]

    # add a human baseline:
    hbb_5utr_sequence = 'ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACC'
    df = df.append({
        'accession': 'HBB-201',
        'ensembl_gene_id': 'ENSG00000139618',
    df.to_csv(output_file, sep='\t', index=False)


    print("Done ! Now you can run:\n\tRiboScanner predict\n\t--input ./example_data/input.txt\n\t--column_sequence utr5\n\t--output ./output.txt")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, help='Input CSV file')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output TSV file')
    args = parser.parse_args()

    if not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    main(args.input, args.output)
