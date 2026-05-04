import argparse
import pandas as pd
import csv
from collections import defaultdict
import os


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

def extract_hl_te(tsv_file):
    hl_te_dict = defaultdict(list)
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            accession = row['bovine_id']
            hl_human = row['human_half_life']
            hl_mouse = row['mouse_half_life']
            te_human = row['human_te']
            te_mouse = row['mouse_te']
            hl_te_dict[accession].append(hl_human)
            hl_te_dict[accession].append(hl_mouse)
            hl_te_dict[accession].append(te_human)
            hl_te_dict[accession].append(te_mouse)
    return hl_te_dict

def extract_CDS_UTRs(input_file):
    df = pd.read_csv(input_file)
    # trim the columns to only include uniprot_accession, cds, utr5, utr3
    df = df[['uniprot_accession', 'cds', 'utr5', 'utr3']]
    # convert the sequences to be uppercaps
    df['cds'] = df['cds'].str.upper()
    df['utr5'] = df['utr5'].str.upper()
    df['utr3'] = df['utr3'].str.upper()

    # trim rows with no utr5 and utr3
    df = df[(df['utr5'] != '') & (df['utr3'] != '')]

    return df


def main(input_file, hl_te, output_file):
    # combine the input file, half life and transcript efficiency, and get the length and GC content for cds, utr5, and utr3
    df = extract_CDS_UTRs(input_file)

    if hl_te:
        hl_te_dict = extract_hl_te(hl_te)
        df['hl_human'] = df['uniprot_accession'].map(hl_te_dict).str[0]
        df['hl_mouse'] = df['uniprot_accession'].map(hl_te_dict).str[1]
        df['te_human'] = df['uniprot_accession'].map(hl_te_dict).str[2]
        df['te_mouse'] = df['uniprot_accession'].map(hl_te_dict).str[3]

    df['cds_length'], df['cds_gc'] = zip(*df['cds'].map(extract_length_gc))
    df['utr5_length'], df['utr5_gc'] = zip(*df['utr5'].map(extract_length_gc))
    df['utr3_length'], df['utr3_gc'] = zip(*df['utr3'].map(extract_length_gc))

    df.to_csv(output_file, index=False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract length and GC content from a CSV summary file')
    parser.add_argument('-i', '--input_file', required=True, type=str, help='Input CSV file')
    parser.add_argument('-o', '--output_file', default='data/seq_length_gc.csv', type=str, help='Output CSV file (default: length_gc.csv)')
    parser.add_argument('--hl-te', type=str, help='Half life and Transcript efficiency from TSV file')

    args = parser.parse_args()

    # check that the output directory exists
    if not os.path.exists(os.path.dirname(args.output_file)):
        os.makedirs(os.path.dirname(args.output_file))

    main(args.input_file, args.hl_te, args.output_file)