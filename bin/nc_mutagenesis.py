import argparse
from Bio import SeqIO
import os


# script to mutate each position of a sequence to all the alternative bases while keeping the original
def mutate(seq):
    for i in range(len(seq)):
        original = seq[i]
        if original not in 'ATGC':
            yield seq
        else:
            for base in 'ATGC':
                if base != original:
                    yield seq[:i] + base + seq[i+1:], str(original) + str(i+1) + str(base)

def main(input_file, output_file):
    with open(input_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            original_seq = record.seq
            original_id = record.id  # Store the original ID
            SeqIO.write(record, output_file, 'fasta') # write the original sequence
            for mutated_seq in mutate(original_seq):
                record.id = f'{original_id}_{mutated_seq[1]}'  # Use original_id, not the modified one
                record.description = ''
                record.seq = mutated_seq[0]
                with open(output_file, 'a') as f:
                    SeqIO.write(record, f, 'fasta') # allow writing on top of existing file


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mutate sequence to all possible bases')
    parser.add_argument('-i', '--input_file', required=True, type=str, help='Input FASTA file')
    parser.add_argument('-o', '--output_file', default='data/nc_mutagenesis.fasta', type=str, help='Output FASTA file (default: nc_mutagenesis.fasta)')

    args = parser.parse_args()

    # check that the output directory exists
    if not os.path.exists(os.path.dirname(args.output_file)):
        os.makedirs(os.path.dirname(args.output_file))

    main(args.input_file, args.output_file)
