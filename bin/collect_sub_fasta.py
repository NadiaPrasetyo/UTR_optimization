import argparse
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

def main(input_file, output_file, accessions):
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        if record.id.split(" ")[0] in accessions:
            records.append(record)

    SeqIO.write(records, output_file, "fasta")

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, help="Input FASTA file")
    parser.add_argument("-o", "--output_file", type=str, help="Output file")
    parser.add_argument("--accessions", nargs="+", type=str, default=["NC_024120.1", "KC663628.1","PP346853.1","NC_033793.1","KY369300.1","KY369299.1","PQ368551.1","PP228890.1","PP228891.1","PP228897.1","KF961187.1","NC_024768.1","KF979335.1","NC_024769.1","KF979336.1","NC_023858.1","KF961188.1","OP169446.1","OR364988.1","NC_038957.1","KC811837.1","PQ368550.1","NC_039235.1","KF961186.1"], help="Accessions to include")
    args = parser.parse_args()

    if not os.path.exists(os.path.dirname(args.output_file)):
        os.makedirs(os.path.dirname(args.output_file))

    main(args.input_file, args.output_file, args.accessions)