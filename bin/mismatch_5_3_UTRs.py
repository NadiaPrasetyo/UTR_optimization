import argparse
import os
import subprocess
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-dir", required=True, help="Input directory that contains 3' and 5' UTR sequences")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    records = []
    for record in SeqIO.parse(args.input, "fasta"):
        seq = record.seq