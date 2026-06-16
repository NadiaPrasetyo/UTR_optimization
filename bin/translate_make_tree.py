import argparse
import os
import subprocess
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from Bio import Phylo
import sys

def visualize_tree(tree_file, highlight_accessions=None):
    # Increase recursion limit for large trees
    sys.setrecursionlimit(10000)

    tree = Phylo.read(tree_file, "newick")
    fig, ax = plt.subplots(figsize=(12, 20))  # taller figure for large trees

    # Color highlighted clades
    if highlight_accessions:
        for clade in tree.find_clades():
            if clade.name and any(acc in clade.name for acc in highlight_accessions):
                clade.color = "red"

    Phylo.draw(tree, axes=ax, do_show=False)
    plt.tight_layout()
    plt.savefig(tree_file.replace(".tree", ".pdf"), bbox_inches="tight", height=250)  # save to file too


def translate_sequences(input_file):
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        seq = record.seq.translate(to_stop=True)
        record.id = record.id.split(" ")[0]
        record.description = ""
        record.seq = seq
        records.append(record)
    return records

def main():
    parser = argparse.ArgumentParser(description="Translate, align, and visualize sequences")
    parser.add_argument("-i", "--input",    required=True,  help="Input FASTA file with full genome sequences to be translated")
    parser.add_argument("-o", "--output",   required=True,  help="Output alignment and tree file prefix")
    parser.add_argument("--highlight",      nargs="+",      help="Accession IDs to highlight in the tree", default=[])
    args = parser.parse_args()

    # Ensure output directory exists
    out_dir = os.path.dirname(args.output)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Translate sequences
    records = translate_sequences(args.input)

    # Write temporary FASTA file
    temp_fa = "temp.fa"
    with open(temp_fa, "w") as fh:
        for record in records:
            fh.write(f">{record.id}\n{record.seq}\n")

    # Run MAFFT to align sequences and output tree
    alignment_file = f"{args.output}.aln"
    tree_file      = f"{args.output}.tree"
    try:
        with open(alignment_file, "w") as aln_out:
            subprocess.run(
                ["mafft", "--auto", "--treeout", temp_fa],
                stdout=aln_out,
                check=True
            )
        # MAFFT writes the tree as <input>.tree, so rename it
        mafft_tree = f"{temp_fa}.tree"
        if os.path.exists(mafft_tree):
            shutil.move(mafft_tree, tree_file)

        print(f"Alignment written to {alignment_file}")
        print(f"Tree written to {tree_file}")

    finally:
        # Clean up temp file even if something fails
        if os.path.exists(temp_fa):
            os.remove(temp_fa)

    # Visualize tree with optional highlighting
    if os.path.exists(tree_file):
        visualize_tree(tree_file, highlight_accessions=args.highlight)
    else:
        print(f"Warning: tree file not found at {tree_file}")

if __name__ == "__main__":
    main()