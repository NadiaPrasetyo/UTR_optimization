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
    """
    Visualize a Newick phylogenetic tree and save it as an SVG.
    Args:
        tree_file: Path to the .tree (Newick format) file.
        highlight_accessions: Optional list of accession strings to highlight
                              with a turquoise label background.
    """
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_svg import FigureCanvasSVG

    sys.setrecursionlimit(10000)
    print(f"Visualizing tree: {tree_file}")

    tree = Phylo.read(tree_file, "newick")

    # Rename tips only — internal nodes have name=None, skip those
    # Also guard against names that don't have enough underscore-delimited parts
    for clade in tree.find_clades():
        if clade.name:
            parts = clade.name.split("_")
            if len(parts) >= 3:
                if parts[1] == "NC":
                    clade.name = "NC_" + parts[2] + "." + parts[3]
                else:
                    clade.name = parts[1] + "." + parts[2]

    # Build a set of names to highlight for fast lookup
    highlight_set = set()
    if highlight_accessions:
        for clade in tree.find_clades():
            if clade.name and any(acc in clade.name for acc in highlight_accessions):
                clade.color = "red"
                highlight_set.add(clade.name)

    num_tips = tree.count_terminals()
    fig_height = max(20, num_tips * 0.3)
    print(f"Tree has {num_tips} tips — setting figure height to {fig_height:.0f} inches")

    fig, ax = plt.subplots(figsize=(16, fig_height))
    Phylo.draw(tree, axes=ax, do_show=False)

    # Add turquoise background boxes behind highlighted tip labels
    if highlight_set:
        for text in ax.texts:
            if text.get_text().strip() in highlight_set:
                text.set_backgroundcolor("#00CED1")
                text.get_bbox_patch().set_boxstyle("round,pad=0.15")
                text.get_bbox_patch().set_edgecolor("none")
                text.set_color("black")

    plt.tight_layout()

    output_path = os.path.splitext(tree_file)[0] + ".svg"
    canvas = FigureCanvasSVG(fig)
    with open(output_path, "w") as f:
        canvas.print_svg(f)
    plt.close(fig)

    print(f"Tree saved to: {output_path}")

def deduplicate_aa_seq(records):
    seen = set()
    unique_records = []
    for record in records:
        if record.id not in seen:
            unique_records.append(record)
            seen.add(record.id)
    return unique_records

def translate_sequences(input_file):
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        seq = record.seq.translate(to_stop=False)
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
    parser.add_argument("--no-translate",   action="store_true", help="Skip translation step")
    parser.add_argument("--tree",           help="Visualize existing tree")
    args = parser.parse_args()

    # Ensure output directory exists
    out_dir = os.path.dirname(args.output)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    if args.tree:
        visualize_tree(args.tree, highlight_accessions=args.highlight)
        return
    # Translate sequences
    if not args.no_translate:
        records = translate_sequences(args.input)
    else:
        records = list(SeqIO.parse(args.input, "fasta"))

    # Deduplicate sequences
    records = deduplicate_aa_seq(records)

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
