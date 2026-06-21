import argparse
import os
import re
import subprocess
from Bio import SeqIO


def build_mafft(input_file, output_file, verbose=False):
    """Align sequences using MAFFT and return the path to the alignment file."""
    aln_file = f"{output_file}.aln"
    print(f"Running MAFFT on {input_file}...")
    with open(aln_file, "w") as fh:
        subprocess.run(
            ["mafft", "--auto", input_file],
            stdout=fh,
            check=True,
            stderr=None if verbose else subprocess.DEVNULL,
        )
    print(f"  Alignment saved to: {aln_file}")
    return aln_file


def run_rnaalifold(input_file, output_file, verbose=False):
    """
    Run RNAalifold on a FASTA/CLUSTAL MSA.

    Returns the consensus dot-bracket structure string (without the energy term)
    so it can be embedded into a Stockholm file for R-scape.
    """
    print(f"Running RNAalifold on {input_file}...")
    result = subprocess.run(
        ["RNAalifold", "--aln", input_file],
        capture_output=True,
        text=True,
        check=True,
    )
    if verbose and result.stderr:
        print(result.stderr)

    fold_file = f"{output_file}.rnaalifold"
    with open(fold_file, "w") as fh:
        fh.write(result.stdout)
    print(f"  RNAalifold output saved to: {fold_file}")

    # Extract the dot-bracket structure from the last output line.
    # RNAalifold stdout format (last line):
    #   <dot-bracket>  (<energy> = <energy>)
    last_line = result.stdout.strip().splitlines()[-1]
    match = re.match(r"^([.()[\]{}<>]+)", last_line)
    if not match:
        raise ValueError(
            f"Could not parse dot-bracket structure from RNAalifold output:\n{last_line}"
        )
    structure = match.group(1)
    print(f"  Consensus structure: {structure}")
    return structure


def fasta_aln_to_stockholm(aln_file, sto_file, ss_cons=None, verbose=False):
    """
    Convert a FASTA-format MSA to Stockholm format and optionally embed a
    consensus secondary structure annotation (#=GC SS_cons).

    Writes a single unwrapped block (all sequences on one line each) so that
    the #=GC SS_cons annotation sits unambiguously in the same block as the
    sequence rows -- multi-block wrapped output from esl-reformat causes R-scape
    to reject the SS_cons line with "unexpected #=GC SS_cons; earlier block(s)
    in different order".

    Parameters
    ----------
    aln_file : str
        Path to the FASTA alignment produced by MAFFT.
    sto_file : str
        Destination path for the Stockholm file.
    ss_cons  : str or None
        Dot-bracket consensus structure from RNAalifold.  When provided its
        length is validated against the alignment length.
        Pass None to produce a bare Stockholm file (e.g. for --cacofold).
    """
    records = list(SeqIO.parse(aln_file, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {aln_file}")

    aln_len = len(records[0].seq)

    if ss_cons is not None and len(ss_cons) != aln_len:
        raise ValueError(
            f"Structure length ({len(ss_cons)}) does not match "
            f"alignment length ({aln_len}). "
            "RNAalifold may have been run on a different alignment."
        )

    # Pad names so sequence columns are aligned (purely cosmetic)
    col_width = max(len(r.id) for r in records)
    if ss_cons is not None:
        col_width = max(col_width, len("#=GC SS_cons"))

    print(f"  Writing Stockholm MSA to {sto_file}...")
    with open(sto_file, "w") as fh:
        fh.write("# STOCKHOLM 1.0\n\n")
        for r in records:
            fh.write(f"{r.id:<{col_width}}  {r.seq}\n")
        if ss_cons is not None:
            fh.write(f"{'#=GC SS_cons':<{col_width}}  {ss_cons}\n")
        fh.write("//\n")

    print(f"  Stockholm MSA saved to: {sto_file}")
    return sto_file


def run_rscape(input_file, output_file, use_cacofold=False, verbose=False):
    """
    Run R-scape covariation analysis.

    Parameters
    ----------
    input_file   : str
        Stockholm MSA.  When use_cacofold=False this must already contain a
        #=GC SS_cons line; when use_cacofold=True R-scape predicts its own
        CaCoFold structure and writes a new *.cacofold.sto file.
    use_cacofold : bool
        Pass --cacofold so R-scape predicts structure internally rather than
        reading it from the #=GC SS_cons annotation.
    """
    print(f"Running R-scape on {input_file}...")
    out_dir = os.path.dirname(output_file) or "."
    cmd = [
        "R-scape",
        "--outdir",  out_dir,
        "--outname", os.path.basename(output_file),
    ]
    if use_cacofold:
        cmd.append("--cacofold")
    else:
        # Analyse the structure already in the Stockholm file
        cmd.append("-s")
    cmd.append(input_file)

    os.makedirs(out_dir, exist_ok=True)
    subprocess.run(
        cmd,
        check=True,
        stdout=None if verbose else subprocess.DEVNULL,
        stderr=None,   # always show R-scape stderr so errors are visible
    )

    if use_cacofold:
        cacofold_sto = f"{output_file}.cacofold.sto"
        print(f"  R-scape (CaCoFold) results prefix : {output_file}")
        print(f"  Annotated Stockholm file           : {cacofold_sto}")
    else:
        print(f"  R-scape results prefix: {output_file}")


def run_nhmmer(query, subject, output_file, verbose=False):
    """Build a profile HMM from query sequences and search a target database."""
    print(f"Running nhmmer — query: {query}  subject: {subject}")
    tmp_hmm_profile = query.replace(".fasta", ".hmm")

    subprocess.run(
        ["hmmbuild", tmp_hmm_profile, query],
        check=True,
        stdout=None if verbose else subprocess.DEVNULL,
        stderr=None if verbose else subprocess.DEVNULL,
    )
    subprocess.run(
        [
            "nhmmer",
            "-o",       f"{output_file}.txt",
            "--tblout", f"{output_file}.tblout",
            "-A",       f"{output_file}.sto",
            tmp_hmm_profile,
            subject,
        ],
        check=True,
        stdout=None if verbose else subprocess.DEVNULL,
        stderr=None if verbose else subprocess.DEVNULL,
    )

    if not verbose:
        print("  Removing temporary HMM profile...")
        os.remove(tmp_hmm_profile)

    print(f"  nhmmer results saved to: {output_file}.txt / .tblout / .sto")


def main(input_file, output_file, use_cacofold=False, verbose=False):
    # ------------------------------------------------------------------ #
    # 1. Split query (first sequence) from the rest                       #
    # ------------------------------------------------------------------ #
    records = list(SeqIO.parse(input_file, "fasta"))
    if len(records) < 2:
        raise ValueError("Input FASTA must contain at least two sequences.")

    print(f"Loaded {len(records)} sequences from {input_file}")
    print(f"  Query   : {records[0].id}")
    print(f"  Subjects: {len(records) - 1} remaining sequences")

    query_file = input_file.replace(".fasta", "_query.fasta")
    with open(query_file, "w") as fh:
        fh.write(f">{records[0].id}\n{records[0].seq}\n")

    # ------------------------------------------------------------------ #
    # 2. nhmmer: profile search                                           #
    # ------------------------------------------------------------------ #
    run_nhmmer(query_file, input_file, output_file, verbose)

    if not verbose:
        print("  Removing temporary query file...")
        os.remove(query_file)

    # ------------------------------------------------------------------ #
    # 3. MAFFT: multiple sequence alignment                               #
    # ------------------------------------------------------------------ #
    aln_file = build_mafft(input_file, output_file, verbose)

    # ------------------------------------------------------------------ #
    # 4. RNAalifold + Stockholm conversion  OR  R-scape CaCoFold         #
    # ------------------------------------------------------------------ #
    sto_file = f"{output_file}.mafft.sto"

    if use_cacofold:
        # Bare Stockholm -- R-scape predicts structure internally via CaCoFold
        print("CaCoFold mode: converting alignment to Stockholm (no SS_cons)...")
        fasta_aln_to_stockholm(aln_file, sto_file, ss_cons=None, verbose=verbose)
    else:
        # RNAalifold -> dot-bracket -> annotated Stockholm for R-scape -s
        ss_cons = run_rnaalifold(aln_file, output_file, verbose)
        fasta_aln_to_stockholm(aln_file, sto_file, ss_cons=ss_cons, verbose=verbose)

    # ------------------------------------------------------------------ #
    # 5. R-scape: covariation analysis                                    #
    # ------------------------------------------------------------------ #
    run_rscape(sto_file, output_file, use_cacofold=use_cacofold, verbose=verbose)

    # ------------------------------------------------------------------ #
    # Summary                                                             #
    # ------------------------------------------------------------------ #
    print("\nPipeline complete! Output files:")
    for ext in (".aln", ".mafft.sto", ".txt", ".tblout"):
        print(f"  {output_file}{ext}")
    if not use_cacofold:
        print(f"  {output_file}.rnaalifold")
    else:
        print(f"  {output_file}.cacofold.sto  (R-scape annotated Stockholm)")
    print(f"  {output_file}.*  (remaining R-scape output files)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="RNA homology search and structural analysis pipeline."
    )
    parser.add_argument("input_file",  help="Input FASTA file with sequences")
    parser.add_argument("output_file", help="Output file prefix (no extension)")
    parser.add_argument(
        "--cacofold",
        action="store_true",
        help=(
            "Use R-scape's built-in CaCoFold for structure prediction instead of "
            "RNAalifold. R-scape will write a *.cacofold.sto file with the "
            "predicted structure annotated in Stockholm format."
        ),
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Keep temporary files and show tool output",
    )
    args = parser.parse_args()
    
    if not os.path.exists(os.path.dirname(args.output_file)):
        os.makedirs(os.path.dirname(args.output_file))
        
    main(args.input_file, args.output_file, args.cacofold, args.verbose)