#!/usr/bin/env python3
"""
compare_af3_candidates.py
==========================

Compare two RNA (or DNA/RNA) structure files in mmCIF format:

  1. Runs `x3dna-dssr` on each .cif to derive secondary structure
     (dot-bracket notation, base pairs, helices, stems, motifs).
  2. Uses Biopython to parse both structures, align matching residues
     (via a pairwise sequence alignment of the DSSR-derived sequences)
     and compute RMSD with Bio.PDB.Superimposer after least-squares
     superposition.
  3. Uses PyMOL (via the `pymol2` API) to independently superimpose the
     two structures (sequence-independent, `super` command), render a
     side-by-side / overlay image, and save a .pse session so the
     result can be inspected interactively.
  4. Prints a side-by-side dot-bracket comparison and a text summary
     report; writes everything to an output directory.

Requirements
------------
  - x3dna-dssr executable on PATH (licensed binary, see https://x3dna.org)
  - PyMOL with the `pymol2` python module (open-source PyMOL via conda:
        conda install -c conda-forge pymol-open-source
    or the commercial PyMOL, which also ships pymol2)
  - Biopython  (pip install biopython)

Usage
-----
    python3 compare_af3_candidates.py structA.cif structB.cif \
        --outdir results/ \
        --atom "C1'" \
        --chain1 A --chain2 A

Run with -h for all options.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

try:
    from Bio.PDB import MMCIFParser, Superimposer
    from Bio.Align import PairwiseAligner
except ImportError:
    print("ERROR: Biopython is required. Install with: pip install biopython",
          file=sys.stderr)
    sys.exit(1)


# --------------------------------------------------------------------------
# Small utilities
# --------------------------------------------------------------------------

def check_executable(name: str) -> Optional[str]:
    """Return the resolved path of an executable, or None if missing."""
    return shutil.which(name)


def run(cmd: List[str], **kwargs) -> subprocess.CompletedProcess:
    """Run a subprocess, raising a helpful error if it fails."""
    proc = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed ({' '.join(cmd)}):\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    return proc


# --------------------------------------------------------------------------
# Step 1: x3dna-dssr wrapper
# --------------------------------------------------------------------------

@dataclass
class DssrResult:
    cif_path: str
    json_path: str
    txt_path: str = ""
    dbn: str = ""
    sequence: str = ""
    chain_id: str = ""
    nt_ids: List[str] = field(default_factory=list)   # e.g. "A.U1"
    num_pairs: int = 0
    num_stems: int = 0
    num_helices: int = 0
    num_hairpins: int = 0
    num_multiplets: int = 0
    pseudoknot_order: int = 0


def run_dssr(cif_path: str, outdir: str, dssr_bin: str = "x3dna-dssr") -> DssrResult:
    """Run x3dna-dssr on a single .cif file (JSON mode) and return its JSON path."""
    if check_executable(dssr_bin) is None:
        raise FileNotFoundError(
            f"'{dssr_bin}' not found on PATH. Install/license x3dna-dssr "
            f"(https://x3dna.org) or pass --dssr-path."
        )

    os.makedirs(outdir, exist_ok=True)
    base = os.path.splitext(os.path.basename(cif_path))[0]
    json_path = os.path.join(outdir, f"{base}.dssr.json")

    # NB: x3dna-dssr writes auxiliary dssr-*.pdb files into the CWD, so we
    # run it from within outdir to keep things tidy.
    cmd = [dssr_bin, f"-i={os.path.abspath(cif_path)}", "--format=mmcif",
           "--json", f"-o={os.path.abspath(json_path)}"]
    run(cmd, cwd=outdir)

    return DssrResult(cif_path=cif_path, json_path=json_path)


def run_dssr_text_summary(cif_path: str, outdir: str,
                           dssr_bin: str = "x3dna-dssr") -> str:
    """
    Run x3dna-dssr in its default (plain-text) mode -- i.e. *without*
    --json -- and return the full text summary. This is the same kind
    of report shown in the DSSR examples: it has a stable, documented
    "Secondary structures in dot-bracket notation (dbn) as a whole and
    per chain" section that is much more reliable to parse than
    reverse-engineering the JSON schema (which varies across DSSR
    builds/versions), and it correctly represents pseudoknots with
    multiple bracket types (), [], {}, etc.
    """
    if check_executable(dssr_bin) is None:
        raise FileNotFoundError(
            f"'{dssr_bin}' not found on PATH. Install/license x3dna-dssr "
            f"(https://x3dna.org) or pass --dssr-path."
        )

    os.makedirs(outdir, exist_ok=True)
    base = os.path.splitext(os.path.basename(cif_path))[0]
    txt_path = os.path.join(outdir, f"{base}.dssr.txt")

    cmd = [dssr_bin, f"-i={os.path.abspath(cif_path)}", "--format=mmcif",
           f"-o={os.path.abspath(txt_path)}"]
    run(cmd, cwd=outdir)

    with open(txt_path) as fh:
        return fh.read(), txt_path


def parse_dbn_block(text: str) -> Tuple[str, str]:
    """
    Extract (sequence, dot-bracket) from the plain-text DSSR summary's
    "Secondary structures in dot-bracket notation" section. Each entry
    there looks like:

        >label nts=N [whole]
        SEQUENCE
        DOT.BRACKET.STRING

    with one entry for the whole complex ("[whole]") and one per chain
    ("[chain]"). We prefer the "[whole]" entry (the correct choice even
    for single-chain structures, since some builds omit "[whole]" and
    only emit "[chain]" -- handled by the fallback below).
    """
    entries = re.findall(
        r"^>(\S+)\s+nts=(\d+)[^\[\n]*\[(\w+)\][^\n]*\n(\S+)\n(\S+)",
        text, re.MULTILINE)
    if not entries:
        return "", ""

    for _label, _nts, tag, seq, dbn in entries:
        if tag == "whole":
            return seq, dbn

    # Fallback: no explicit "[whole]" tag found (e.g. some DSSR builds
    # only emit a "[chain]" entry for single-chain structures) -- use
    # the first entry found.
    _label, _nts, _tag, seq, dbn = entries[0]
    return seq, dbn


def parse_pseudoknot_order(text: str) -> int:
    """Return the pseudoknot order reported by DSSR, or 0 if none."""
    m = re.search(r"contains (\d+)-order pseudoknot", text)
    return int(m.group(1)) if m else 0


def parse_dssr_json(result: DssrResult) -> DssrResult:
    """Populate a DssrResult with summary counts from the JSON output."""
    with open(result.json_path) as fh:
        data = json.load(fh)

    nts = data.get("nts", [])
    result.nt_ids = [nt.get("nt_id", "") for nt in nts]
    if result.nt_ids:
        result.chain_id = result.nt_ids[0].split(".")[0]

    result.num_pairs = len(data.get("pairs", []))
    result.num_stems = len(data.get("stems", []))
    result.num_helices = len(data.get("helices", []))
    result.num_hairpins = len(data.get("hairpins", []))
    result.num_multiplets = len(data.get("multiplets", []))

    return result


def get_dssr_summary(cif_path: str, outdir: str, dssr_bin: str) -> DssrResult:
    """
    Run x3dna-dssr twice: once in --json mode (for structured counts of
    pairs/stems/helices/etc.) and once in its default plain-text mode
    (for a reliable, pseudoknot-aware dot-bracket + sequence). Merges
    both into a single DssrResult.
    """
    json_result = run_dssr(cif_path, outdir, dssr_bin)
    json_result = parse_dssr_json(json_result)

    text, txt_path = run_dssr_text_summary(cif_path, outdir, dssr_bin)
    json_result.txt_path = txt_path
    json_result.sequence, json_result.dbn = parse_dbn_block(text)
    json_result.pseudoknot_order = parse_pseudoknot_order(text)

    if not json_result.sequence:
        # last-resort fallback: build the sequence directly from the
        # JSON nts list if even the text-based parse came up empty
        with open(json_result.json_path) as fh:
            data = json.load(fh)
        json_result.sequence = "".join(
            nt.get("nt_code", "?") for nt in data.get("nts", []))

    return json_result


# --------------------------------------------------------------------------
# Step 2: Biopython structural alignment + RMSD
# --------------------------------------------------------------------------

@dataclass
class RmsdResult:
    n_matched: int
    n_total_1: int
    n_total_2: int
    rmsd: float
    atom_name: str


def load_structure(cif_path: str, struct_id: str):
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure(struct_id, cif_path)


def get_chain(structure, chain_id: Optional[str]):
    model = next(iter(structure))
    if chain_id:
        return model[chain_id]
    # default: first chain in the model
    return next(iter(model))


def residue_sequence(chain) -> Tuple[List, str]:
    """Return (list of residue objects, one-letter sequence string)."""
    residues = [res for res in chain if res.id[0] == " "]  # skip het/waters
    one_letter = {
        "A": "A", "C": "C", "G": "G", "U": "U", "T": "T",
        "DA": "A", "DC": "C", "DG": "G", "DT": "T",
    }
    seq = "".join(one_letter.get(res.resname.strip(), "N") for res in residues)
    return residues, seq


def align_sequences(seq1: str, seq2: str) -> List[Tuple[int, int]]:
    """
    Pairwise-align two sequences and return a list of (i, j) index pairs
    for columns where both sequences have an aligned, identical residue.
    Indices are 0-based positions into seq1 / seq2 respectively.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.match_score = 2
    aligner.mismatch_score = -1

    alignment = aligner.align(seq1, seq2)[0]
    pairs = []
    for (s1, e1), (s2, e2) in zip(*alignment.aligned):
        for offset in range(e1 - s1):
            i, j = s1 + offset, s2 + offset
            if seq1[i] == seq2[j]:
                pairs.append((i, j))
    return pairs


def compute_rmsd(cif1: str, cif2: str,
                  chain1: Optional[str], chain2: Optional[str],
                  atom_name: str = "C1'") -> RmsdResult:
    """
    Superimpose structure 2 onto structure 1 using matched residues
    (from a sequence alignment) and a chosen representative atom per
    residue (default C1', a good sequence-independent backbone marker
    for nucleic acids; use "CA" for proteins).
    """
    s1 = load_structure(cif1, "struct1")
    s2 = load_structure(cif2, "struct2")

    c1 = get_chain(s1, chain1)
    c2 = get_chain(s2, chain2)

    res1, seq1 = residue_sequence(c1)
    res2, seq2 = residue_sequence(c2)

    matched_idx = align_sequences(seq1, seq2)

    atoms1, atoms2 = [], []
    for i, j in matched_idx:
        r1, r2 = res1[i], res2[j]
        if atom_name in r1 and atom_name in r2:
            atoms1.append(r1[atom_name])
            atoms2.append(r2[atom_name])

    if len(atoms1) < 3:
        raise RuntimeError(
            f"Only {len(atoms1)} matched '{atom_name}' atoms found between "
            f"the two chains -- need at least 3 for superposition. Check "
            f"--chain1/--chain2/--atom options."
        )

    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    sup.apply(s2.get_atoms())  # move structure 2 onto structure 1 in place

    return RmsdResult(
        n_matched=len(atoms1),
        n_total_1=len(res1),
        n_total_2=len(res2),
        rmsd=sup.rms,
        atom_name=atom_name,
    ), s2


def save_superimposed(structure, out_path: str):
    """Write the (now superimposed) structure 2 to disk as a new .cif."""
    from Bio.PDB.mmcifio import MMCIFIO
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(out_path)


# --------------------------------------------------------------------------
# Step 3: PyMOL cross-check + visualization
# --------------------------------------------------------------------------

# This is written out as a standalone .py file and executed *inside*
# PyMOL's own embedded interpreter via `pymol -cq -r <script> -- <args>`.
# That works with any command-line PyMOL install (open-source, SBGrid
# module, commercial, etc.) without requiring the separate `pymol2`
# python package to be importable in the calling environment.
_PYMOL_WORKER_SCRIPT = r'''
import sys
import json
from pymol import cmd

cif1, cif2, png_out, pse_out, json_out = sys.argv[1:6]

cmd.load(cif1, "structA")
cmd.load(cif2, "structB")

cmd.hide("everything")
cmd.show("cartoon")
cmd.color("skyblue", "structA")
cmd.color("salmon", "structB")
cmd.bg_color("black")
cmd.set("cartoon_ring_mode", 3)   # nice nucleic-acid ring rendering
cmd.set("cartoon_ring_finder", 1)

# Force a fully opaque black background in the rendered PNG (instead of a
# transparent one), and make sure nothing in the scene is see-through by
# zeroing out every transparency-related setting that could apply to the
# rendered representations.
cmd.set("ray_opaque_background", 1)
cmd.set("cartoon_transparency", 0)
cmd.set("surface_transparency", 0)
cmd.set("stick_transparency", 0)
cmd.set("sphere_transparency", 0)
cmd.set("transparency", 0)
cmd.set("ray_shadows", 1)
cmd.set("depth_cue", 0)  # no fog fade-to-background on distant atoms

# 'super' is sequence-independent structural superposition -- appropriate
# even if numbering/sequence differ slightly between the two inputs.
result = cmd.super("structB", "structA")
rmsd_pymol = result[0] if result else None

cmd.orient()
cmd.zoom(buffer=5)
cmd.ray(1600, 1200)
cmd.png(png_out, dpi=300)
cmd.save(pse_out)

with open(json_out, "w") as fh:
    json.dump({"rmsd": rmsd_pymol}, fh)
'''


def pymol_super_and_render(cif1: str, cif2: str, outdir: str,
                            png_out: str, pse_out: str,
                            pymol_bin: str = "pymol") -> Optional[float]:
    """
    Independently superimpose the two structures in PyMOL (run headlessly
    as a subprocess: `pymol -cq -r <worker script>`) using the
    sequence-independent `super` command, save an overlay image (PNG,
    opaque black background) and a PyMOL session (.pse) for interactive
    inspection.

    Returns the RMSD reported by PyMOL's `super`, or None if the `pymol`
    executable isn't available (in which case this step is skipped and
    only the Biopython RMSD is reported).
    """
    if check_executable(pymol_bin) is None:
        print(f"WARNING: '{pymol_bin}' not found on PATH -- skipping PyMOL "
              f"visualization step. Pass --pymol-path if it's installed "
              f"under a different name/location (e.g. an SBGrid module).",
              file=sys.stderr)
        return None

    os.makedirs(outdir, exist_ok=True)

    worker_path = os.path.join(outdir, "_pymol_worker.py")
    json_out = os.path.join(outdir, "_pymol_result.json")
    with open(worker_path, "w") as fh:
        fh.write(_PYMOL_WORKER_SCRIPT)

    cmd_list = [
        pymol_bin, "-cq", "-r", worker_path, "--",
        os.path.abspath(cif1), os.path.abspath(cif2),
        os.path.abspath(png_out), os.path.abspath(pse_out),
        os.path.abspath(json_out),
    ]
    run(cmd_list)

    if not os.path.exists(json_out):
        print("WARNING: PyMOL ran but produced no result file -- check "
              f"{worker_path} output manually.", file=sys.stderr)
        return None

    with open(json_out) as fh:
        data = json.load(fh)
    return data.get("rmsd")


# --------------------------------------------------------------------------
# Step 4: dot-bracket comparison / reporting
# --------------------------------------------------------------------------

def print_dbn_comparison(name1: str, r1: DssrResult, name2: str, r2: DssrResult):
    print(f"\n{'-'*78}")
    print("Secondary structure (dot-bracket notation)")
    print(f"{'-'*78}")
    pk1 = f", pseudoknot order {r1.pseudoknot_order}" if r1.pseudoknot_order else ""
    pk2 = f", pseudoknot order {r2.pseudoknot_order}" if r2.pseudoknot_order else ""
    print(f"{name1} ({len(r1.sequence)} nt, {r1.num_pairs} bp, "
          f"{r1.num_stems} stems, {r1.num_helices} helices, "
          f"{r1.num_hairpins} hairpins, {r1.num_multiplets} multiplets{pk1})")
    print(r1.sequence)
    print(r1.dbn)
    print()
    print(f"{name2} ({len(r2.sequence)} nt, {r2.num_pairs} bp, "
          f"{r2.num_stems} stems, {r2.num_helices} helices, "
          f"{r2.num_hairpins} hairpins, {r2.num_multiplets} multiplets{pk2})")
    print(r2.sequence)
    print(r2.dbn)

    if len(r1.dbn) == len(r2.dbn) and r1.dbn and r2.dbn:
        diff = "".join("." if a == b else "^" for a, b in zip(r1.dbn, r2.dbn))
        n_diff = diff.count("^")
        print(f"\nPer-position agreement ('^' marks a differing bracket/dot; "
              f"{n_diff}/{len(diff)} positions differ):")
        print(diff)
    else:
        print("\n(Sequences differ in length -- position-by-position dot-bracket "
              "diff skipped; the residue-level RMSD alignment above still "
              "handles length differences via pairwise sequence alignment.)")


def write_report(path: str, name1, name2, dssr1: DssrResult, dssr2: DssrResult,
                  rmsd_bio: RmsdResult, rmsd_pymol: Optional[float]):
    with open(path, "w") as fh:
        fh.write("RNA structure comparison report\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(f"Structure 1: {name1}\n")
        fh.write(f"Structure 2: {name2}\n\n")

        fh.write("-- Secondary structure (x3dna-dssr) --\n")
        pk1 = f", pseudoknot order {dssr1.pseudoknot_order}" if dssr1.pseudoknot_order else ""
        pk2 = f", pseudoknot order {dssr2.pseudoknot_order}" if dssr2.pseudoknot_order else ""
        fh.write(f"{name1}: {len(dssr1.sequence)} nt, {dssr1.num_pairs} bp, "
                  f"{dssr1.num_stems} stems, {dssr1.num_helices} helices{pk1}\n")
        fh.write(f"  seq: {dssr1.sequence}\n  dbn: {dssr1.dbn}\n")
        fh.write(f"{name2}: {len(dssr2.sequence)} nt, {dssr2.num_pairs} bp, "
                  f"{dssr2.num_stems} stems, {dssr2.num_helices} helices{pk2}\n")
        fh.write(f"  seq: {dssr2.sequence}\n  dbn: {dssr2.dbn}\n\n")

        fh.write("-- RMSD (Biopython Superimposer, sequence-aligned atoms) --\n")
        fh.write(f"  atom used       : {rmsd_bio.atom_name}\n")
        fh.write(f"  matched residues: {rmsd_bio.n_matched} "
                  f"(of {rmsd_bio.n_total_1} / {rmsd_bio.n_total_2} total)\n")
        fh.write(f"  RMSD            : {rmsd_bio.rmsd:.3f} A\n\n")

        fh.write("-- RMSD (PyMOL 'super', sequence-independent) --\n")
        if rmsd_pymol is not None:
            fh.write(f"  RMSD            : {rmsd_pymol:.3f} A\n")
        else:
            fh.write("  (skipped -- pymol2 module not available)\n")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Compare two RNA/DNA .cif structures: DSSR secondary "
                     "structure, Biopython RMSD, and PyMOL visualization.")
    ap.add_argument("cif1", help="First .cif file")
    ap.add_argument("cif2", help="Second .cif file")
    ap.add_argument("--outdir", default="rna_comparison_out",
                     help="Directory for all outputs (default: %(default)s)")
    ap.add_argument("--chain1", default=None,
                     help="Chain ID to use in cif1 (default: first chain)")
    ap.add_argument("--chain2", default=None,
                     help="Chain ID to use in cif2 (default: first chain)")
    ap.add_argument("--atom", default="C1'",
                     help="Representative atom per residue for RMSD "
                          "(default: %(default)s; use CA for protein)")
    ap.add_argument("--dssr-path", default="x3dna-dssr",
                     help="Path to the x3dna-dssr executable "
                          "(default: %(default)s, i.e. must be on PATH)")
    ap.add_argument("--pymol-path", default="pymol",
                     help="Path/name of the PyMOL executable "
                          "(default: %(default)s, i.e. must be on PATH -- "
                          "works with any command-line PyMOL install, "
                          "including SBGrid modules)")
    ap.add_argument("--skip-pymol", action="store_true",
                     help="Skip the PyMOL visualization/RMSD step")
    ap.add_argument("--skip-dssr", action="store_true",
                     help="Skip the x3dna-dssr secondary-structure step")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    name1 = os.path.splitext(os.path.basename(args.cif1))[0]
    name2 = os.path.splitext(os.path.basename(args.cif2))[0]

    # ---- Step 1: DSSR ----
    dssr1 = dssr2 = None
    if not args.skip_dssr:
        print(f"[1/3] Running x3dna-dssr on {args.cif1} and {args.cif2} ...")
        try:
            dssr1 = get_dssr_summary(args.cif1, os.path.join(args.outdir, "dssr"),
                                      args.dssr_path)
            dssr2 = get_dssr_summary(args.cif2, os.path.join(args.outdir, "dssr"),
                                      args.dssr_path)
            print_dbn_comparison(name1, dssr1, name2, dssr2)
        except (FileNotFoundError, RuntimeError) as exc:
            print(f"WARNING: DSSR step failed/skipped: {exc}", file=sys.stderr)
            dssr1 = dssr2 = None
    else:
        print("[1/3] Skipping x3dna-dssr step (--skip-dssr).")

    # ---- Step 2: Biopython RMSD ----
    print(f"\n[2/3] Computing RMSD with Biopython "
          f"(atom='{args.atom}') ...")
    rmsd_result, superimposed_s2 = compute_rmsd(
        args.cif1, args.cif2, args.chain1, args.chain2, args.atom)
    print(f"  Matched residues : {rmsd_result.n_matched} "
          f"(of {rmsd_result.n_total_1} / {rmsd_result.n_total_2})")
    print(f"  RMSD             : {rmsd_result.rmsd:.3f} A")

    superimposed_path = os.path.join(args.outdir, f"{name2}_superimposed.cif")
    save_superimposed(superimposed_s2, superimposed_path)
    print(f"  Superimposed structure 2 written to: {superimposed_path}")

    # ---- Step 3: PyMOL visualization + cross-check RMSD ----
    rmsd_pymol = None
    if not args.skip_pymol:
        print(f"\n[3/3] Rendering overlay with PyMOL ...")
        png_out = os.path.join(args.outdir, f"{name1}_vs_{name2}.png")
        pse_out = os.path.join(args.outdir, f"{name1}_vs_{name2}.pse")
        rmsd_pymol = pymol_super_and_render(args.cif1, args.cif2, args.outdir,
                                             png_out, pse_out, args.pymol_path)
        if rmsd_pymol is not None:
            print(f"  PyMOL 'super' RMSD: {rmsd_pymol:.3f} A")
            print(f"  Image  : {png_out}")
            print(f"  Session: {pse_out}")
    else:
        print("\n[3/3] Skipping PyMOL step (--skip-pymol).")

    # ---- Report ----
    report_path = os.path.join(args.outdir, "comparison_report.txt")
    if dssr1 and dssr2:
        write_report(report_path, name1, name2, dssr1, dssr2,
                     rmsd_result, rmsd_pymol)
        print(f"\nFull report written to: {report_path}")
    else:
        print("\n(DSSR data unavailable -- report limited to RMSD values above.)")


if __name__ == "__main__":
    main()