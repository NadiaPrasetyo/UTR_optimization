#!/usr/bin/env python3
"""
compare_af3_candidates.py
==========================

Compare two RNA (or DNA/RNA) structure files in mmCIF format. This copy is
used as a library by mutate_3utr_ga.py (biopython_rmsd + get_dssr_summary +
check_executable are imported directly), and remains usable stand-alone
exactly as documented below.

  1. Runs `x3dna-dssr` on each .cif to derive secondary structure
     (dot-bracket notation, base pairs, helices, stems, motifs).
  2. Uses Biopython to parse both structures and compute RMSD three
     independent ways:
       a. "Biopython RMSD"    -- global (Needleman-Wunsch) pairwise
                                  sequence alignment of the DSSR-derived
                                  sequences establishes 1:1 residue
                                  correspondence end-to-end, then
                                  Bio.PDB.Superimposer performs a single
                                  rigid-body (Kabsch) least-squares fit.
       b. "Smith-Waterman 3D" -- the same Superimposer-based fit, but
                                  residue correspondence comes from a
                                  *local* (Smith-Waterman) alignment
                                  instead, so only the best-matching
                                  local region drives the superposition.
       c. "PyMOL"             -- an independent, sequence-*independent*
                                  structural superposition via PyMOL's
                                  `super` command (see step 3).
  3. Uses PyMOL (via a headless `pymol -cq` subprocess) to independently
     superimpose the two structures (sequence-independent, `super`
     command), render a side-by-side / overlay image, and save a .pse
     session so the result can be inspected interactively. It also
     renders each structure on its own, colored by per-residue pLDDT
     confidence, and saves each of those as its own standalone PyMOL
     session (see "pLDDT confidence coloring" below).
  4. Prints a side-by-side dot-bracket comparison and a text summary
     report; writes everything to an output directory.

pLDDT confidence coloring
--------------------------
If a `molecule_confidences.json` file (AlphaFold3-style: top-level
"atom_chain_ids" and "atom_plddts" arrays, one entry per atom, in the
same order as the atoms appear in the corresponding .cif) is found next
to a structure (or passed explicitly via --confidences1/--confidences2),
its per-atom pLDDT values are:

  - stamped onto the B-factor column of the loaded structure (both for
    the Biopython-saved superimposed .cif and for the PyMOL render), and
  - used to color each structure with the standard discrete AlphaFold
    confidence scheme:

        Very high  (pLDDT > 90)              blue
        Confident  (90 >= pLDDT > 70)         turquoise
        Low        (70 >= pLDDT > 50)         yellow
        Very low   (pLDDT <= 50)              red

  A per-structure confidence breakdown (mean pLDDT + % of atoms in each
  bucket) is also printed and written to the report.

  Each confidence-colored structure is saved as its own standalone
  PyMOL session (`<name>_plddt.pse`) rather than a plain .cif. Plain
  mmCIF has no field for a viewer's display color -- only numeric data
  like B-factor survives -- so a bare .cif reopened in PyMOL would NOT
  show the coloring. A .pse session, by contrast, stores the actual
  object state (colors, representations, camera), so opening it in
  PyMOL shows exactly the coloring used for the render, and you're free
  to move it around, edit it further, or re-export from there.

Requirements
------------
  - x3dna-dssr executable on PATH (licensed binary, see https://x3dna.org)
  - PyMOL with the `pymol` executable on PATH (open-source PyMOL via
    conda:
        conda install -c conda-forge pymol-open-source
    or the commercial PyMOL)
  - Biopython  (pip install biopython)

Usage
-----
    python3 compare_af3_candidates.py structA.cif structB.cif \
        --outdir results/ \
        --atom "C1'" \
        --chain1 A --chain2 A \
        --confidences1 structA_dir/molecule_confidences.json \
        --confidences2 structB_dir/molecule_confidences.json

If --confidences1/--confidences2 are omitted, the script looks for a
file named `molecule_confidences.json` alongside each .cif (as well as
a couple of common AF3-style naming variants) and uses it automatically
if found.

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
    import numpy as np
except ImportError:
    print("ERROR: Biopython and numpy are required. "
          "Install with: pip install biopython numpy",
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
# pLDDT confidence loading / coloring
# --------------------------------------------------------------------------

PLDDT_BUCKETS = [
    ("Very high (pLDDT > 90)",        lambda v: v > 90,             "blue"),
    ("Confident (70 < pLDDT <= 90)",  lambda v: 70 < v <= 90,       "cyan"),
    ("Low (50 < pLDDT <= 70)",        lambda v: 50 < v <= 70,       "yellow"),
    ("Very low (pLDDT <= 50)",        lambda v: v <= 50,            "red"),
]


def plddt_color(value: float) -> str:
    for _label, predicate, color in PLDDT_BUCKETS:
        if predicate(value):
            return color
    return "gray"


def find_confidences_json(cif_path: str, explicit: Optional[str]) -> Optional[str]:
    if explicit:
        if not os.path.exists(explicit):
            raise FileNotFoundError(f"--confidences path not found: {explicit}")
        return explicit

    cif_dir = os.path.dirname(os.path.abspath(cif_path))
    base = os.path.splitext(os.path.basename(cif_path))[0].replace("_model", "")
    candidates = [
        os.path.join(cif_dir, "molecule_confidences.json"),
        os.path.join(cif_dir, f"{base}_confidences.json"),
        os.path.join(cif_dir, "confidences.json"),
    ]
    for c in candidates:
        if os.path.exists(c):
            print(f"Found molecule_confidences.json: {c}")
            return c
    return None


def load_plddt_confidences(path: str) -> Tuple[List[str], List[float]]:
    with open(path) as fh:
        data = json.load(fh)
    chain_ids = data.get("atom_chain_ids", [])
    plddts = data.get("atom_plddts", [])
    if len(chain_ids) != len(plddts):
        raise ValueError(
            f"{path}: atom_chain_ids ({len(chain_ids)}) and atom_plddts "
            f"({len(plddts)}) lengths do not match.")
    return chain_ids, [float(v) for v in plddts]


def summarize_plddt(plddts: List[float]) -> Dict:
    if not plddts:
        return {}
    n = len(plddts)
    counts = {label: 0 for label, _pred, _color in PLDDT_BUCKETS}
    for v in plddts:
        for label, predicate, _color in PLDDT_BUCKETS:
            if predicate(v):
                counts[label] += 1
                break
    return {
        "n_atoms": n,
        "mean": sum(plddts) / n,
        "buckets": {label: (c, 100.0 * c / n) for label, c in counts.items()},
    }


def format_plddt_summary(name: str, summary: Dict) -> str:
    if not summary:
        return f"{name}: (no confidences loaded)"
    lines = [f"{name}: mean pLDDT = {summary['mean']:.1f} "
             f"({summary['n_atoms']} atoms)"]
    for label, (count, pct) in summary["buckets"].items():
        lines.append(f"    {label:<28} {count:>6d}  ({pct:5.1f}%)")
    return "\n".join(lines)


def apply_plddt_bfactors(structure, plddts: List[float]) -> bool:
    atoms = list(structure.get_atoms())
    if len(atoms) != len(plddts):
        print(f"WARNING: confidences file has {len(plddts)} atoms but "
              f"the structure has {len(atoms)} atoms -- skipping B-factor "
              f"assignment (atom-order correspondence not reliable).",
              file=sys.stderr)
        return False
    for atom, plddt in zip(atoms, plddts):
        atom.set_bfactor(float(plddt))
    return True


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
    nt_ids: List[str] = field(default_factory=list)
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

    cmd = [dssr_bin, f"-i={os.path.abspath(cif_path)}", "--format=mmcif",
           "--json", f"-o={os.path.abspath(json_path)}"]
    run(cmd, cwd=outdir)

    return DssrResult(cif_path=cif_path, json_path=json_path)


def run_dssr_text_summary(cif_path: str, outdir: str,
                           dssr_bin: str = "x3dna-dssr"):
    """
    Run x3dna-dssr in its default (plain-text) mode and return the full
    text summary, which has a stable "Secondary structures in dot-bracket
    notation" section that correctly represents pseudoknots.
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
    entries = re.findall(
        r"^>(\S+)\s+nts=(\d+)[^\[\n]*\[(\w+)\][^\n]*\n(\S+)\n(\S+)",
        text, re.MULTILINE)
    if not entries:
        return "", ""

    for _label, _nts, tag, seq, dbn in entries:
        if tag == "whole":
            return seq, dbn

    _label, _nts, _tag, seq, dbn = entries[0]
    return seq, dbn


def parse_pseudoknot_order(text: str) -> int:
    m = re.search(r"contains (\d+)-order pseudoknot", text)
    return int(m.group(1)) if m else 0


def parse_dssr_json(result: DssrResult) -> DssrResult:
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
    Run x3dna-dssr twice (structured JSON counts + plain-text dot-bracket)
    and merge both into a single DssrResult, including result.dbn — the
    dot-bracket secondary structure string used for reporting.
    """
    json_result = run_dssr(cif_path, outdir, dssr_bin)
    json_result = parse_dssr_json(json_result)

    text, txt_path = run_dssr_text_summary(cif_path, outdir, dssr_bin)
    json_result.txt_path = txt_path
    json_result.sequence, json_result.dbn = parse_dbn_block(text)
    json_result.pseudoknot_order = parse_pseudoknot_order(text)

    if not json_result.sequence:
        with open(json_result.json_path) as fh:
            data = json.load(fh)
        json_result.sequence = "".join(
            nt.get("nt_code", "?") for nt in data.get("nts", []))

    return json_result


# --------------------------------------------------------------------------
# Step 2: RMSD calculations
# --------------------------------------------------------------------------

@dataclass
class AlignmentSummary:
    method: str
    rmsd: float
    identity_pct: float
    n_equivalent: int
    seq_len_1: int
    seq_len_2: int
    modeled_1: int
    modeled_2: int
    note: str = ""


def load_structure(cif_path: str, struct_id: str):
    parser = MMCIFParser(QUIET=True)
    return parser.get_structure(struct_id, cif_path)


def get_chain(structure, chain_id: Optional[str]):
    model = next(iter(structure))
    if chain_id:
        return model[chain_id]
    return next(iter(model))


def residue_sequence(chain) -> Tuple[List, str]:
    residues = [res for res in chain if res.id[0] == " "]
    one_letter = {
        "A": "A", "C": "C", "G": "G", "U": "U", "T": "T",
        "DA": "A", "DC": "C", "DG": "G", "DT": "T",
    }
    seq = "".join(one_letter.get(res.resname.strip(), "N") for res in residues)
    return residues, seq


def _aligned_atom_pairs(res1, seq1, res2, seq2, alignment, atom_name):
    idx_pairs = []
    for (a_start, a_end), (b_start, b_end) in zip(*alignment.aligned):
        for offset in range(a_end - a_start):
            idx_pairs.append((a_start + offset, b_start + offset))

    atoms1, atoms2, n_identical = [], [], 0
    for i, j in idx_pairs:
        r1, r2 = res1[i], res2[j]
        if atom_name not in r1 or atom_name not in r2:
            continue
        atoms1.append(r1[atom_name])
        atoms2.append(r2[atom_name])
        if seq1[i] == seq2[j]:
            n_identical += 1

    return atoms1, atoms2, n_identical


def biopython_rmsd(cif1: str, cif2: str,
                    chain1: Optional[str], chain2: Optional[str],
                    atom_name: str = "C1'",
                    plddts1: Optional[List[float]] = None,
                    plddts2: Optional[List[float]] = None):
    """
    'Biopython RMSD': a global (Needleman-Wunsch) pairwise alignment of
    the two chains' sequences establishes 1:1 residue correspondence
    end-to-end, then a single rigid-body (Kabsch) superposition
    (Bio.PDB.Superimposer) is fit over every aligned residue pair.
    Returns (AlignmentSummary, Superimposer, structure2).
    """
    s1 = load_structure(cif1, "bp_struct1")
    s2 = load_structure(cif2, "bp_struct2")

    if plddts1:
        apply_plddt_bfactors(s1, plddts1)
    if plddts2:
        apply_plddt_bfactors(s2, plddts2)

    c1 = get_chain(s1, chain1)
    c2 = get_chain(s2, chain2)
    res1, seq1 = residue_sequence(c1)
    res2, seq2 = residue_sequence(c2)

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(seq1, seq2)[0]

    atoms1, atoms2, n_identical = _aligned_atom_pairs(
        res1, seq1, res2, seq2, alignment, atom_name)

    if len(atoms1) < 3:
        raise RuntimeError(
            f"Biopython RMSD: only {len(atoms1)} aligned '{atom_name}' "
            f"atoms found between the two chains -- need at least 3. "
            f"Check --chain1/--chain2/--atom options."
        )

    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)

    n_equiv = len(atoms1)
    modeled_1 = sum(1 for r in res1 if atom_name in r)
    modeled_2 = sum(1 for r in res2 if atom_name in r)

    summary = AlignmentSummary(
        method="Biopython RMSD",
        rmsd=sup.rms,
        identity_pct=100.0 * n_identical / n_equiv,
        n_equivalent=n_equiv,
        seq_len_1=len(res1),
        seq_len_2=len(res2),
        modeled_1=modeled_1,
        modeled_2=modeled_2,
    )
    return summary, sup, s2


def smith_waterman_3d(cif1: str, cif2: str,
                       chain1: Optional[str], chain2: Optional[str],
                       atom_name: str = "C1'",
                       plddts1: Optional[List[float]] = None,
                       plddts2: Optional[List[float]] = None):
    s1 = load_structure(cif1, "sw_struct1")
    s2 = load_structure(cif2, "sw_struct2")

    if plddts1:
        apply_plddt_bfactors(s1, plddts1)
    if plddts2:
        apply_plddt_bfactors(s2, plddts2)

    c1 = get_chain(s1, chain1)
    c2 = get_chain(s2, chain2)
    res1, seq1 = residue_sequence(c1)
    res2, seq2 = residue_sequence(c2)

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(seq1, seq2)[0]

    atoms1, atoms2, n_identical = _aligned_atom_pairs(
        res1, seq1, res2, seq2, alignment, atom_name)

    if len(atoms1) < 3:
        raise RuntimeError(
            f"Smith-Waterman 3D: only {len(atoms1)} aligned '{atom_name}' "
            f"atoms found between the two chains -- need at least 3. "
            f"Check --chain1/--chain2/--atom options."
        )

    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)

    n_equiv = len(atoms1)
    modeled_1 = sum(1 for r in res1 if atom_name in r)
    modeled_2 = sum(1 for r in res2 if atom_name in r)

    summary = AlignmentSummary(
        method="Smith-Waterman 3D",
        rmsd=sup.rms,
        identity_pct=100.0 * n_identical / n_equiv,
        n_equivalent=n_equiv,
        seq_len_1=len(res1),
        seq_len_2=len(res2),
        modeled_1=modeled_1,
        modeled_2=modeled_2,
    )
    return summary, sup, s2


def save_superimposed(structure, out_path: str):
    from Bio.PDB.mmcifio import MMCIFIO
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(out_path)


def format_alignment_table(summaries: List[AlignmentSummary]) -> str:
    cols = [
        ("Method", 22, "<"), ("RMSD(A)", 8, ">"), ("Ident%", 7, ">"),
        ("EquivRes", 9, ">"), ("SeqLen1", 8, ">"), ("SeqLen2", 8, ">"),
        ("Modeled1", 9, ">"), ("Modeled2", 9, ">"),
    ]
    header = "".join(f"{name:{align}{w}}" for name, w, align in cols)
    lines = [header, "-" * len(header)]
    for s in summaries:
        row = (
            f"{s.method:<22}{s.rmsd:>8.3f}{s.identity_pct:>7.1f}"
            f"{s.n_equivalent:>9d}{s.seq_len_1:>8d}{s.seq_len_2:>8d}"
            f"{s.modeled_1:>9d}{s.modeled_2:>9d}"
        )
        lines.append(row)
    return "\n".join(lines)


# --------------------------------------------------------------------------
# Step 4: dot-bracket comparison / reporting (unused by the GA, kept for
# stand-alone use of this script)
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
              "diff skipped; the residue-level RMSD alignments above still "
              "handle length differences via pairwise sequence alignment.)")


def write_report(path: str, name1, name2, dssr1: DssrResult, dssr2: DssrResult,
                  alignment_summaries: List[AlignmentSummary],
                  plddt_summary1: Optional[Dict] = None,
                  plddt_summary2: Optional[Dict] = None):
    with open(path, "w") as fh:
        fh.write("RNA structure comparison report\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(f"Structure 1: {name1}\n")
        fh.write(f"Structure 2: {name2}\n\n")

        if dssr1 and dssr2:
            fh.write("-- Secondary structure (x3dna-dssr) --\n")
            pk1 = f", pseudoknot order {dssr1.pseudoknot_order}" if dssr1.pseudoknot_order else ""
            pk2 = f", pseudoknot order {dssr2.pseudoknot_order}" if dssr2.pseudoknot_order else ""
            fh.write(f"{name1}: {len(dssr1.sequence)} nt, {dssr1.num_pairs} bp, "
                      f"{dssr1.num_stems} stems, {dssr1.num_helices} helices{pk1}\n")
            fh.write(f"  seq: {dssr1.sequence}\n  dbn: {dssr1.dbn}\n")
            fh.write(f"{name2}: {len(dssr2.sequence)} nt, {dssr2.num_pairs} bp, "
                      f"{dssr2.num_stems} stems, {dssr2.num_helices} helices{pk2}\n")
            fh.write(f"  seq: {dssr2.sequence}\n  dbn: {dssr2.dbn}\n\n")

        fh.write("-- RMSD summary (Needleman-Wunsch global, Smith-Waterman 3D local) --\n")
        fh.write("(Ident% over equivalent residues; Modeled = residues with a "
                  "usable representative atom)\n\n")
        fh.write(format_alignment_table(alignment_summaries) + "\n")
        for s in alignment_summaries:
            if s.note:
                fh.write(f"  [{s.method}] {s.note}\n")
        fh.write("\n")

        fh.write("-- pLDDT confidence (molecule_confidences.json) --\n")
        fh.write("  Bins: Very high (>90) blue | Confident (70-90) turquoise | "
                  "Low (50-70) yellow | Very low (<50) red\n\n")
        if plddt_summary1:
            fh.write(format_plddt_summary(name1, plddt_summary1) + "\n")
        else:
            fh.write(f"{name1}: (no confidences loaded)\n")
        if plddt_summary2:
            fh.write(format_plddt_summary(name2, plddt_summary2) + "\n")
        else:
            fh.write(f"{name2}: (no confidences loaded)\n")


# --------------------------------------------------------------------------
# Main (stand-alone CLI use; the GA imports the functions above directly
# and does not call this)
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Compare two RNA/DNA .cif structures: DSSR secondary "
                     "structure and Needleman-Wunsch / Smith-Waterman 3D "
                     "Biopython RMSD (no PyMOL).")
    ap.add_argument("cif1", help="First .cif file")
    ap.add_argument("cif2", help="Second .cif file")
    ap.add_argument("--outdir", default="rna_comparison_out",
                     help="Directory for all outputs (default: %(default)s)")
    ap.add_argument("--chain1", default=None)
    ap.add_argument("--chain2", default=None)
    ap.add_argument("--atom", default="C1'")
    ap.add_argument("--confidences1", default=None)
    ap.add_argument("--confidences2", default=None)
    ap.add_argument("--skip-plddt-coloring", action="store_true")
    ap.add_argument("--dssr-path", default="x3dna-dssr")
    ap.add_argument("--skip-dssr", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    name1 = os.path.splitext(os.path.basename(args.cif1))[0]
    name2 = os.path.splitext(os.path.basename(args.cif2))[0]

    conf_path1 = conf_path2 = None
    plddts1 = plddts2 = None
    plddt_summary1 = plddt_summary2 = None
    if not args.skip_plddt_coloring:
        try:
            conf_path1 = find_confidences_json(args.cif1, args.confidences1)
            if conf_path1:
                _chain_ids1, plddts1 = load_plddt_confidences(conf_path1)
                plddt_summary1 = summarize_plddt(plddts1)
        except (FileNotFoundError, ValueError) as exc:
            print(f"WARNING: could not load confidences for {name1}: {exc}", file=sys.stderr)
        try:
            conf_path2 = find_confidences_json(args.cif2, args.confidences2)
            if conf_path2:
                _chain_ids2, plddts2 = load_plddt_confidences(conf_path2)
                plddt_summary2 = summarize_plddt(plddts2)
        except (FileNotFoundError, ValueError) as exc:
            print(f"WARNING: could not load confidences for {name2}: {exc}", file=sys.stderr)

    dssr1 = dssr2 = None
    if not args.skip_dssr:
        try:
            dssr1 = get_dssr_summary(args.cif1, os.path.join(args.outdir, "dssr"), args.dssr_path)
            dssr2 = get_dssr_summary(args.cif2, os.path.join(args.outdir, "dssr"), args.dssr_path)
            print_dbn_comparison(name1, dssr1, name2, dssr2)
        except (FileNotFoundError, RuntimeError) as exc:
            print(f"WARNING: DSSR step failed/skipped: {exc}", file=sys.stderr)
            dssr1 = dssr2 = None

    bp_summary, bp_sup, bp_s2 = biopython_rmsd(
        args.cif1, args.cif2, args.chain1, args.chain2, args.atom,
        plddts1=plddts1, plddts2=plddts2)
    alignment_summaries = [bp_summary]
    print(format_alignment_table([bp_summary]))

    sw3d_summary, sw_sup, sw_s2 = smith_waterman_3d(
        args.cif1, args.cif2, args.chain1, args.chain2, args.atom,
        plddts1=plddts1, plddts2=plddts2)
    alignment_summaries.append(sw3d_summary)
    print(format_alignment_table([sw3d_summary]))

    sw_sup.apply(sw_s2.get_atoms())
    superimposed_path = os.path.join(args.outdir, f"{name2}_superimposed.cif")
    save_superimposed(sw_s2, superimposed_path)
    print(f"  Superimposed structure 2 (Smith-Waterman 3D fit) written to: {superimposed_path}")

    report_path = os.path.join(args.outdir, f"{name1}_vs_{name2}_comparison_report.txt")
    write_report(report_path, name1, name2, dssr1, dssr2,
                 alignment_summaries, plddt_summary1, plddt_summary2)
    print(f"\nFull report written to: {report_path}")


if __name__ == "__main__":
    main()