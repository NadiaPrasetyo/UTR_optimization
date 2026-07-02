#!/usr/bin/env python3
"""
compare_af3_candidates.py

Compare a set of AlphaFold3-predicted candidate structures against a single
AlphaFold3-predicted "patented" reference structure, using AF3's standard
output directory layout (see AF3 docs: <job_name>_model.cif, <job_name>_data.json).

For every candidate this computes:

  1. RMSD          - structural distance to the patented structure, after a
                      sequence-guided global superposition (Kabsch algorithm
                      over matched, non-gap aligned residues).
  2. Seq. difference - percent sequence difference (100 - percent identity)
                      from the patented sequence, from the same alignment.
  3. cmscore        - bit score of the candidate sequence against an Infernal
                      covariance model (.cm) of the patented RNA family.
                      Requires Infernal (`cmsearch`) and a --cm-model file.
  4. MFE difference - |MFE(candidate) - MFE(patented)|, minimum free energy
                      of the secondary structure computed with ViennaRNA's
                      RNAfold ("similar to the base" = closeness to the
                      patented sequence's own MFE).

All four metrics are min-max normalized across the candidate set and
combined into a single composite score with EQUAL weights (0.25 each by
default, configurable). Composite score is scaled to [0, 1], where 1 =
most similar to the patented reference across all four metrics, 0 = least
similar. If a metric can't be computed for every candidate (e.g. Infernal
or RNAfold isn't installed, or no --cm-model was given), it is dropped and
the remaining weights are renormalized automatically, with a warning.

Requirements
------------
- Python packages: biopython  (pip install biopython --break-system-packages)
- Optional external tools (only needed for the corresponding metric):
    * RNAfold   (ViennaRNA)   -> for MFE
    * cmsearch  (Infernal)    -> for cmscore (also need a .cm model file)

Expected input layout
----------------------
This is designed to plug directly into AF3 output directories:

  patented-dir/
      <job_name>_model.cif
      <job_name>_data.json      # AF3's echo of the input, has the sequence(s)

  candidates-dir/
      candidate_1/
          candidate_1_model.cif
          candidate_1_data.json
      candidate_2/
          ...

(this is exactly what `prepare_af3_jobs.py` + AF3 produce: candidates-dir ==
 <work-dir>/af3_output/, with one subdirectory per job name.)

Usage
-----
python3 compare_af3_candidates.py \\
    --patented-dir  /path/to/af3_output/patented_job \\
    --candidates-dir /path/to/af3_output \\
    --exclude patented_job \\
    --cm-model /path/to/family.cm \\
    --out comparison.csv
"""

import argparse
import csv
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio.PDB import MMCIFParser, Superimposer
from Bio.Align import PairwiseAligner, substitution_matrices

# --------------------------------------------------------------------------
# Residue -> one-letter code
# --------------------------------------------------------------------------
THREE_TO_ONE = {
    # RNA
    "A": "A", "C": "C", "G": "G", "U": "U",
    # DNA
    "DA": "A", "DC": "C", "DG": "G", "DT": "T",
    # Protein (20 standard)
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# preferred representative atom per residue class, in priority order
NUCLEIC_ATOM_PRIORITY = ["C1'", "C1*", "P"]
PROTEIN_ATOM_PRIORITY = ["CA"]


def resname_one_letter(resname: str):
    return THREE_TO_ONE.get(resname.strip().upper())


def representative_atom(residue):
    resname = residue.get_resname().strip().upper()
    priorities = PROTEIN_ATOM_PRIORITY if resname in (
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    ) else NUCLEIC_ATOM_PRIORITY
    for atom_name in priorities:
        if residue.has_id(atom_name):
            return residue[atom_name]
    # fallback: first atom in the residue
    atoms = list(residue.get_atoms())
    return atoms[0] if atoms else None


# --------------------------------------------------------------------------
# Structure loading
# --------------------------------------------------------------------------
def find_model_cif(af3_output_dir: Path):
    """Find the top-ranking <job_name>_model.cif in an AF3 output dir."""
    candidates = sorted(af3_output_dir.glob("*_model.cif"))
    if not candidates:
        # fall back to any .cif directly in the dir (not in seed-*/sample-* subdirs)
        candidates = sorted(af3_output_dir.glob("*.cif"))
    if not candidates:
        raise FileNotFoundError(f"No *_model.cif found in {af3_output_dir}")
    return candidates[0]


def find_data_json(af3_output_dir: Path):
    matches = sorted(af3_output_dir.glob("*_data.json"))
    return matches[0] if matches else None


def load_sequences_from_data_json(data_json_path: Path):
    """Returns list of (chain_id, sequence, molecule_type)."""
    data = json.loads(data_json_path.read_text())
    out = []
    for entry in data.get("sequences", []):
        for mol_type, info in entry.items():
            chain_id = info.get("id")
            seq = info.get("sequence")
            if isinstance(chain_id, list):  # AF3 allows id lists for identical copies
                for cid in chain_id:
                    out.append((cid, seq, mol_type))
            else:
                out.append((chain_id, seq, mol_type))
    return out


def load_chains_from_cif(cif_path: Path):
    """Returns dict chain_id -> {"seq": str, "residues": [Bio.PDB residue,...]}"""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(cif_path.stem, str(cif_path))
    model = next(structure.get_models())
    chains = {}
    for chain in model:
        seq_chars = []
        residues = []
        for residue in chain:
            code = resname_one_letter(residue.get_resname())
            if code is None:
                continue  # skip waters/heteroatoms/ligands we can't map
            seq_chars.append(code)
            residues.append(residue)
        if residues:
            chains[chain.id] = {"seq": "".join(seq_chars), "residues": residues}
    return chains


# --------------------------------------------------------------------------
# Sequence alignment
# --------------------------------------------------------------------------
def make_aligner(is_protein: bool):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    if is_protein:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
    else:
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -6
        aligner.extend_gap_score = -0.5
    return aligner


def align_chains(seq_a, seq_b):
    """Global-align two sequences. Returns (aligned_a, aligned_b) as equal-length strings with '-' gaps."""
    looks_protein = any(c not in "ACGUT-" for c in (seq_a + seq_b).upper())
    aligner = make_aligner(looks_protein)
    alignment = aligner.align(seq_a, seq_b)[0]
    # Index the Alignment object directly (not str() parsing, which line-wraps
    # for long sequences and breaks under naive splitting).
    aligned_a, aligned_b = str(alignment[0]), str(alignment[1])
    return aligned_a, aligned_b


def match_residue_pairs(residues_a, seq_a, residues_b, seq_b):
    """
    Align seq_a vs seq_b, and return parallel lists of matched (non-gap,
    non-mismatch-agnostic - all aligned columns with a residue on both
    sides count) residue objects from residues_a and residues_b, plus
    counts of (matches, aligned_columns) for identity calculation.
    """
    aligned_a, aligned_b = align_chains(seq_a, seq_b)
    ia = ib = 0
    matched_a, matched_b = [], []
    n_match, n_aligned_cols = 0, 0
    for ca, cb in zip(aligned_a, aligned_b):
        a_has = ca != "-"
        b_has = cb != "-"
        if a_has and b_has:
            matched_a.append(residues_a[ia])
            matched_b.append(residues_b[ib])
            n_aligned_cols += 1
            if ca == cb:
                n_match += 1
        elif a_has and not b_has:
            n_aligned_cols += 1
        elif b_has and not a_has:
            n_aligned_cols += 1
        if a_has:
            ia += 1
        if b_has:
            ib += 1
    return matched_a, matched_b, n_match, n_aligned_cols


def greedy_match_chains(patented_chains, candidate_chains):
    """
    Pair up patented chain IDs with candidate chain IDs by best sequence
    identity (simple greedy matching - fine for the common monomer /
    small-complex case).
    """
    pairs = []
    remaining_cand = dict(candidate_chains)
    for p_id, p_info in patented_chains.items():
        best_id, best_score = None, -1
        for c_id, c_info in remaining_cand.items():
            aligned_a, aligned_b = align_chains(p_info["seq"], c_info["seq"])
            score = sum(1 for a, b in zip(aligned_a, aligned_b) if a == b and a != "-")
            if score > best_score:
                best_score, best_id = score, c_id
        if best_id is not None:
            pairs.append((p_id, best_id))
            del remaining_cand[best_id]
    return pairs


# --------------------------------------------------------------------------
# RMSD + sequence identity across (possibly multi-chain) structures
# --------------------------------------------------------------------------
def compute_rmsd_and_identity(patented_chains, candidate_chains):
    chain_pairs = greedy_match_chains(patented_chains, candidate_chains)
    if not chain_pairs:
        return None, None

    all_atoms_p, all_atoms_c = [], []
    total_match, total_cols = 0, 0

    for p_id, c_id in chain_pairs:
        p_info, c_info = patented_chains[p_id], candidate_chains[c_id]
        matched_res_p, matched_res_c, n_match, n_cols = match_residue_pairs(
            p_info["residues"], p_info["seq"], c_info["residues"], c_info["seq"]
        )
        total_match += n_match
        total_cols += n_cols
        for res_p, res_c in zip(matched_res_p, matched_res_c):
            atom_p = representative_atom(res_p)
            atom_c = representative_atom(res_c)
            if atom_p is not None and atom_c is not None:
                all_atoms_p.append(atom_p)
                all_atoms_c.append(atom_c)

    if len(all_atoms_p) < 3:
        return None, (100.0 * (1 - total_match / total_cols) if total_cols else None)

    sup = Superimposer()
    sup.set_atoms(all_atoms_p, all_atoms_c)
    rmsd = sup.rms
    seq_diff_pct = 100.0 * (1 - total_match / total_cols) if total_cols else None
    return rmsd, seq_diff_pct


# --------------------------------------------------------------------------
# MFE via ViennaRNA RNAfold
# --------------------------------------------------------------------------
_RNAFOLD_WARNED = False


def compute_mfe(sequence: str):
    global _RNAFOLD_WARNED
    if shutil.which("RNAfold") is None:
        if not _RNAFOLD_WARNED:
            print("WARNING: RNAfold not found on PATH - skipping MFE metric.", file=sys.stderr)
            _RNAFOLD_WARNED = True
        return None
    rna_seq = sequence.upper().replace("T", "U")
    try:
        result = subprocess.run(
            ["RNAfold", "--noPS"],
            input=rna_seq + "\n",
            capture_output=True,
            text=True,
            timeout=300,
        )
        # Output line 2 looks like: ....((((...)))).... (-12.30)
        for line in result.stdout.splitlines():
            if "(" in line and ")" in line and any(ch.isdigit() for ch in line):
                mfe_str = line.split("(")[-1].split(")")[0].strip()
                return float(mfe_str)
    except Exception as exc:
        print(f"WARNING: RNAfold failed: {exc}", file=sys.stderr)
    return None


# --------------------------------------------------------------------------
# cmscore via Infernal cmsearch
# --------------------------------------------------------------------------
_CMSEARCH_WARNED = False


def compute_cmscore(sequence: str, cm_model: Path, job_name: str, tmp_dir: Path):
    global _CMSEARCH_WARNED
    if cm_model is None:
        return None
    if shutil.which("cmsearch") is None:
        if not _CMSEARCH_WARNED:
            print("WARNING: cmsearch (Infernal) not found on PATH - skipping cmscore metric.", file=sys.stderr)
            _CMSEARCH_WARNED = True
        return None

    fasta_path = tmp_dir / f"{job_name}.fa"
    tblout_path = tmp_dir / f"{job_name}.tbl"
    fasta_path.write_text(f">{job_name}\n{sequence}\n")
    try:
        subprocess.run(
            ["cmsearch", "--tblout", str(tblout_path), "--noali", "-E", "10.0",
             str(cm_model), str(fasta_path)],
            capture_output=True, text=True, timeout=600, check=True,
        )
    except Exception as exc:
        print(f"WARNING: cmsearch failed for {job_name}: {exc}", file=sys.stderr)
        return None

    if not tblout_path.exists():
        return None
    best_score = None
    for line in tblout_path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        # Infernal --tblout column 15 (0-indexed 14) is the bit score
        try:
            score = float(fields[14])
        except (IndexError, ValueError):
            continue
        if best_score is None or score > best_score:
            best_score = score
    return best_score


# --------------------------------------------------------------------------
# Normalization / composite score
# --------------------------------------------------------------------------
def minmax_normalize(values, lower_is_better):
    present = [v for v in values if v is not None]
    if len(present) < 2 or max(present) == min(present):
        # can't usefully normalize -> return 1.0 for present values, None for missing
        return [(1.0 if v is not None else None) for v in values]
    lo, hi = min(present), max(present)
    norm = []
    for v in values:
        if v is None:
            norm.append(None)
            continue
        n = (v - lo) / (hi - lo)
        norm.append(1 - n if lower_is_better else n)
    return norm


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="Compare AF3 candidate structures/sequences against a patented reference "
                    "(RMSD, sequence difference, cmscore, MFE difference -> equal-weight composite).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--patented-dir", required=True, type=Path,
                    help="AF3 output directory for the patented/reference job.")
    p.add_argument("--candidates-dir", required=True, type=Path,
                    help="Directory containing one AF3 output subdirectory per candidate "
                         "(e.g. <work-dir>/af3_output/).")
    p.add_argument("--exclude", action="append", default=[],
                    help="Candidate subdirectory name(s) to skip (e.g. the patented job itself, "
                         "if it lives inside --candidates-dir too). Can repeat.")
    p.add_argument("--cm-model", type=Path, default=None,
                    help="Infernal .cm covariance model file for cmscore. Optional.")
    p.add_argument("--weights", default="0.25,0.25,0.25,0.25",
                    help="Comma-separated weights for rmsd,seq_diff,cmscore,mfe_diff. "
                         "Default equal weighting: %(default)s")
    p.add_argument("--out", type=Path, default=Path("af3_comparison.csv"),
                    help="Output CSV path. Default: %(default)s")
    args = p.parse_args()

    weight_names = ["rmsd", "seq_diff", "cmscore", "mfe_diff"]
    weights = dict(zip(weight_names, [float(w) for w in args.weights.split(",")]))
    if len(weights) != 4:
        sys.exit("ERROR: --weights must have exactly 4 comma-separated values "
                  "(rmsd,seq_diff,cmscore,mfe_diff)")

    # --- load patented reference ---
    patented_cif = find_model_cif(args.patented_dir)
    patented_chains = load_chains_from_cif(patented_cif)
    patented_data_json = find_data_json(args.patented_dir)
    if patented_data_json:
        patented_full_seq = "".join(seq for _, seq, _ in load_sequences_from_data_json(patented_data_json))
    else:
        patented_full_seq = "".join(info["seq"] for info in patented_chains.values())
    patented_mfe = compute_mfe(patented_full_seq)
    print(f"Patented reference: {patented_cif.name}  "
          f"({sum(len(c['seq']) for c in patented_chains.values())} residues, "
          f"MFE={patented_mfe})")

    # --- discover candidate job dirs ---
    candidate_dirs = sorted(
        d for d in args.candidates_dir.iterdir()
        if d.is_dir() and d.name not in args.exclude and d.resolve() != args.patented_dir.resolve()
    )
    if not candidate_dirs:
        sys.exit(f"ERROR: no candidate subdirectories found in {args.candidates_dir}")

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for cand_dir in candidate_dirs:
            job_name = cand_dir.name
            try:
                cand_cif = find_model_cif(cand_dir)
            except FileNotFoundError as e:
                print(f"WARNING: {e} - skipping {job_name}", file=sys.stderr)
                continue
            cand_chains = load_chains_from_cif(cand_cif)
            data_json = find_data_json(cand_dir)
            if data_json:
                cand_full_seq = "".join(seq for _, seq, _ in load_sequences_from_data_json(data_json))
            else:
                cand_full_seq = "".join(info["seq"] for info in cand_chains.values())

            rmsd, seq_diff_pct = compute_rmsd_and_identity(patented_chains, cand_chains)
            mfe = compute_mfe(cand_full_seq)
            mfe_diff = abs(mfe - patented_mfe) if (mfe is not None and patented_mfe is not None) else None
            cmscore = compute_cmscore(cand_full_seq, args.cm_model, job_name, tmp_dir)

            rows.append({
                "job_name": job_name,
                "rmsd": rmsd,
                "seq_diff_pct": seq_diff_pct,
                "cmscore": cmscore,
                "mfe": mfe,
                "mfe_diff": mfe_diff,
            })
            print(f"  {job_name}: RMSD={rmsd}, seq_diff={seq_diff_pct}%, "
                  f"cmscore={cmscore}, MFE={mfe} (diff={mfe_diff})")

    if not rows:
        sys.exit("ERROR: no candidates could be scored.")

    # --- normalize each metric across the candidate set ---
    norm_rmsd = minmax_normalize([r["rmsd"] for r in rows], lower_is_better=True)
    norm_seqdiff = minmax_normalize([r["seq_diff_pct"] for r in rows], lower_is_better=True)
    norm_cmscore = minmax_normalize([r["cmscore"] for r in rows], lower_is_better=False)
    norm_mfediff = minmax_normalize([r["mfe_diff"] for r in rows], lower_is_better=True)

    for i, r in enumerate(rows):
        r["norm_rmsd"] = norm_rmsd[i]
        r["norm_seq_diff"] = norm_seqdiff[i]
        r["norm_cmscore"] = norm_cmscore[i]
        r["norm_mfe_diff"] = norm_mfediff[i]

        present_weights = {}
        for key, norm_key in [("rmsd", "norm_rmsd"), ("seq_diff", "norm_seq_diff"),
                               ("cmscore", "norm_cmscore"), ("mfe_diff", "norm_mfe_diff")]:
            if r[norm_key] is not None:
                present_weights[key] = weights[key]
        total_w = sum(present_weights.values())
        if total_w == 0:
            r["composite_score"] = None
            continue
        score = 0.0
        for key, norm_key in [("rmsd", "norm_rmsd"), ("seq_diff", "norm_seq_diff"),
                               ("cmscore", "norm_cmscore"), ("mfe_diff", "norm_mfe_diff")]:
            if key in present_weights:
                score += (present_weights[key] / total_w) * r[norm_key]
        r["composite_score"] = score

    rows.sort(key=lambda r: (r["composite_score"] is None, -(r["composite_score"] or 0)))

    fieldnames = ["job_name", "rmsd", "seq_diff_pct", "cmscore", "mfe", "mfe_diff",
                  "norm_rmsd", "norm_seq_diff", "norm_cmscore", "norm_mfe_diff", "composite_score"]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"\nWrote comparison table for {len(rows)} candidate(s) to: {args.out}")
    print("composite_score: 1.0 = most similar to patented reference across all available "
          "metrics, 0.0 = least similar. (equal weights unless overridden with --weights)")
    print("\nTop candidate(s):")
    for r in rows[:5]:
        cs = r["composite_score"]
        print(f"  {r['job_name']}: composite_score={cs:.3f}" if cs is not None
              else f"  {r['job_name']}: composite_score=N/A")


if __name__ == "__main__":
    main()