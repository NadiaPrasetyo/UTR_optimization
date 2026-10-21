
import sys
import os
import json
from pymol import cmd, stored

(cif1, cif2, png_out, pse_out, json_out,
 conf1, conf2, plddt_png1, plddt_png2,
 plddt_pse1, plddt_pse2) = sys.argv[1:12]

cmd.load(cif1, "structA")
cmd.load(cif2, "structB")

cmd.hide("everything")
cmd.show("cartoon")
cmd.bg_color("white")
cmd.set("cartoon_ring_mode", 3)   # nice nucleic-acid ring rendering
cmd.set("cartoon_ring_finder", 1)
cmd.set("ray_opaque_background", 1)


def assign_plddt_bfactors(obj_name, conf_path):
    """
    Try to stamp per-atom pLDDT values from molecule_confidences.json
    (top-level "atom_plddts", in file/atom order) onto the B-factor
    column of the already-loaded PyMOL object, in atom order.

    If the JSON's atom count doesn't match what PyMOL loaded (this can
    happen even when Biopython's atom count matched fine elsewhere in
    this script, e.g. because PyMOL's mmCIF reader handles altlocs,
    symmetry/assembly records, or waters differently), this does NOT
    give up on coloring: AlphaFold-family mmCIF files already carry
    per-atom pLDDT in the B_iso_or_equiv column, which PyMOL loads into
    `b` automatically -- so we fall back to using whatever is already
    there. Coloring is only truly unavailable if there's no JSON *and*
    no usable (non-zero) B-factor signal on disk.

    Returns (has_signal: bool, source: str) where source is one of
    "json", "cif_bfactor", or "none".
    """
    n_atoms = cmd.count_atoms(obj_name)
    used_json = False

    if conf_path and conf_path != "None":
        with open(conf_path) as fh:
            data = json.load(fh)
        plddts = data.get("atom_plddts", [])
        if len(plddts) == n_atoms:
            stored.plddts = list(plddts)
            cmd.alter(obj_name, "b = stored.plddts.pop(0)")
            used_json = True
        else:
            sys.stderr.write(
                f"WARNING: {conf_path} has {len(plddts)} atoms but "
                f"PyMOL loaded {n_atoms} atoms for {obj_name} -- "
                f"falling back to the B-factor values already present "
                f"in the loaded structure instead of the JSON.\n")

    b_vals = [a.b for a in cmd.get_model(obj_name).atom]
    has_signal = any(b != 0.0 for b in b_vals)

    if used_json:
        source = "json"
    elif has_signal:
        source = "cif_bfactor"
    else:
        source = "none"
        sys.stderr.write(
            f"WARNING: no usable pLDDT signal (JSON or B-factor) found "
            f"for {obj_name} -- skipping confidence coloring.\n")

    return has_signal, source


def color_by_plddt(obj_name):
    """
    Discrete AlphaFold-style pLDDT confidence coloring:
      Very high (pLDDT > 90)             blue
      Confident (70 < pLDDT <= 90)       turquoise (aquamarine)
      Low       (50 < pLDDT <= 70)       yellow (yelloworange)
      Very low  (pLDDT <= 50)            red (deepsalmon)
    Colored from low to high so each later, narrower selection overrides
    the broader one below it. Uses PyMOL's own builtin color names
    (rather than a custom set_color palette), so the exact same names
    resolve identically whether the session is built here or re-created
    by hand later.
    """
    cmd.color("deepsalmon", obj_name)
    cmd.color("yelloworange", f"{obj_name} and b > 50")
    cmd.color("aquamarine", f"{obj_name} and b > 70")
    cmd.color("slate", f"{obj_name} and b > 90")


has_plddt_a, plddt_source_a = assign_plddt_bfactors("structA", conf1)
has_plddt_b, plddt_source_b = assign_plddt_bfactors("structB", conf2)
sys.stderr.write(f"pLDDT coloring source -- structA: {plddt_source_a}, "
                  f"structB: {plddt_source_b}\n")

# -- Render 1: overlay, colored by structure identity (for RMSD comparison) --
cmd.color("skyblue", "structA")
cmd.color("salmon", "structB")

# 'super' is sequence-independent structural superposition -- appropriate
# even if numbering/sequence differ slightly between the two inputs.
result = cmd.super("structB", "structA")
rmsd_pymol = result[0] if result else None
n_aligned_pymol = result[1] if result else None

cmd.orient()
cmd.zoom(buffer=5)
cmd.ray(1600, 1200)
cmd.png(png_out, dpi=300)
cmd.save(pse_out)

# -- Render 2 & 3: each structure individually, colored by pLDDT confidence --
# Each is also saved as its own standalone PyMOL session (.pse). Unlike
# a plain .cif -- which has no field for a viewer's display color, only
# numeric data like B-factor -- a .pse session stores the actual object
# state (colors, representations, camera), so reopening it in PyMOL
# reproduces the coloring exactly, and it can still be freely edited,
# re-colored, or re-exported from there.
if has_plddt_a:
    color_by_plddt("structA")
    cmd.disable("structB")
    cmd.orient("structA")
    cmd.zoom("structA", buffer=5)
    cmd.ray(1600, 1200)
    cmd.png(plddt_png1, dpi=300)
    cmd.save(plddt_pse1, "structA")
    cmd.enable("structB")

if has_plddt_b:
    color_by_plddt("structB")
    cmd.disable("structA")
    cmd.orient("structB")
    cmd.zoom("structB", buffer=5)
    cmd.ray(1600, 1200)
    cmd.png(plddt_png2, dpi=300)
    cmd.save(plddt_pse2, "structB")
    cmd.enable("structA")

with open(json_out, "w") as fh:
    json.dump({
        "rmsd": rmsd_pymol,
        "n_aligned": n_aligned_pymol,
        "plddt_colored_A": has_plddt_a,
        "plddt_colored_B": has_plddt_b,
        "plddt_source_A": plddt_source_a,
        "plddt_source_B": plddt_source_b,
    }, fh)
