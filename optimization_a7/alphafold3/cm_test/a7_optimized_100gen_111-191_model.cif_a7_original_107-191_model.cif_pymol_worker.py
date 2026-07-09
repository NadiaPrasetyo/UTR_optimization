
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
cmd.set("ray_opaque_background", 1)

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

with open(json_out, "w") as fh:
    json.dump({"rmsd": rmsd_pymol, "n_aligned": n_aligned_pymol}, fh)
