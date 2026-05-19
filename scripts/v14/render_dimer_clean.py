#!/usr/bin/env python3
"""Phase 14 — re-render the dimer-aware structure images cleanly for the PPTX.

Fixes from v1:
  - Overview: don't label every residue (clutter); only catalytic + clamp residues
  - Active-site close-up: keep cartoon + sticks of the catalytic dyad + clamp;
    NO surface (mixed surface+sticks looked messy)
  - Cavity carve-out: replace with a proper FPocket cavity 18 surface
    coloured differently from the rest of the protein
"""
import subprocess, shutil
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "14_inhibitor_design" / "presentation" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
PYMOL = shutil.which("pymol") or "/opt/homebrew/bin/pymol"

APO_PDB = REPO / "06f_receptor_fixed" / "dimer_noH.pdb"
HOLO_PDB = REPO / "03_structure" / "1hvy.pdb"
POCKET18 = REPO / "14_inhibitor_design" / "04_allosteric" / "apo_for_fpocket_out" / "pockets" / "pocket18_atm.pdb"

# 1. Clean dimer overview — cartoon only, two chains in different muted greys,
#    catalytic residues as coloured spheres on Cα (no overlapping text labels)
overview_pml = f"""
reinitialize
load {APO_PDB}, prot
bg_color white
hide everything
show cartoon, prot
color grey80, prot and chain A
color grey60, prot and chain B
# catalytic + clamp residues highlighted as small spheres (NO labels — they collided in v1)
select cat_a, prot and chain A and resi 195+196+175+176+215
show spheres, cat_a and name CA
color salmon, cat_a and name CA
set sphere_scale, 1.5, cat_a and name CA
# single off-image legend label
set ray_opaque_background, 1
orient prot
zoom prot, 8
turn x, 10
ray 1600, 1200
png {OUT / 'dimer_overview_clean.png'}, dpi=160
"""

# 2. Active-site close-up — chain A only, cartoon + sticks of the dyad + clamp,
#    NO surface
closeup_pml = f"""
reinitialize
load {APO_PDB}, prot
bg_color white
hide everything
remove chain B
show cartoon, prot
color grey70, prot
select cat, prot and resi 195+196+175+176+215+226+258
show sticks, cat and not name H*
util.cnc cat
color salmon, cat and elem C
# only 3 labels — the ones the legend will refer to
label cat and resi 195 and name CA, "C195"
label cat and resi 215 and name CA, "R215"
label cat and resi 258 and name CA, "Y258"
set label_position, (3, 3, 3)
set label_size, 18
set label_color, black
set label_outline_color, white
set label_font_id, 7
set ray_opaque_background, 1
orient cat
zoom cat, 6
ray 1600, 1200
png {OUT / 'dimer_activesite_clean.png'}, dpi=160
"""

# 3. Cavity 18 carve-out — apo cartoon as transparent grey, FPocket cavity 18
#    residues as wheat surface, allosteric loop 181-197 ∩ cavity in red
cavity_pml = f"""
reinitialize
load {APO_PDB}, prot
load {POCKET18}, pocket18
bg_color white
hide everything, pocket18
show cartoon, prot
color grey80, prot
set cartoon_transparency, 0.4
# residues that belong to cavity 18
select cav, byres prot within 4.0 of pocket18
show surface, cav
set transparency, 0.10, cav
color wheat, cav
# the allosteric-loop subset (181-197)
select loop, cav and resi 189-197
color salmon, loop
# show sticks of the most-conserved contact residues only
select hot, prot and resi 54+87+190+191+196+200+201 and chain B
show sticks, hot and not name H*
color forest, hot and elem C
# only label two key residues to avoid stacking
label hot and resi 196 and name CA, "H196"
label hot and resi 200 and name CA, "F200"
set label_position, (4, 4, 4)
set label_size, 18
set label_color, black
set label_outline_color, white
set label_font_id, 7
set ray_opaque_background, 1
orient cav
zoom cav, 4
ray 1600, 1200
png {OUT / 'cavity18_carve_clean.png'}, dpi=160
"""

# 4. Bonus: holo state showing dUMP + cofactor crystal pose (real, clean)
holo_pml = f"""
reinitialize
load {HOLO_PDB}, full
bg_color white
hide everything
show cartoon, full and polymer
color grey80, full and chain A
color grey60, full and chain B
# dUMP (UMP residue)
select dump, full and resn UMP
show sticks, dump and not name H*
util.cnc dump
color cyan, dump and elem C
# raltitrexed cofactor (resname D16)
select cofac, full and resn D16
show sticks, cofac and not name H*
util.cnc cofac
color hotpink, cofac and elem C
# catalytic residues
select cat, full and chain A and resi 195+196+215 and polymer
show sticks, cat and not name H*
color salmon, cat and elem C
label cat and resi 195 and name CA, "Cys195"
label dump and name N3, "dUMP"
label cofac and name N1 and chain A, "D16"
set label_size, 16
set label_color, black
set label_outline_color, white
set label_font_id, 7
set ray_opaque_background, 1
orient dump or cofac
zoom dump or cofac, 5
ray 1600, 1200
png {OUT / 'holo_dump_cofactor_clean.png'}, dpi=160
"""

for name, script in [("overview", overview_pml),
                     ("closeup", closeup_pml),
                     ("cavity", cavity_pml),
                     ("holo", holo_pml)]:
    pml_path = OUT / f"_render_{name}.pml"
    pml_path.write_text(script)
    print(f"Rendering {name}…")
    r = subprocess.run([PYMOL, "-cq", str(pml_path)], capture_output=True, timeout=240)
    if r.returncode != 0:
        print(f"  ✗ pymol returned {r.returncode}: {r.stderr.decode()[:300]}")
    else:
        png = OUT / {
            "overview": "dimer_overview_clean.png",
            "closeup":  "dimer_activesite_clean.png",
            "cavity":   "cavity18_carve_clean.png",
            "holo":     "holo_dump_cofactor_clean.png",
        }[name]
        print(f"  ✓ {png}  ({png.stat().st_size/1024:.0f} KB)" if png.exists() else "  ✗ no png")
