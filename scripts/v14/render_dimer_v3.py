#!/usr/bin/env python3
"""Phase 14 — v3 dimer renders, addressing the user's screenshot feedback:

  - Dimer overview: cartoon only, NO spheres (the v2 spheres at Cα were tiny and
    looked like blood spots). Use proper sticks ONLY for catalytic residues with
    label_position offsets so labels don't overlap.
  - Active-site close-up: orient so labels are visible; use 3D position offsets
    to spread labels.
  - Cavity-18 carve-out: previous version showed a dark black hole — replaced
    with a clean wheat MS surface + protein cartoon visible behind.
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

# 1. CLEAN dimer overview — cartoon + small thin sticks for catalytic dyad on chain A.
#    No labels in the image (they always collide); legend will be in the slide.
overview_pml = f"""
reinitialize
load {APO_PDB}, prot
bg_color white
hide everything
show cartoon, prot
color grey80, prot and chain A
color grey60, prot and chain B
# catalytic + clamp sticks (chain A only, thin so they don't crowd)
select cat, prot and chain A and resi 195+196+215
show sticks, cat and not name H*
color salmon, cat and elem C
util.cnc cat
set stick_radius, 0.20
set ray_opaque_background, 1
orient prot
zoom prot, 10
turn x, 12
ray 1600, 1200
png {OUT / 'dimer_overview_v3.png'}, dpi=160
"""

# 2. Active-site close-up — chain A, cartoon transparent, key residues as sticks
#    with strategically-placed labels (3D position offsets via label_position)
closeup_pml = f"""
reinitialize
load {APO_PDB}, prot
bg_color white
hide everything
remove chain B
show cartoon, prot
color grey75, prot
set cartoon_transparency, 0.25
select cat, prot and resi 195+196+215+258
show sticks, cat and not name H*
util.cnc cat
color salmon, cat and elem C
set stick_radius, 0.22
# strategically-spaced labels
label cat and resi 195 and name CA, "C195"
label cat and resi 196 and name CA, "H196"
label cat and resi 215 and name CA, "R215"
label cat and resi 258 and name CA, "Y258"
# offset labels so they don't overlap
set label_position, (3, 3, 3)
set label_size, 22
set label_color, black
set label_outline_color, white
set label_font_id, 7
set label_shadow_mode, 2
set ray_opaque_background, 1
orient cat
zoom cat, 8
ray 1600, 1200
png {OUT / 'dimer_activesite_v3.png'}, dpi=160
"""

# 3. Cavity 18 carve-out — CLEAN wheat-coloured pocket surface
#    Previous version had a dark hole because surface_quality was low.
cavity_pml = f"""
reinitialize
load {APO_PDB}, prot
load {POCKET18}, pocket18
bg_color white
hide everything
show cartoon, prot
color grey85, prot
# cavity-18 residues
select cav, byres prot within 4.0 of pocket18
# OPAQUE wheat surface
show surface, cav
color wheat, cav
set surface_quality, 2
set surface_color, wheat, cav
set transparency, 0.0, cav
set ambient, 0.55
set light_count, 4
set spec_count, 0
set spec_reflect, 0
set ray_shadow, 0
set two_sided_lighting, 0
# loop 181-197 ∩ cavity in salmon
select cavloop, cav and resi 189-197
color salmon, cavloop
# 2 key sticks for orientation, well-placed labels with offset
select hot, prot and chain B and resi 196+200
show sticks, hot and not name H*
util.cnc hot
color forest, hot and elem C
set stick_radius, 0.24
label hot and resi 200 and name CA, "F200"
label hot and resi 196 and name CA, "H196"
set label_position, (5, 3, 0)
set label_size, 22
set label_color, black
set label_outline_color, white
set label_font_id, 7
set label_shadow_mode, 2
set ray_opaque_background, 1
# orient so cavity is FRONT-facing (zoom on the surface centroid + small rotation)
center cav
zoom cav, 4
ray 1600, 1200
png {OUT / 'cavity18_carve_v3.png'}, dpi=160
"""

for name, script in [("overview", overview_pml),
                     ("closeup", closeup_pml),
                     ("cavity", cavity_pml)]:
    pml_path = OUT / f"_render_{name}_v3.pml"
    pml_path.write_text(script)
    print(f"Rendering {name}…")
    r = subprocess.run([PYMOL, "-cq", str(pml_path)], capture_output=True, timeout=240)
    png = OUT / f"dimer_{name}_v3.png" if name != "cavity" else OUT / "cavity18_carve_v3.png"
    if name == "closeup": png = OUT / "dimer_activesite_v3.png"
    if name == "overview": png = OUT / "dimer_overview_v3.png"
    if png.exists():
        print(f"  ✓ {png}  ({png.stat().st_size/1024:.0f} KB)")
    else:
        print(f"  ✗ failed: {r.stderr.decode()[:200]}")
