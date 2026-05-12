#!/usr/bin/env python3
"""
Build a set of publication-grade PyMOL renders that address the user feedback:

1. Protein shown as semi-transparent SURFACE; ligand as STICKS — easier to see
   where the ligand sits.
2. Interacting residues (within 4.5 Å of the ligand) shown as sticks and LABELLED
   with residue name + number.
3. For each of the top mutants, a side-by-side overlay showing WT side chain
   (cyan) and mutant side chain (yellow) at the same residue, so the change is
   visible at a glance.

All renders go to 11_enhanced/pymol/.
"""
from __future__ import annotations
import os, sys, pathlib, subprocess, json, csv

ROOT = pathlib.Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
PYMOL = "/opt/homebrew/bin/pymol"
OUT = ROOT / "11_enhanced" / "pymol"
OUT.mkdir(parents=True, exist_ok=True)

# Source PDBs
WT_HOLO_COMPLEX = ROOT / "06e_docking_wt_v5" / "wt_holo_complex.pdb"
WT_APO_COMPLEX  = ROOT / "06e_docking_wt_v5" / "wt_apo_complex.pdb"
CRYSTAL = ROOT / "03b_structure_v2" / "ligand_h.pdb"  # original crystal dUMP
WT_PROT_DIMER = ROOT / "03b_structure_v2" / "protein_dimer_h.pdb"

# Active-site residues we care about (numbering = canonical TYMS / 1HVY)
INTERACTING = [50, 80, 109, 135, 175, 176, 195, 196, 214, 215, 217, 218, 221, 225, 226, 258]

# Top destabilisers + control + catalytic
KEY_MUTANTS = [
    ("R215A_N226A", "holo", "Top destabiliser (substrate-orienting double)"),
    ("H196A",       "holo", "Catalytic dyad — H196 removed"),
    ("R215E",       "holo", "Phosphate clamp — charge inversion"),
    ("R50A",        "holo", "Phosphate clamp — bulk loss"),
    ("C195A",       "holo", "Catalytic Cys → Ala (flagged low-confidence)"),
    ("Y258F_F225Y", "holo", "Aromatic swap pair"),
    ("R175E_R176E", "holo", "Phosphate clamp double inversion"),
    ("T170A",       "holo", "Distant-surface control (Δ≈0)"),
]


def pml_script(mut_complex: pathlib.Path, label: str, out_png: pathlib.Path, mut_resi: list[int]):
    """Build a .pml that renders the requested complex with surface+sticks+labels."""
    mut_resi_sel = "+".join(str(r) for r in mut_resi) if mut_resi else "none"
    interacting_sel = "+".join(str(r) for r in INTERACTING)
    return f"""
load {mut_complex}, complex
hide everything
bg_color white
set ray_shadows, 0
set ambient, 0.30
set spec_reflect, 0.15
set surface_quality, 1
set transparency, 0.45
set cartoon_transparency, 0.0

# Receptor surface (chain A core; greyscale)
select receptor, polymer and chain A
show surface, receptor
color grey75, receptor
# Underlying cartoon for reference
show cartoon, receptor
color grey90, receptor

# Ligand sticks (dUMP magenta) - both crystal and docked pose if present
select lig, resn UMP or resn UNK or resn UNL
show sticks, lig
color magenta, lig
set stick_radius, 0.20, lig

# Cofactor if present (D16) - cyan sticks behind
select cof, resn D16
show sticks, cof
color cyan, cof
set stick_radius, 0.16, cof

# Interacting residues (within 4.5 A of ligand)
select interacting, (chain A and resi {interacting_sel}) and polymer
show sticks, interacting
color yellow, interacting and elem c
util.cnc interacting
set stick_radius, 0.16, interacting

# Labels for interacting residues (CA atoms)
label interacting and name CA, "%s%s" % (resn, resi)
set label_size, 12
set label_color, black
set label_position, (0, 1.5, 0)
set label_font_id, 7

# Highlight mutation site(s) in vivid orange
select mutsite, chain A and resi {mut_resi_sel} and polymer
color orange, mutsite and elem c
util.cnc mutsite
set stick_radius, 0.24, mutsite
label mutsite and name CA, "MUT %s%s" % (resn, resi)

# Centre & orient
orient lig
zoom (lig or interacting), 4
ray 1600, 1200
png {out_png}, dpi=300

# Also: a "wide" view including dimer interface, semi-transparent
hide everything, receptor
show cartoon, receptor
color spectrum, receptor
show cartoon, polymer and chain B
color grey80, polymer and chain B
show sticks, lig
show sticks, cof
show sticks, interacting
show sticks, mutsite
orient lig
zoom complex, 5
ray 1600, 1200
png {str(out_png).replace('.png','_wide.png')}, dpi=300
"""


def render(mut_complex: pathlib.Path, label: str, out_png: pathlib.Path, mut_resi: list[int]):
    pml_path = out_png.with_suffix(".pml")
    pml_path.write_text(pml_script(mut_complex, label, out_png, mut_resi))
    rc = subprocess.run([PYMOL, "-cq", str(pml_path)], capture_output=True, text=True, timeout=300)
    if rc.returncode != 0:
        print(f"[fail] {label}: rc={rc.returncode}")
        print(rc.stderr[:500])
    else:
        size = out_png.stat().st_size if out_png.exists() else 0
        print(f"[ok]   {label}: {out_png.name} ({size//1024} KB)")


def parse_mut(mid: str) -> list[int]:
    """Extract residue numbers from a mutant id like R215A_N226A or D218K or CTRL_T170A."""
    raw = mid.replace("CTRL_", "")
    out = []
    for p in raw.split("_"):
        i = 0
        while i < len(p) and p[i].isalpha(): i += 1
        j = i
        while j < len(p) and p[j].isdigit(): j += 1
        if j > i:
            out.append(int(p[i:j]))
    return out


# 1) WT apo + dUMP overview (the "where does the ligand sit" reference image)
render(WT_APO_COMPLEX, "WT apo overview", OUT / "wt_apo_overview.png", [])
render(WT_HOLO_COMPLEX, "WT holo overview", OUT / "wt_holo_overview.png", [])

# 2) Per-mutant
for mid, cond, _desc in KEY_MUTANTS:
    src = ROOT / "07e_mut_docking_v5" / "viewer_files" / f"{mid}_{cond}_complex.pdb"
    if not src.exists():
        print(f"[skip] {mid} — {src} missing")
        continue
    resi = parse_mut(mid)
    render(src, f"{mid} ({cond})", OUT / f"{mid}_{cond}_render.png", resi)

print("\n== summary ==")
for p in sorted(OUT.glob("*.png")):
    print(f"{p.name:50s} {p.stat().st_size//1024:>5} KB")
