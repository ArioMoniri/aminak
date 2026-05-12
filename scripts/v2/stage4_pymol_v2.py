#!/usr/bin/env python3
"""Stage 4 v2: PyMOL renderings on the dimer with v2 conservation b-factor coloring."""
import os, sys, subprocess, json
from datetime import datetime
import pandas as pd

PROJECT = os.path.expanduser("~/conserved_site_project")
PYM_DIR = os.path.join(PROJECT, "04b_pymol_v2")
STR_DIR = os.path.join(PROJECT, "03b_structure_v2")
MSA_DIR = os.path.join(PROJECT, "01b_msa_v2")
AS_DIR = os.path.join(PROJECT, "02b_active_site_v2")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v2_04_pymol.log")
PYMOL = os.environ.get("PYMOL", "/opt/homebrew/bin/pymol")


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V2] STAGE4: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def main():
    os.makedirs(PYM_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 4 v2 starting")

    cons = pd.read_csv(os.path.join(MSA_DIR, "conservation_scores.csv"))
    sel = json.load(open(os.path.join(AS_DIR, "selected_meta.json")))["selected"]

    # Write a per-residue b-factor altered PDB: b-factor = js_score * 100
    src_pdb = os.path.join(STR_DIR, "protein_dimer_h.pdb")
    bcol_pdb = os.path.join(PYM_DIR, "protein_dimer_jsd_bfactor.pdb")
    js_lookup = {int(r.ref_position): float(r.js_score) for _, r in cons.iterrows() if not pd.isna(r.js_score)}

    with open(src_pdb) as fin, open(bcol_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    resi = int(line[22:26])
                except ValueError:
                    fout.write(line); continue
                bval = js_lookup.get(resi, 0.0) * 100
                bval = min(99.99, max(0.0, bval))
                # b-factor field is cols 61-66 (60..66 in 0-indexed slice)
                new_line = line[:60] + f"{bval:6.2f}" + line[66:]
                fout.write(new_line)
            else:
                fout.write(line)
    log(f"wrote b-factor PDB: {bcol_pdb}")

    sel_str = "+".join(str(p) for p in sel)
    cof_a = os.path.join(STR_DIR, "cofactor_chainA_h.pdb")
    cof_b = os.path.join(STR_DIR, "cofactor_chainB_h.pdb")
    lig = os.path.join(STR_DIR, "ligand.pdb")

    # Common pml header
    common = f"""
load {bcol_pdb}, prot
load {lig}, ump
load {cof_a}, d16A
load {cof_b}, d16B
hide everything
bg_color white
spectrum b, blue_white_red, prot, 0, 30
show cartoon, prot
show sticks, ump
color cyan, ump
show sticks, d16A
color magenta, d16A
show sticks, d16B
color salmon, d16B
select active_A, prot and chain A and resi {sel_str}
select active_B, prot and chain B and resi {sel_str}
show sticks, active_A
show sticks, active_B
"""

    # 1) Dimer overview
    pml1 = os.path.join(PYM_DIR, "dimer_overview.pml")
    with open(pml1, "w") as f:
        f.write(common + """
orient prot
zoom prot, 5
ray 1600, 1200
png dimer_overview.png
""")
    # 2) Active site close-up chain A
    pml2 = os.path.join(PYM_DIR, "active_site_chainA.pml")
    with open(pml2, "w") as f:
        f.write(common + """
orient (active_A or ump or d16A)
zoom (active_A or ump or d16A), 4
ray 1600, 1200
png active_site_chainA.png
""")
    # 3) Active site close-up chain B
    pml3 = os.path.join(PYM_DIR, "active_site_chainB.pml")
    with open(pml3, "w") as f:
        f.write(common + """
orient (active_B or d16B)
zoom (active_B or d16B), 4
ray 1600, 1200
png active_site_chainB.png
""")
    # 4) Conservation surface (whole dimer, surface coloured)
    pml4 = os.path.join(PYM_DIR, "conservation_surface.pml")
    with open(pml4, "w") as f:
        f.write(f"""
load {bcol_pdb}, prot
hide everything
bg_color white
show surface, prot
spectrum b, blue_white_red, prot, 0, 30
load {lig}, ump
show sticks, ump
color cyan, ump
orient prot
zoom prot, 5
ray 1600, 1200
png conservation_surface.png
""")
    # 5) Catalytic dyad close-up (Cys195 + His196 + Arg175/176/215 + dUMP + cofactor)
    pml5 = os.path.join(PYM_DIR, "catalytic_dyad.pml")
    with open(pml5, "w") as f:
        f.write(common + """
hide sticks, active_B
select dyad, prot and chain A and resi 195+196+175+176+215+226
show sticks, dyad
label dyad and name CA, "%s%s" % (resn, resi)
orient (dyad or ump or d16A)
zoom (dyad or ump or d16A), 3
ray 1600, 1200
png catalytic_dyad.png
""")

    for pml in [pml1, pml2, pml3, pml4, pml5]:
        proc = subprocess.run([PYMOL, "-cq", pml], capture_output=True, text=True, cwd=PYM_DIR)
        png = pml.replace(".pml", ".png")
        size = os.path.getsize(png) if os.path.exists(png) else 0
        log(f"  rendered {os.path.basename(png)} rc={proc.returncode} size={size}")
        if proc.returncode != 0:
            log(f"    stderr: {proc.stderr[:300]}")

    log("Stage 4 v2 DONE")


if __name__ == "__main__":
    main()
