#!/usr/bin/env python3
"""Stage 4: Generate PyMOL renders via headless subprocess calls."""
import os, sys, subprocess, json
from datetime import datetime
import pandas as pd

PROJECT = os.path.expanduser("~/conserved_site_project")
PYMOL_DIR = os.path.join(PROJECT, "04_pymol")
STR_DIR = os.path.join(PROJECT, "03_structure")
MSA_DIR = os.path.join(PROJECT, "01_msa")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
PYMOL = os.environ.get("PYMOL", "/opt/homebrew/bin/pymol")

def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] STAGE4: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

def run_pml(script, name):
    pml_path = os.path.join(PYMOL_DIR, f"{name}.pml")
    with open(pml_path, "w") as f:
        f.write(script)
    log(f"running pymol script {name}")
    proc = subprocess.run([PYMOL, "-cq", pml_path], capture_output=True, text=True, cwd=PYMOL_DIR)
    log(f"pymol {name}: rc={proc.returncode}")
    if proc.returncode != 0 or proc.stderr:
        log(f"  stderr: {proc.stderr[:300]}")
    return proc

def main():
    os.makedirs(PYMOL_DIR, exist_ok=True)
    selected = json.load(open(os.path.join(PROJECT, "02_active_site/selected_meta.json")))["selected"]
    sel_str = "+".join(str(p) for p in selected)
    log(f"Selected residues: {selected}")

    prot = os.path.join(STR_DIR, "protein_h.pdb")
    lig = os.path.join(STR_DIR, "ligand.pdb")

    # Make a B-factor-colored copy with conservation scores
    cons = pd.read_csv(os.path.join(MSA_DIR, "conservation_scores.csv"))
    score_map = {int(r.ref_position): float(r.js_score) if pd.notna(r.js_score) else 0.0
                 for r in cons.itertuples()}
    # Read protein and rewrite B factors
    prot_cons = os.path.join(STR_DIR, "protein_cons.pdb")
    with open(prot, "r") as f, open(prot_cons, "w") as g:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    resi = int(line[22:26])
                    score = score_map.get(resi, 0.0)
                    bf = f"{score*100:6.2f}"
                    line = line[:60] + bf + line[66:]
                except ValueError:
                    pass
            g.write(line)
    log(f"wrote {prot_cons}")

    # Render 1: overview
    s1 = f"""
load {prot}, prot
load {lig}, lig
hide everything
show surface, prot
color gray80, prot
show sticks, lig
color cyan, lig
select active, prot and resi {sel_str}
show sticks, active
color orange, active
set transparency, 0.35
bg_color white
orient
zoom prot, 4
ray 1600, 1200
png 01_overview.png
"""
    run_pml(s1, "01_overview")

    # Render 2: closeup
    s2 = f"""
load {prot}, prot
load {lig}, lig
hide everything
show cartoon, prot
color gray70, prot
show sticks, lig
color cyan, lig
select active, prot and resi {sel_str}
show sticks, active
color orange, active
label active and name CA, "%s%s" % (resn, resi)
set label_size, 16
set label_color, black
bg_color white
orient active
zoom active, 8
ray 1600, 1200
png 02_closeup.png
"""
    run_pml(s2, "02_closeup")

    # Render 3: conservation cartoon
    s3 = f"""
load {prot_cons}, prot
load {lig}, lig
hide everything
show cartoon, prot
spectrum b, blue_white_red, prot
show sticks, lig
color green, lig
bg_color white
orient
zoom prot, 4
ray 1600, 1200
png 03_conservation.png
"""
    run_pml(s3, "03_conservation")

    # Render 4: cavity
    s4 = f"""
load {prot}, prot
load {lig}, lig
hide everything
show surface, prot
color gray80, prot
set surface_cavity_mode, 1
set surface_cavity_radius, 7
set surface_cavity_cutoff, -7
show sticks, lig
color cyan, lig
bg_color white
orient
zoom prot, 4
ray 1600, 1200
png 04_cavity.png
"""
    run_pml(s4, "04_cavity")

    # Verify all 4 PNGs exist and are nontrivial
    for n in ["01_overview", "02_closeup", "03_conservation", "04_cavity"]:
        p = os.path.join(PYMOL_DIR, f"{n}.png")
        if os.path.exists(p) and os.path.getsize(p) > 5000:
            log(f"OK {n}.png ({os.path.getsize(p)} bytes)")
        else:
            log(f"FAIL {n}.png (size={os.path.getsize(p) if os.path.exists(p) else 0})")
    log("Stage 4 DONE")

if __name__ == "__main__":
    main()
