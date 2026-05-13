#!/usr/bin/env python3
"""Task G: publication-quality PyMOL renders for WT_holo + 6 key mutants.

Each render uses the corresponding *_holo_complex.pdb (receptor + ligand + cofactor)
from 06e (WT) or 07e (mutants).
"""
from __future__ import annotations
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

import os
PROJECT = Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
OUT = PROJECT / "12_phase7" / "07_pub_renders"
LOG = PROJECT / "logs" / "v7_task_g.log"
PIPELOG = PROJECT / "pipeline.log"
PYMOL = "/opt/homebrew/bin/pymol"

WT_COMPLEX = PROJECT / "06e_docking_wt_v5" / "wt_holo_complex.pdb"
MUT_DIR = PROJECT / "07e_mut_docking_v5" / "viewer_files"
VINA_CSV = PROJECT / "07e_mut_docking_v5" / "mutant_results_v5.csv"

TARGETS = [
    "WT_holo",
    "R215A_N226A",
    "H196A",
    "R215E",
    "R50A",
    "C195A",
    "R175E_R176E",
]


def log(msg: str) -> None:
    line = f"[V7][taskG] {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as fh: fh.write(line+"\n")
    with PIPELOG.open("a") as fh: fh.write(line+"\n")


def parse_mut_residues(name: str) -> list[int]:
    out = []
    for tok in name.split("_"):
        m = re.match(r"^[A-Z](\d+)[A-Z]$", tok)
        if m:
            out.append(int(m.group(1)))
    return out


def render(target: str, dvina: float | None) -> Path:
    if target == "WT_holo":
        cpx = WT_COMPLEX
        mut_resi = []
        title = f"WT (holo) | top affinity = {dvina:+.2f} kcal/mol" if dvina is not None else "WT (holo)"
    else:
        cpx = MUT_DIR / f"{target}_holo_complex.pdb"
        mut_resi = parse_mut_residues(target)
        title = f"{target} | dVina = {dvina:+.2f} kcal/mol" if dvina is not None else target
    if not cpx.exists():
        log(f"missing {cpx}")
        return Path("")

    out_png = OUT / f"{target}_pub.png"
    sel_mut = " or ".join([f"resi {r}" for r in mut_resi]) if mut_resi else "none"

    # ligand resn UMP/DUP variants — check in v5 it's UMP per top_pose.  Use resn UMP.
    pml = f"""
load {cpx}, cpx
remove resn HOH
bg_color white

# detect ligand residue id
hide everything
# Receptor
select receptor, polymer.protein
show cartoon, receptor
color grey80, receptor
set cartoon_transparency, 0.55, receptor

# Active site (within 5 A of UMP/dUMP)
select ligand, resn UMP+DUP+UPN+DUM
show sticks, ligand
color cyan, ligand
util.cnc ligand

# Cofactor D16
select cof, resn D16+CB3+THF+MTX
show sticks, cof
color teal, cof
util.cnc cof

# Active site residues within 5 A of ligand
select active5, byres (receptor within 5 of ligand)
show sticks, active5
color slate, active5
util.cnc active5
label active5 and name CA, "%s%s*" % (oneletter, resi)

# Mutated residue (sticks, pink)
select mut_res, ({sel_mut}) and receptor
show sticks, mut_res
color magenta, mut_res
util.cnc mut_res

# H-bonds and close contacts: distances <3.5 A between ligand and active5
distance hbonds, (ligand and not elem H), (active5 and not elem H), 3.5
color yellow, hbonds
hide labels, hbonds

# Title via PyMOL set_title
set label_size, 14
set label_color, black
set label_outline_color, white

# View
orient ligand, animate=0
zoom ligand, 7
set ray_shadows, 0
set ambient, 0.35
set ray_trace_mode, 1
set antialias, 2
viewport 1600, 1200

# Header strip via cmd.set_title (only shows in some pymol forks; use a label too)
pseudoatom title_anchor, pos=[0,0,0]
label title_anchor, "{title}"
hide labels, title_anchor

ray 1600, 1200
png {out_png}
"""
    pml_path = OUT / f"_{target}.pml"
    pml_path.write_text(pml)
    proc = subprocess.run([PYMOL, "-cq", str(pml_path)], capture_output=True, text=True, timeout=600)
    log(f"{target} render rc={proc.returncode}, png exists={out_png.exists()}")
    if proc.returncode != 0:
        (OUT / f"{target}_FAIL.log").write_text(proc.stdout + "\n--STDERR--\n" + proc.stderr)
    return out_png


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    log("Task G start")
    vina = pd.read_csv(VINA_CSV, comment="#")
    holo = vina[vina.condition == "holo"]
    dvina = {row.mutant: row.delta_vina_vs_wt for _, row in holo.iterrows()}
    dvina["WT_holo"] = float(holo[holo.mutant == "WT"].top_affinity.iloc[0])

    rendered = []
    for t in TARGETS:
        d = dvina.get(t.replace("_holo",""), None)
        if t == "WT_holo": d = dvina.get("WT_holo")
        png = render(t, d)
        if png and png.exists():
            rendered.append(str(png))

    (OUT / "summary.json").write_text(json.dumps({"rendered": rendered}, indent=2))
    log(f"Task G done: {len(rendered)} images")
    return 0


if __name__ == "__main__":
    sys.exit(main())
