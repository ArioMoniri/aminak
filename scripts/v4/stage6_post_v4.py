#!/usr/bin/env python3
"""Re-run atom-name restoration for the canonical wt_{apo,holo}_top.pdbqt
files in 06d_docking_wt_v4/ using the corrected restore_atom_names function.
Also recompute named-pose RMSD and write back into wt_{apo,holo}.json.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common_v4 import (restore_atom_names, native_heavy, rmsd_top, PROJECT)

DOCK_DIR = os.path.join(PROJECT, "06d_docking_wt_v4")
STR_V2 = os.path.join(PROJECT, "03b_structure_v2")
LIG_DIR = os.path.join(PROJECT, "05b_ligand_v2")

native_lig = os.path.join(STR_V2, "ligand.pdb")
ref_ligand_pdb = os.path.join(STR_V2, "ligand_h.pdb")
lig_pdbqt = os.path.join(LIG_DIR, "dump.pdbqt")

for cond in ("apo", "holo"):
    top_pdbqt = os.path.join(DOCK_DIR, f"wt_{cond}_top.pdbqt")
    top_pdb = os.path.join(DOCK_DIR, f"wt_{cond}_top_pose.pdb")
    restore_atom_names(top_pdbqt, ref_ligand_pdb, top_pdb,
                       reference_dock_pdbqt=lig_pdbqt)
    pose = native_heavy(top_pdb)
    crystal = native_heavy(native_lig)
    rmsd = rmsd_top(pose, crystal)
    print(f"{cond}: named-pose RMSD vs crystal = {rmsd:.3f}")
    json_path = os.path.join(DOCK_DIR, f"wt_{cond}.json")
    with open(json_path) as f:
        d = json.load(f)
    d["rmsd_top_named"] = rmsd
    with open(json_path, "w") as f:
        json.dump(d, f, indent=2, default=str)
    print(f"  updated {json_path}")
