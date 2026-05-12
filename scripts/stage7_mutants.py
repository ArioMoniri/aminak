#!/usr/bin/env python3
"""Stage 7: Build mutation panel, mutate via PyMOL wizard, re-dock with vina."""
import os, sys, subprocess, json, re, math, csv
from datetime import datetime
import numpy as np
from Bio.PDB import PDBParser

sys.path.insert(0, os.path.join(os.path.expanduser("~/conserved_site_project"), "scripts"))
from stage5_6_dock_wt import (parse_vina_stdout, prepare_receptor, parse_pdbqt_models,
                                 native_heavy_coords, rmsd_to_native, compute_centroid)

PROJECT = os.path.expanduser("~/conserved_site_project")
MUT_DIR = os.path.join(PROJECT, "07_mut_docking")
WT_DIR = os.path.join(PROJECT, "06_docking_wt")
LIG_DIR = os.path.join(PROJECT, "05_ligand")
STR_DIR = os.path.join(PROJECT, "03_structure")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
PYMOL = os.environ.get("PYMOL", "/opt/homebrew/bin/pymol")
VINA = os.environ.get("VINA", "/opt/homebrew/bin/vina")
OBABEL = os.environ.get("OBABEL", "/Users/ario/conserved_site_project/.venv/bin/obabel")
VENV_PY = os.path.join(PROJECT, ".venv/bin/python")

def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] STAGE7: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

# 3-letter codes
AA1to3 = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
          "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
          "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}

# Chemically opposite-class single mutations (besides Ala-scan)
# Selected per spec hints
def opposite_for(aa):
    return {
        "C": "S",  # Cys→Ser (kill nucleophile, retain H-bond geometry)
        "H": "F",  # His→Phe (kill proton shuttle, hydrophobic)
        "F": "D",  # Phe→Asp (hydrophobic→charged)
        "Y": "A",  # Tyr→Ala (already covered) — pick Y→F? but both hydroxyl variants. Use Y→D.
        "G": "W",  # Gly→Trp (small→bulky)
        "D": "K",  # Asp→Lys (charge swap)
        "E": "K",
        "N": "D",  # Asn→Asp (gain charge)
        "R": "E",  # Arg→Glu (charge swap)
    }.get(aa, "A")

def pymol_mutate(in_pdb, out_pdb, resi, new_aa3):
    """Use PyMOL mutagenesis wizard to mutate single residue."""
    pml = f"""
load {in_pdb}, prot
wizard mutagenesis
refresh_wizard
cmd.get_wizard().set_mode("{new_aa3}")
cmd.get_wizard().do_select("/prot//A/{resi}/")
cmd.get_wizard().apply()
set_wizard
save {out_pdb}, prot
"""
    pml_path = out_pdb.replace(".pdb", ".pml")
    with open(pml_path, "w") as f:
        f.write(pml)
    proc = subprocess.run([PYMOL, "-cq", pml_path], capture_output=True, text=True)
    return proc

def pymol_mutate_double(in_pdb, out_pdb, mutations):
    """Apply multiple single mutations in sequence."""
    cmds = [f"load {in_pdb}, prot", "wizard mutagenesis", "refresh_wizard"]
    for resi, new_aa3 in mutations:
        cmds += [
            f'cmd.get_wizard().set_mode("{new_aa3}")',
            f'cmd.get_wizard().do_select("/prot//A/{resi}/")',
            "cmd.get_wizard().apply()",
        ]
    cmds += ["set_wizard", f"save {out_pdb}, prot"]
    pml_path = out_pdb.replace(".pdb", ".pml")
    with open(pml_path, "w") as f:
        f.write("\n".join(cmds))
    proc = subprocess.run([PYMOL, "-cq", pml_path], capture_output=True, text=True)
    return proc

def add_h_obabel(in_pdb, out_pdb):
    return subprocess.run([OBABEL, in_pdb, "-O", out_pdb, "-h"],
                          capture_output=True, text=True)

def dock(rec_pdbqt, lig_pdbqt, out_pdbqt, log_file, centroid):
    cx, cy, cz = centroid
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", "22", "--size_y", "22", "--size_z", "22",
           "--exhaustiveness", "16", "--num_modes", "20", "--seed", "42",
           "--out", out_pdbqt, "--cpu", "4"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_file, "w") as f:
        f.write("STDOUT\n" + proc.stdout + "\nSTDERR\n" + proc.stderr)
    return proc

def find_surface_residue(prot_pdb, centroid, min_dist=18.0):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("p", prot_pdb)
    centroid = np.array(centroid)
    candidates = []
    for res in s[0]["A"]:
        if res.id[0] != " ": continue
        if "CA" not in res: continue
        d = np.linalg.norm(res["CA"].get_coord() - centroid)
        if d >= min_dist:
            # we want a surface residue — pick one with hydrophilic/charged side chain
            rn = res.get_resname()
            if rn in ("LYS","GLU","ASP","ARG","GLN","ASN","SER","THR"):
                candidates.append((d, res.id[1], rn))
    candidates.sort(reverse=True)
    return candidates[0] if candidates else None

def make_screenshot(prot_pdb, top_pose_pdb, mutated_resi_list, out_png):
    sel_str = "+".join(str(p) for p in mutated_resi_list)
    pml = f"""
load {prot_pdb}, prot
load {top_pose_pdb}, pose
hide everything
show cartoon, prot
color gray70, prot
select mut, prot and resi {sel_str}
show sticks, mut
color red, mut
show sticks, pose
color magenta, pose
bg_color white
orient pose
zoom pose, 8
ray 1200, 900
png {out_png}
"""
    pml_path = out_png.replace(".png", ".pml")
    with open(pml_path, "w") as f:
        f.write(pml)
    subprocess.run([PYMOL, "-cq", pml_path], capture_output=True, text=True)

def main():
    os.makedirs(MUT_DIR, exist_ok=True)
    log("Stage 7 starting")

    selected = json.load(open(os.path.join(PROJECT, "02_active_site/selected_meta.json")))["selected"]
    log(f"Active site residues: {selected}")
    wt_result = json.load(open(os.path.join(WT_DIR, "wt_result.json")))
    wt_top = wt_result["top_affinity"]
    centroid = wt_result["centroid"]

    # parse residue letters from active_site_residues.csv
    import pandas as pd
    asdf = pd.read_csv(os.path.join(PROJECT, "02_active_site/active_site_residues.csv"))
    res_letter = {int(r.ref_position): r.residue for r in asdf.itertuples()}

    prot_h = os.path.join(STR_DIR, "protein_h.pdb")
    prot_chainA = os.path.join(STR_DIR, "protein_chainA.pdb")  # use unH for mutagenesis
    lig_pdbqt = os.path.join(LIG_DIR, "ligand.pdbqt")
    native_lig = os.path.join(STR_DIR, "ligand.pdb")
    native_atoms = native_heavy_coords(native_lig)

    # Build mutation list
    mutations = []
    for p in selected:
        aa = res_letter.get(p, "A")
        # always Ala (unless already Ala)
        if aa != "A":
            mutations.append({"id": f"{aa}{p}A", "type": "single_ala", "muts": [(p, "ALA")]})
        # opposite
        opp = opposite_for(aa)
        if opp != "A" and opp != aa:
            mutations.append({"id": f"{aa}{p}{opp}", "type": f"single_opposite",
                              "muts": [(p, AA1to3[opp])]})

    # 5 biologically meaningful doubles
    # Pick pairs of nearby/mechanistically related residues
    doubles = [
        ("C195A_H196A", [(195,"ALA"),(196,"ALA")], "double_catalytic_dyad"),
        ("R175E_R176E", [(175,"GLU"),(176,"GLU")], "double_arg_clamp_swap"),
        ("D218N_N226D", [(218,"ASN"),(226,"ASP")], "double_charge_swap"),
        ("Y258F_F225Y", [(258,"PHE"),(225,"TYR")], "double_aromatic_swap"),
        ("C195S_H196N", [(195,"SER"),(196,"ASN")], "double_polar_neutral"),
    ]
    for did, muts, dtype in doubles:
        mutations.append({"id": did, "type": dtype, "muts": muts})

    # 1 control: distant surface residue → Ala
    surf = find_surface_residue(prot_chainA, centroid, min_dist=18.0)
    if surf:
        d, ri, rn = surf
        aa1 = {v:k for k,v in AA1to3.items()}.get(rn, "X")
        mutations.append({"id": f"CTRL_{aa1}{ri}A", "type": "control_surface",
                          "muts": [(ri, "ALA")], "ctrl_dist_A": d})
        log(f"Control: {aa1}{ri}A (distance {d:.1f} Å from active centroid)")

    log(f"Total mutations to test: {len(mutations)}")

    rows = []
    # WT row first
    rows.append({
        "mutation_id": "WT", "type": "wildtype", "residues_changed": "",
        "top_affinity": wt_top, "delta_vs_wt": 0.0,
        "mean_top3": wt_result["mean_top3"],
        "rmsd_top_to_native": wt_result["rmsd_top_to_native"],
        "screenshot_path": os.path.join(WT_DIR, "wt_topdock.png"),
    })

    for i, m in enumerate(mutations):
        mid = m["id"]; muts = m["muts"]
        log(f"[{i+1}/{len(mutations)}] {mid} ({m['type']})")
        mut_pdb = os.path.join(MUT_DIR, f"{mid}.pdb")
        mut_h = os.path.join(MUT_DIR, f"{mid}_h.pdb")
        rec_pdbqt = os.path.join(MUT_DIR, f"{mid}.pdbqt")
        out_pdbqt = os.path.join(MUT_DIR, f"{mid}_poses.pdbqt")
        log_path = os.path.join(MUT_DIR, f"{mid}_vina.log")
        screenshot = os.path.join(MUT_DIR, f"{mid}.png")

        # mutate
        if len(muts) == 1:
            proc = pymol_mutate(prot_chainA, mut_pdb, muts[0][0], muts[0][1])
        else:
            proc = pymol_mutate_double(prot_chainA, mut_pdb, muts)
        if proc.returncode != 0 or not os.path.exists(mut_pdb):
            log(f"  mutation failed: rc={proc.returncode} stderr={proc.stderr[:200]}")
            continue
        # add H
        add_h_obabel(mut_pdb, mut_h)
        # prepare receptor
        ok = prepare_receptor(mut_h, rec_pdbqt, mid)
        if not ok:
            log(f"  receptor prep failed for {mid}")
            continue
        # dock
        proc = dock(rec_pdbqt, lig_pdbqt, out_pdbqt, log_path, centroid)
        if proc.returncode != 0:
            log(f"  vina failed for {mid}")
            continue
        with open(log_path) as f:
            txt = f.read().split("STDERR")[0]
        affs = parse_vina_stdout(txt)
        if not affs:
            log(f"  no affinities parsed for {mid}")
            continue
        top = affs[0]; m3 = float(np.mean(affs[:3]))
        # rmsd top
        poses = parse_pdbqt_models(out_pdbqt)
        rmsd = rmsd_to_native(poses[0], native_atoms) if poses else float("nan")
        # screenshot
        # split top pose
        top_pose_pdb = os.path.join(MUT_DIR, f"{mid}_top.pdb")
        top_pose_pdbqt = os.path.join(MUT_DIR, f"{mid}_top.pdbqt")
        with open(out_pdbqt) as f, open(top_pose_pdbqt, "w") as g:
            in_first = False; done = False
            for line in f:
                if line.startswith("MODEL 1") and not done:
                    in_first = True
                if in_first and not done:
                    g.write(line)
                if in_first and line.startswith("ENDMDL"):
                    done = True; break
        subprocess.run([OBABEL, top_pose_pdbqt, "-O", top_pose_pdb], capture_output=True, text=True)
        make_screenshot(mut_h, top_pose_pdb, [r for r,_ in muts], screenshot)

        rows.append({
            "mutation_id": mid, "type": m["type"],
            "residues_changed": ",".join(f"{r}->{a}" for r,a in muts),
            "top_affinity": top, "delta_vs_wt": top - wt_top,
            "mean_top3": m3, "rmsd_top_to_native": rmsd,
            "screenshot_path": screenshot,
        })
        log(f"  {mid}: top={top:.2f} Δ={top-wt_top:+.2f} rmsd={rmsd:.2f}")

    # Write CSV
    out_csv = os.path.join(MUT_DIR, "results_full.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    log(f"Wrote {out_csv} with {len(rows)} rows")
    log("Stage 7 DONE")

if __name__ == "__main__":
    main()
