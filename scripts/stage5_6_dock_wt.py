#!/usr/bin/env python3
"""Stage 5: Ligand PDBQT prep. Stage 6: WT docking with vina."""
import os, sys, subprocess, json, re, math
from datetime import datetime
import numpy as np
from Bio.PDB import PDBParser

PROJECT = os.path.expanduser("~/conserved_site_project")
LIG_DIR = os.path.join(PROJECT, "05_ligand")
DOCK_DIR = os.path.join(PROJECT, "06_docking_wt")
STR_DIR = os.path.join(PROJECT, "03_structure")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
OBABEL = os.environ.get("OBABEL", "/Users/ario/conserved_site_project/.venv/bin/obabel")
VINA = os.environ.get("VINA", "/opt/homebrew/bin/vina")
PYMOL = os.environ.get("PYMOL", "/opt/homebrew/bin/pymol")
VENV_PY = os.path.join(PROJECT, ".venv/bin/python")
MK_PREP_REC = os.path.join(PROJECT, ".venv/bin/mk_prepare_receptor.py")
MK_PREP_LIG = os.path.join(PROJECT, ".venv/bin/mk_prepare_ligand.py")

def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] STAGE5_6: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

def stage5():
    os.makedirs(LIG_DIR, exist_ok=True)
    log("Stage 5: ligand prep")
    lig_h = os.path.join(STR_DIR, "ligand_h.pdb")
    lig_pdbqt = os.path.join(LIG_DIR, "ligand.pdbqt")

    # Try obabel first
    proc = subprocess.run([OBABEL, lig_h, "-O", lig_pdbqt,
                           "--partialcharge", "gasteiger", "-p", "7.4"],
                          capture_output=True, text=True)
    log(f"obabel ligand→pdbqt rc={proc.returncode}")
    if proc.returncode != 0 or not os.path.exists(lig_pdbqt) or os.path.getsize(lig_pdbqt) < 100:
        log(f"obabel main failed: {proc.stderr[:300]}; trying meeko fallback")
        # meeko fallback
        if os.path.exists(MK_PREP_LIG):
            proc = subprocess.run([VENV_PY, MK_PREP_LIG, "-i", lig_h, "-o", lig_pdbqt],
                                  capture_output=True, text=True)
            log(f"meeko ligand rc={proc.returncode} stderr={proc.stderr[:200]}")
    # Validate atom types
    with open(lig_pdbqt) as f:
        content = f.read()
    has_unk = " UNK " in content
    log(f"ligand.pdbqt size={os.path.getsize(lig_pdbqt)}, has_UNK={has_unk}")
    log("Stage 5 DONE")
    return lig_pdbqt

def prepare_receptor(prot_in, prot_out, label):
    """Try meeko first, fall back to obabel."""
    if os.path.exists(MK_PREP_REC):
        proc = subprocess.run([VENV_PY, MK_PREP_REC, "-i", prot_in, "-o", prot_out, "--no-flexible"],
                              capture_output=True, text=True)
        log(f"meeko receptor {label} rc={proc.returncode}")
        if proc.returncode == 0 and os.path.exists(prot_out) and os.path.getsize(prot_out) > 1000:
            return True
        log(f"meeko failed: {proc.stderr[:300]}")
    proc = subprocess.run([OBABEL, prot_in, "-O", prot_out, "-xr"],
                          capture_output=True, text=True)
    log(f"obabel receptor {label} rc={proc.returncode}")
    if proc.returncode != 0 or not os.path.exists(prot_out) or os.path.getsize(prot_out) < 1000:
        log(f"obabel receptor failed: {proc.stderr[:300]}")
        return False
    return True

def compute_centroid(pdb_path, residues):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("p", pdb_path)
    coords = []
    for chain in s[0]:
        for res in chain:
            if res.id[0] == " " and res.id[1] in residues and "CA" in res:
                coords.append(res["CA"].get_coord())
    arr = np.array(coords)
    return arr.mean(axis=0), len(arr)

def parse_vina_stdout(stdout):
    """Parse the affinity table from vina stdout."""
    affinities = []
    in_table = False
    for line in stdout.split("\n"):
        if re.match(r"^\s*-+\+-+\+-+\+-+\s*$", line) or re.match(r"^-----+\+-----", line):
            in_table = True
            continue
        if in_table:
            # row format: "   1       -7.729          0          0"
            m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$", line)
            if m:
                affinities.append(float(m.group(2)))
            elif line.strip() == "":
                if affinities:
                    break
            elif not re.match(r"^\s*\d", line):
                if affinities:
                    break
    return affinities

def parse_pdbqt_models(pdbqt_path):
    """Parse vina output PDBQT into a list of model atom-coord dicts."""
    models = []
    cur = None
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                cur = []
            elif line.startswith("ENDMDL"):
                if cur is not None:
                    models.append(cur)
                cur = None
            elif line.startswith(("ATOM  ", "HETATM")) and cur is not None:
                try:
                    name = line[12:16].strip()
                    elem = line[76:78].strip() if len(line) >= 78 else name[0]
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    if elem != "H" and not name.startswith("H"):
                        cur.append((name, np.array([x, y, z])))
                except Exception:
                    pass
    return models

def native_heavy_coords(lig_pdb):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("l", lig_pdb)
    coords = []
    for atom in s.get_atoms():
        if atom.element != "H":
            coords.append((atom.get_name(), atom.get_coord()))
    return coords

def rmsd_to_native(pose_atoms, native_atoms):
    # Match by atom name first; fall back to nearest-pair if names differ
    nat_map = {n: c for n, c in native_atoms}
    matched = []
    for n, c in pose_atoms:
        if n in nat_map:
            matched.append((c, nat_map[n]))
    if len(matched) < max(3, len(native_atoms) // 2):
        # nearest-pair greedy match using full distance matrix
        pose_arr = np.array([c for _, c in pose_atoms])
        nat_arr = np.array([c for _, c in native_atoms])
        if len(pose_arr) == 0 or len(nat_arr) == 0:
            return float("nan")
        # use min(pose, native) pairs
        n = min(len(pose_arr), len(nat_arr))
        used_n = set()
        pairs = []
        for i in range(len(pose_arr)):
            best = None; best_d = 1e9
            for j in range(len(nat_arr)):
                if j in used_n: continue
                d = np.linalg.norm(pose_arr[i] - nat_arr[j])
                if d < best_d:
                    best_d = d; best = j
            if best is not None:
                used_n.add(best)
                pairs.append((pose_arr[i], nat_arr[best]))
            if len(pairs) >= n: break
        matched = pairs
    if not matched:
        return float("nan")
    sq = [np.sum((a-b)**2) for a, b in matched]
    return float(math.sqrt(sum(sq) / len(sq)))

def stage6(lig_pdbqt):
    os.makedirs(DOCK_DIR, exist_ok=True)
    log("Stage 6: WT docking")

    selected = json.load(open(os.path.join(PROJECT, "02_active_site/selected_meta.json")))["selected"]
    prot_h = os.path.join(STR_DIR, "protein_h.pdb")
    rec_pdbqt = os.path.join(DOCK_DIR, "protein_wt.pdbqt")
    ok = prepare_receptor(prot_h, rec_pdbqt, "WT")
    if not ok:
        log("FATAL: receptor prep failed")
        sys.exit(1)

    centroid, n = compute_centroid(prot_h, set(selected))
    log(f"centroid from {n} CA atoms: {centroid}")
    cx, cy, cz = float(centroid[0]), float(centroid[1]), float(centroid[2])

    out_pdbqt = os.path.join(DOCK_DIR, "wt_poses.pdbqt")
    log_file = os.path.join(DOCK_DIR, "vina_wt.log")
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", "22", "--size_y", "22", "--size_z", "22",
           "--exhaustiveness", "16", "--num_modes", "20", "--seed", "42",
           "--out", out_pdbqt, "--cpu", "4"]
    log(f"vina cmd: {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_file, "w") as f:
        f.write("STDOUT\n" + proc.stdout + "\nSTDERR\n" + proc.stderr)
    if proc.returncode != 0:
        log(f"vina rc={proc.returncode}, stderr={proc.stderr[:400]}")
        sys.exit(1)

    affs = parse_vina_stdout(proc.stdout)
    log(f"vina affinities: {affs[:5]}")
    top_aff = affs[0] if affs else float("nan")
    mean_top3 = float(np.mean(affs[:3])) if len(affs) >= 3 else float("nan")

    # RMSD top pose vs native dUMP
    poses = parse_pdbqt_models(out_pdbqt)
    native = native_heavy_coords(os.path.join(STR_DIR, "ligand.pdb"))
    rmsd = rmsd_to_native(poses[0], native) if poses else float("nan")
    log(f"top affinity {top_aff} kcal/mol, mean_top3 {mean_top3}, rmsd_to_native {rmsd:.2f} A")

    # Save WT result
    result = {"top_affinity": top_aff, "mean_top3": mean_top3, "rmsd_top_to_native": rmsd,
              "n_poses": len(poses), "centroid": [cx, cy, cz], "all_affinities": affs}
    with open(os.path.join(DOCK_DIR, "wt_result.json"), "w") as f:
        json.dump(result, f, indent=2)

    # PyMOL screenshot
    # split top pose to its own pdb for visualization
    top_pose_pdbqt = os.path.join(DOCK_DIR, "wt_top_pose.pdbqt")
    with open(out_pdbqt) as f, open(top_pose_pdbqt, "w") as g:
        in_first = False; done = False
        for line in f:
            if line.startswith("MODEL 1"):
                in_first = True
            if in_first and not done:
                g.write(line)
            if in_first and line.startswith("ENDMDL"):
                done = True
                break
    top_pose_pdb = os.path.join(DOCK_DIR, "wt_top_pose.pdb")
    subprocess.run([OBABEL, top_pose_pdbqt, "-O", top_pose_pdb], capture_output=True, text=True)

    pml = os.path.join(DOCK_DIR, "wt_view.pml")
    sel_str = "+".join(str(p) for p in selected)
    with open(pml, "w") as f:
        f.write(f"""
load {prot_h}, prot
load {top_pose_pdb}, pose
load {os.path.join(STR_DIR, 'ligand.pdb')}, native
hide everything
show cartoon, prot
color gray70, prot
select active, prot and resi {sel_str}
show sticks, active
color orange, active
show sticks, pose
color magenta, pose
show sticks, native
color cyan, native
bg_color white
orient pose
zoom pose, 8
ray 1600, 1200
png wt_topdock.png
""")
    subprocess.run([PYMOL, "-cq", pml], capture_output=True, text=True, cwd=DOCK_DIR)
    log(f"wrote wt_topdock.png ({os.path.getsize(os.path.join(DOCK_DIR, 'wt_topdock.png'))} bytes)")
    log("Stage 6 DONE")

def main():
    lig_pdbqt = stage5()
    stage6(lig_pdbqt)

if __name__ == "__main__":
    main()
