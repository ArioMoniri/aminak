#!/usr/bin/env python3
"""Stage 6 v2: WT docking on the dimer in two conditions:
  - apo  (cofactor D16 removed)
  - holo (cofactor D16 retained)
Tighter box (18^3 A), exhaustiveness 32, parse from REMARK VINA RESULT."""
import os, sys, subprocess, json, re, math
from datetime import datetime
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

PROJECT = os.path.expanduser("~/conserved_site_project")
DOCK_DIR = os.path.join(PROJECT, "06b_docking_wt_v2")
LIG_DIR = os.path.join(PROJECT, "05b_ligand_v2")
STR_DIR = os.path.join(PROJECT, "03b_structure_v2")
AS_DIR = os.path.join(PROJECT, "02b_active_site_v2")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v2_06_dock_wt.log")
OBABEL = os.environ.get("OBABEL", "/Users/ario/conserved_site_project/.venv/bin/obabel")
VINA = os.environ.get("VINA", "/opt/homebrew/bin/vina")
PYMOL = os.environ.get("PYMOL", "/opt/homebrew/bin/pymol")
VENV_PY = os.path.join(PROJECT, ".venv/bin/python")
MK_PREP_REC = os.path.join(PROJECT, ".venv/bin/mk_prepare_receptor.py")


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V2] STAGE6: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def prepare_receptor(prot_in, prot_out, label):
    """Try meeko first, fall back to obabel."""
    if os.path.exists(MK_PREP_REC):
        proc = subprocess.run([VENV_PY, MK_PREP_REC, "-i", prot_in, "-o", prot_out, "--no-flexible"],
                             capture_output=True, text=True)
        log(f"  meeko receptor {label} rc={proc.returncode}")
        if proc.returncode == 0 and os.path.exists(prot_out) and os.path.getsize(prot_out) > 1000:
            return True
        log(f"  meeko stderr: {proc.stderr[:300]}")
    proc = subprocess.run([OBABEL, prot_in, "-O", prot_out, "-xr"], capture_output=True, text=True)
    log(f"  obabel receptor {label} rc={proc.returncode}")
    return proc.returncode == 0 and os.path.exists(prot_out) and os.path.getsize(prot_out) > 1000


def make_holo_dimer_with_cofactor(prot_dimer_h, cof_a_h, cof_b_h, out_pdb):
    """Concatenate dimer + cofactors into one PDB for holo docking."""
    with open(out_pdb, "w") as out:
        for src in [prot_dimer_h, cof_a_h, cof_b_h]:
            with open(src) as f:
                for line in f:
                    if line.startswith(("ATOM  ", "HETATM", "TER")):
                        out.write(line)
            out.write("TER\n")
        out.write("END\n")
    return out_pdb


def compute_centroid_chainA(pdb_path, residues):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("p", pdb_path)
    coords = []
    for chain in s[0]:
        if chain.id != "A": continue
        for res in chain:
            if res.id[0] == " " and res.id[1] in residues and "CA" in res:
                coords.append(res["CA"].get_coord())
    arr = np.array(coords)
    return arr.mean(axis=0), len(arr)


def parse_vina_result_pdbqt(pdbqt):
    """Parse REMARK VINA RESULT lines from output PDBQT (correct ordering)."""
    affs = []
    if not os.path.exists(pdbqt):
        return affs
    with open(pdbqt) as f:
        for line in f:
            m = re.match(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", line)
            if m:
                affs.append(float(m.group(1)))
    return affs


def parse_pdbqt_models(pdbqt):
    """Yield list of (atom_name, coord) per MODEL."""
    models = []
    cur = []
    with open(pdbqt) as f:
        for line in f:
            if line.startswith("MODEL"):
                cur = []
            elif line.startswith("ENDMDL"):
                models.append(cur)
                cur = []
            elif line.startswith(("ATOM", "HETATM")):
                try:
                    nm = line[12:16].strip()
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    cur.append((nm, np.array([x, y, z])))
                except Exception:
                    pass
    return models


def native_heavy_coords(pdb):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("l", pdb)
    coords = []
    for atom in s.get_atoms():
        if atom.element != "H":
            coords.append((atom.get_name(), atom.get_coord()))
    return coords


def rmsd_to_native(pose_atoms, native_atoms):
    nat_map = {n: c for n, c in native_atoms}
    matched = []
    for n, c in pose_atoms:
        if n in nat_map:
            matched.append((c, nat_map[n]))
    if len(matched) < max(3, len(native_atoms) // 2):
        pose_arr = np.array([c for _, c in pose_atoms])
        nat_arr = np.array([c for _, c in native_atoms])
        if len(pose_arr) == 0 or len(nat_arr) == 0:
            return float("nan")
        used = set()
        pairs = []
        for i in range(len(pose_arr)):
            best = None; best_d = 1e9
            for j in range(len(nat_arr)):
                if j in used: continue
                d = np.linalg.norm(pose_arr[i] - nat_arr[j])
                if d < best_d:
                    best_d = d; best = j
            if best is not None:
                used.add(best)
                pairs.append((pose_arr[i], nat_arr[best]))
        matched = pairs
    if not matched:
        return float("nan")
    sq = [np.sum((a-b)**2) for a, b in matched]
    return float(math.sqrt(sum(sq) / len(sq)))


def split_top_pose_pdbqt(out_pdbqt, top_pdbqt):
    """Extract first MODEL block."""
    with open(out_pdbqt) as f, open(top_pdbqt, "w") as g:
        in_first = False
        for line in f:
            if line.startswith("MODEL 1"):
                in_first = True
            if in_first:
                g.write(line)
                if line.startswith("ENDMDL"):
                    break


def dock_one(rec_pdbqt, lig_pdbqt, out_pdbqt, log_file, centroid, label):
    cx, cy, cz = centroid
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", "18", "--size_y", "18", "--size_z", "18",
           "--exhaustiveness", "32", "--num_modes", "20", "--seed", "42",
           "--out", out_pdbqt, "--cpu", "4"]
    log(f"  vina {label}: cmd={' '.join(cmd[:4])} ... center=({cx:.2f},{cy:.2f},{cz:.2f})")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_file, "w") as f:
        f.write("STDOUT\n" + proc.stdout + "\nSTDERR\n" + proc.stderr)
    return proc


def main():
    os.makedirs(DOCK_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 6 v2 starting")

    selected = json.load(open(os.path.join(AS_DIR, "selected_meta.json")))["selected"]
    prot_h = os.path.join(STR_DIR, "protein_dimer_h.pdb")
    cof_a_h = os.path.join(STR_DIR, "cofactor_chainA_h.pdb")
    cof_b_h = os.path.join(STR_DIR, "cofactor_chainB_h.pdb")
    lig_pdbqt = os.path.join(LIG_DIR, "dump.pdbqt")
    native_lig = os.path.join(STR_DIR, "ligand.pdb")

    centroid, n = compute_centroid_chainA(prot_h, set(selected))
    log(f"chain-A active site centroid from {n} CA: {centroid}")

    # APO receptor: dimer only
    rec_apo_pdbqt = os.path.join(DOCK_DIR, "protein_dimer_apo.pdbqt")
    if not prepare_receptor(prot_h, rec_apo_pdbqt, "APO_dimer"):
        log("FATAL: apo receptor prep failed")
        sys.exit(1)

    # HOLO receptor: dimer + cofactors
    holo_pdb = os.path.join(DOCK_DIR, "protein_dimer_holo.pdb")
    make_holo_dimer_with_cofactor(prot_h, cof_a_h, cof_b_h, holo_pdb)
    rec_holo_pdbqt = os.path.join(DOCK_DIR, "protein_dimer_holo.pdbqt")
    if not prepare_receptor(holo_pdb, rec_holo_pdbqt, "HOLO_dimer"):
        log("WARNING: holo receptor prep failed; falling back to apo for holo run")
        rec_holo_pdbqt = rec_apo_pdbqt

    results = {}
    for cond, rec in [("apo", rec_apo_pdbqt), ("holo", rec_holo_pdbqt)]:
        out_pdbqt = os.path.join(DOCK_DIR, f"wt_{cond}.pdbqt")
        log_file = os.path.join(DOCK_DIR, f"wt_{cond}.log")
        proc = dock_one(rec, lig_pdbqt, out_pdbqt, log_file, centroid, f"WT_{cond}")
        if proc.returncode != 0:
            log(f"  vina WT {cond} rc={proc.returncode}: {proc.stderr[:300]}")
            results[cond] = {"error": proc.stderr[:300]}
            continue
        affs = parse_vina_result_pdbqt(out_pdbqt)
        log(f"  WT {cond} affinities (from PDBQT): {affs[:5]}")
        top = affs[0] if affs else float("nan")
        mean_top3 = float(np.mean(affs[:3])) if len(affs) >= 3 else float("nan")

        poses = parse_pdbqt_models(out_pdbqt)
        native = native_heavy_coords(native_lig)
        rmsd = rmsd_to_native(poses[0], native) if poses else float("nan")

        # Save top pose split
        top_pdbqt = os.path.join(DOCK_DIR, f"wt_{cond}_top.pdbqt")
        split_top_pose_pdbqt(out_pdbqt, top_pdbqt)
        # Convert to PDB for viewer
        top_pdb = os.path.join(DOCK_DIR, f"wt_{cond}_top_pose.pdb")
        subprocess.run([OBABEL, top_pdbqt, "-O", top_pdb], capture_output=True, text=True)
        # Build complex viewer file
        complex_pdb = os.path.join(DOCK_DIR, f"wt_{cond}_complex.pdb")
        with open(complex_pdb, "w") as out:
            with open(prot_h) as f:
                for line in f:
                    if line.startswith(("ATOM  ", "HETATM", "TER")):
                        out.write(line)
            out.write("TER\n")
            with open(top_pdb) as f:
                for line in f:
                    if line.startswith(("ATOM  ", "HETATM")):
                        # rename to ligand chain X
                        new_line = line[:21] + "X" + line[22:]
                        out.write(new_line)
            out.write("TER\nEND\n")

        result = {"top_affinity": top, "mean_top3": mean_top3, "rmsd_top_to_native": rmsd,
                  "n_poses": len(poses), "centroid": centroid.tolist(),
                  "all_affinities": affs,
                  "complex_pdb": complex_pdb, "top_pose_pdb": top_pdb}
        with open(os.path.join(DOCK_DIR, f"wt_{cond}.json"), "w") as f:
            json.dump(result, f, indent=2)
        log(f"  WT {cond}: top={top:.2f}, mean_top3={mean_top3:.2f}, RMSD={rmsd:.2f}")
        results[cond] = result

    # Summary
    with open(os.path.join(DOCK_DIR, "summary.json"), "w") as f:
        json.dump(results, f, indent=2)
    log("Stage 6 v2 DONE")


if __name__ == "__main__":
    main()
