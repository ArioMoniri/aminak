#!/usr/bin/env python3
"""Stage 6 v3: Rebuild WT receptors with non-zero Gasteiger charges and re-dock WT.

Fixes addressed:
- FIX 1: receptor partial charges (try obabel-Gasteiger -> meeko -> pdb2pqr)
- FIX 2: rebuild WT holo reference (centroid on crystal dUMP, exh=96, num_modes=32,
  box 22^3, seeds 42 + 7 sanity check, take best)
- FIX 6: cofactor re-protonated at pH 7.4 (already done in stage3 fallback, redo here)
- FIX 14: echo full Vina parameters into wt_*.json
"""
import os, sys, subprocess, json, re, math, shutil
from datetime import datetime
import numpy as np
from Bio.PDB import PDBParser

PROJECT = os.path.expanduser("~/conserved_site_project")
DOCK_DIR = os.path.join(PROJECT, "06c_docking_wt_v3")
LIG_DIR = os.path.join(PROJECT, "05b_ligand_v2")  # reuse v2 ligand
STR_DIR = os.path.join(PROJECT, "03b_structure_v2")
AS_DIR = os.path.join(PROJECT, "02b_active_site_v2")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v3_06_dock_wt.log")
OBABEL = os.environ.get("OBABEL", "/Users/ario/conserved_site_project/.venv/bin/obabel")
VINA = os.environ.get("VINA", "/opt/homebrew/bin/vina")
VENV_PY = os.path.join(PROJECT, ".venv/bin/python")
MK_PREP_REC = os.path.join(PROJECT, ".venv/bin/mk_prepare_receptor.py")


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V3] STAGE6: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def receptor_max_abs_charge(pdbqt_path):
    if not os.path.exists(pdbqt_path):
        return 0.0
    m = 0.0
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                # PDBQT charge is at columns 67-76 (1-based)
                try:
                    c = float(line[66:76].strip())
                    if abs(c) > m:
                        m = abs(c)
                except (ValueError, IndexError):
                    pass
    return m


def prepare_receptor_with_charges(prot_in, prot_out, label):
    """Try Gasteiger via obabel first, then meeko, then pdb2pqr fallback.
    Returns (ok, method)."""
    # Method 1: obabel + Gasteiger
    proc = subprocess.run(
        [OBABEL, prot_in, "-O", prot_out, "-xr", "-p", "7.4", "--partialcharge", "gasteiger"],
        capture_output=True, text=True
    )
    if proc.returncode == 0 and os.path.exists(prot_out):
        m = receptor_max_abs_charge(prot_out)
        log(f"  [{label}] obabel-gasteiger rc={proc.returncode} max|q|={m:.3f}")
        if m > 0.05:
            return True, "obabel_gasteiger"

    # Method 2: meeko without --no-flexible flag
    if os.path.exists(MK_PREP_REC):
        base = prot_out.rsplit(".", 1)[0]
        proc = subprocess.run([VENV_PY, MK_PREP_REC, "-i", prot_in, "-o", base],
                              capture_output=True, text=True)
        log(f"  [{label}] meeko rc={proc.returncode}")
        if proc.returncode == 0 and os.path.exists(prot_out):
            m = receptor_max_abs_charge(prot_out)
            log(f"  [{label}] meeko max|q|={m:.3f}")
            if m > 0.05:
                return True, "meeko"

    # Method 3: pdb2pqr fallback
    pqr = prot_out.replace(".pdbqt", ".pqr")
    proc = subprocess.run(
        ["pdb2pqr30", "--ff=AMBER", "--with-ph=7.4", prot_in, pqr],
        capture_output=True, text=True
    )
    log(f"  [{label}] pdb2pqr30 rc={proc.returncode}")
    if proc.returncode == 0 and os.path.exists(pqr):
        # Convert to pdbqt with obabel preserving charges from pqr
        proc = subprocess.run([OBABEL, pqr, "-O", prot_out, "-xr"],
                              capture_output=True, text=True)
        m = receptor_max_abs_charge(prot_out)
        log(f"  [{label}] pdb2pqr->pdbqt max|q|={m:.3f}")
        if m > 0.05:
            return True, "pdb2pqr"

    return False, None


def reprotonate_cofactor(in_pdb, out_pdb, label):
    """Re-protonate cofactor at pH 7.4."""
    proc = subprocess.run([OBABEL, in_pdb, "-O", out_pdb, "-h", "-p", "7.4"],
                          capture_output=True, text=True)
    log(f"  [{label}] reprotonate rc={proc.returncode}")
    return proc.returncode == 0


def make_holo_dimer(prot_h, cof_a_h, cof_b_h, out_pdb):
    with open(out_pdb, "w") as out:
        for src in [prot_h, cof_a_h, cof_b_h]:
            with open(src) as f:
                for line in f:
                    if line.startswith(("ATOM  ", "HETATM", "TER")):
                        out.write(line)
            out.write("TER\n")
        out.write("END\n")
    return out_pdb


def crystal_dump_centroid(ligand_pdb):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("l", ligand_pdb)
    coords = []
    for atom in s.get_atoms():
        if atom.element != "H":
            coords.append(atom.get_coord())
    arr = np.array(coords)
    return arr.mean(axis=0)


def parse_vina_pdbqt(pdbqt):
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


def native_heavy(pdb):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("l", pdb)
    coords = []
    for atom in s.get_atoms():
        if atom.element != "H":
            coords.append((atom.get_name(), atom.get_coord()))
    return coords


def rmsd_top(pose_atoms, native_atoms):
    nat_map = {n: c for n, c in native_atoms}
    matched = [(c, nat_map[n]) for n, c in pose_atoms if n in nat_map]
    if len(matched) < max(3, len(native_atoms) // 2):
        # nearest-neighbor matching by index
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
    sq = [np.sum((a - b) ** 2) for a, b in matched]
    return float(math.sqrt(sum(sq) / len(sq)))


def split_top(out_pdbqt, top_pdbqt):
    with open(out_pdbqt) as f, open(top_pdbqt, "w") as g:
        in_first = False
        for line in f:
            if line.startswith("MODEL 1"):
                in_first = True
            if in_first:
                g.write(line)
                if line.startswith("ENDMDL"):
                    break


def dock(rec_pdbqt, lig_pdbqt, out_pdbqt, log_file, centroid, seed, label,
         exh=96, num_modes=32, size=22):
    cx, cy, cz = centroid
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", str(size), "--size_y", str(size), "--size_z", str(size),
           "--exhaustiveness", str(exh), "--num_modes", str(num_modes),
           "--seed", str(seed), "--out", out_pdbqt, "--cpu", "4"]
    log(f"  vina[{label}] seed={seed} exh={exh} center=({cx:.2f},{cy:.2f},{cz:.2f})")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_file, "w") as f:
        f.write("STDOUT\n" + proc.stdout + "\nSTDERR\n" + proc.stderr)
    return proc


def main():
    os.makedirs(DOCK_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 6 v3 starting")

    selected = json.load(open(os.path.join(AS_DIR, "selected_meta.json")))["selected"]
    prot_h = os.path.join(STR_DIR, "protein_dimer_h.pdb")
    cof_a_h = os.path.join(STR_DIR, "cofactor_chainA_h.pdb")
    cof_b_h = os.path.join(STR_DIR, "cofactor_chainB_h.pdb")
    native_lig = os.path.join(STR_DIR, "ligand.pdb")

    # FIX 6: re-protonate cofactor at pH 7.4 (overwrite local v3 copy)
    cof_a_v3 = os.path.join(DOCK_DIR, "cofactor_chainA_h_v3.pdb")
    cof_b_v3 = os.path.join(DOCK_DIR, "cofactor_chainB_h_v3.pdb")
    reprotonate_cofactor(os.path.join(STR_DIR, "cofactor_chainA.pdb"), cof_a_v3, "cofA")
    reprotonate_cofactor(os.path.join(STR_DIR, "cofactor_chainB.pdb"), cof_b_v3, "cofB")

    # FIX 2: centroid on crystal dUMP rather than CA centroid
    centroid = crystal_dump_centroid(native_lig)
    log(f"crystal dUMP centroid: {centroid.tolist()}")

    receptor_methods = {}

    # Apo receptor
    rec_apo = os.path.join(DOCK_DIR, "protein_dimer_apo.pdbqt")
    ok_apo, method_apo = prepare_receptor_with_charges(prot_h, rec_apo, "APO_dimer")
    receptor_methods["apo"] = method_apo
    if not ok_apo:
        log("FATAL: apo receptor prep failed all methods")
        sys.exit(1)
    log(f"  apo receptor OK via {method_apo}, max|q|={receptor_max_abs_charge(rec_apo):.3f}")

    # Holo receptor (with re-protonated cofactor)
    holo_pdb = os.path.join(DOCK_DIR, "protein_dimer_holo.pdb")
    make_holo_dimer(prot_h, cof_a_v3, cof_b_v3, holo_pdb)
    rec_holo = os.path.join(DOCK_DIR, "protein_dimer_holo.pdbqt")
    ok_holo, method_holo = prepare_receptor_with_charges(holo_pdb, rec_holo, "HOLO_dimer")
    receptor_methods["holo"] = method_holo
    if not ok_holo:
        log("WARNING: holo receptor prep failed; falling back to apo for holo run")
        rec_holo = rec_apo
        receptor_methods["holo"] = "fallback_to_apo"
    else:
        log(f"  holo receptor OK via {method_holo}, max|q|={receptor_max_abs_charge(rec_holo):.3f}")

    lig_pdbqt = os.path.join(LIG_DIR, "dump.pdbqt")
    results = {}

    vina_params = {"exhaustiveness": 96, "num_modes": 32, "box_size": 22,
                   "seed_primary": 42, "seed_sanity": 7,
                   "centroid": centroid.tolist(),
                   "centroid_source": "crystal_dump_centroid",
                   "receptor_methods": receptor_methods}

    for cond, rec in [("apo", rec_apo), ("holo", rec_holo)]:
        seed_results = {}
        for seed in [42, 7]:
            tag = f"wt_{cond}_seed{seed}"
            out_pdbqt = os.path.join(DOCK_DIR, f"{tag}.pdbqt")
            log_file = os.path.join(DOCK_DIR, f"{tag}.log")
            proc = dock(rec, lig_pdbqt, out_pdbqt, log_file, centroid, seed,
                        f"WT_{cond}_seed{seed}")
            if proc.returncode != 0:
                log(f"  vina WT {cond} seed{seed} rc={proc.returncode}: {proc.stderr[:200]}")
                seed_results[seed] = None
                continue
            affs = parse_vina_pdbqt(out_pdbqt)
            poses = parse_pdbqt_models(out_pdbqt)
            native = native_heavy(native_lig)
            rmsd = rmsd_top(poses[0], native) if poses else float("nan")
            top = affs[0] if affs else float("nan")
            mean_topk = float(np.mean(affs[:min(3, len(affs))])) if affs else float("nan")
            log(f"  WT {cond} seed{seed}: top={top:.2f} n={len(affs)} rmsd={rmsd:.2f}")
            seed_results[seed] = {"top": top, "mean_topk": mean_topk, "rmsd": rmsd,
                                  "n_modes": len(affs), "all_affinities": affs,
                                  "out_pdbqt": out_pdbqt, "poses": poses, "native": native}

        # Pick best of two seeds: prefer lower RMSD if both have good n_modes; else lower top
        valid = {s: r for s, r in seed_results.items() if r is not None}
        if not valid:
            log(f"  WT {cond}: ALL SEEDS FAILED")
            continue
        # Best = pose with smallest RMSD (closer to crystal)
        best_seed = min(valid.keys(),
                        key=lambda s: (valid[s]["rmsd"] if not math.isnan(valid[s]["rmsd"]) else 1e9,
                                       valid[s]["top"]))
        best = valid[best_seed]
        log(f"  WT {cond} BEST seed = {best_seed} (rmsd={best['rmsd']:.2f}, top={best['top']:.2f})")

        # Promote best to canonical wt_{cond}
        canon_pdbqt = os.path.join(DOCK_DIR, f"wt_{cond}.pdbqt")
        canon_log = os.path.join(DOCK_DIR, f"wt_{cond}.log")
        shutil.copy(best["out_pdbqt"], canon_pdbqt)
        shutil.copy(best["out_pdbqt"].replace(".pdbqt", ".log"), canon_log)

        # Top pose extraction
        top_pdbqt = os.path.join(DOCK_DIR, f"wt_{cond}_top.pdbqt")
        split_top(canon_pdbqt, top_pdbqt)
        top_pdb = os.path.join(DOCK_DIR, f"wt_{cond}_top_pose.pdb")
        subprocess.run([OBABEL, top_pdbqt, "-O", top_pdb], capture_output=True, text=True)

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
                        out.write(line[:21] + "X" + line[22:])
            out.write("TER\nEND\n")

        result = {
            "top_affinity": best["top"],
            "mean_topk": best["mean_topk"],
            "rmsd_top_to_native": best["rmsd"],
            "n_modes": best["n_modes"],
            "n_poses": len(best["poses"]),
            "best_seed": best_seed,
            "all_affinities": best["all_affinities"],
            "all_seed_results": {s: {k: v for k, v in r.items()
                                     if k not in ("poses", "native")}
                                  for s, r in valid.items()},
            "complex_pdb": complex_pdb,
            "top_pose_pdb": top_pdb,
            "vina_params": vina_params,
        }
        with open(os.path.join(DOCK_DIR, f"wt_{cond}.json"), "w") as f:
            json.dump(result, f, indent=2, default=str)
        results[cond] = {k: v for k, v in result.items() if k != "all_seed_results"}

    with open(os.path.join(DOCK_DIR, "summary.json"), "w") as f:
        json.dump({"results": results, "vina_params": vina_params}, f, indent=2, default=str)
    log("Stage 6 v3 DONE")


if __name__ == "__main__":
    main()
