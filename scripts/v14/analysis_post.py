#!/usr/bin/env python3
"""Phase 14 — post-docking analysis: SASA, pose-cluster diversity, water-bridge check (S1).

Run after each strategy's docking completes. Reads results_raw.csv, adds analysis columns,
writes results_analysed.csv and per-pose JSON artefacts.

Usage:  python3 analysis_post.py <strategy_dir>
        e.g.  python3 analysis_post.py 14_inhibitor_design/01_active_site
"""
from __future__ import annotations
import sys, csv, json, glob
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import REPO, RECEPTOR_APO

def parse_pose_pdbqt_heavy(path: Path):
    coords = []
    in_m1 = False
    with path.open() as fh:
        for ln in fh:
            if ln.startswith("MODEL 1"): in_m1 = True; continue
            if ln.startswith("ENDMDL") and in_m1: break
            if in_m1 and (ln.startswith("ATOM") or ln.startswith("HETATM")):
                el = ln[76:78].strip() if len(ln) >= 78 else ""
                # AD type col 78+ — strip "H" types
                if el.startswith("H"): continue
                try:
                    x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
                    coords.append((x, y, z))
                except ValueError: continue
    return np.array(coords) if coords else None


def parse_pose_all_modes_heavy(path: Path):
    """Return dict {model_idx: ndarray(N,3)}."""
    out = {}
    cur_model = None; cur_coords = []
    with path.open() as fh:
        for ln in fh:
            if ln.startswith("MODEL "):
                if cur_model is not None and cur_coords: out[cur_model] = np.array(cur_coords)
                cur_model = int(ln.split()[1]); cur_coords = []
            elif ln.startswith("ENDMDL"):
                if cur_model is not None and cur_coords: out[cur_model] = np.array(cur_coords)
                cur_model = None; cur_coords = []
            elif cur_model is not None and (ln.startswith("ATOM") or ln.startswith("HETATM")):
                el = ln[76:78].strip() if len(ln) >= 78 else ""
                if el.startswith("H"): continue
                try:
                    x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
                    cur_coords.append((x, y, z))
                except ValueError: continue
    return out


def water_bridge_check(pose_pdbqt: Path, xtal_1hvy: Path) -> dict:
    """E1b — strict crystal water-bridge presence check.
    For each crystal water within 3.5 Å of Tyr258:OH AND within 3.5 Å of any pose heavy atom: flag."""
    import MDAnalysis as mda
    if not xtal_1hvy.exists():
        return {"water_bridge_lost": None, "reason": "no_crystal_pdb"}
    pose_heavy = parse_pose_pdbqt_heavy(pose_pdbqt)
    if pose_heavy is None or len(pose_heavy) == 0:
        return {"water_bridge_lost": None, "reason": "no_pose_coords"}
    u = mda.Universe(str(xtal_1hvy))
    waters_o = u.select_atoms("resname HOH and name O")
    # Tyr258 may be in chain A (segid 'A' or none, depending on PDB)
    tyr258 = u.select_atoms("resid 258 and name OH")
    if len(tyr258) == 0:
        return {"water_bridge_lost": None, "reason": "no_tyr258_in_crystal"}
    tyr_pos = tyr258[0].position
    flagged = []
    for w in waters_o:
        d_tyr = float(np.linalg.norm(w.position - tyr_pos))
        if d_tyr > 3.5: continue
        d_pose = float(np.min(np.linalg.norm(pose_heavy - w.position, axis=1)))
        if d_pose <= 3.5:
            flagged.append({"resid": int(w.resid), "d_tyr": round(d_tyr,3), "d_pose": round(d_pose,3)})
    return {"water_bridge_lost": len(flagged) > 0, "water_bridge_residues": flagged}


def pose_cluster_count(pose_pdbqt: Path, eps: float = 2.0) -> dict:
    """DBSCAN on heavy-atom RMSD across all modes in a single pose file."""
    from sklearn.cluster import DBSCAN
    modes = parse_pose_all_modes_heavy(pose_pdbqt)
    if len(modes) < 2:
        return {"n_modes": len(modes), "n_clusters": len(modes)}
    keys = sorted(modes.keys())
    n = len(keys)
    if min(modes[k].shape[0] for k in keys) == 0:
        return {"n_modes": n, "n_clusters": n}
    # min atom count across modes
    m = min(modes[k].shape[0] for k in keys)
    D = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1, n):
            d = float(np.sqrt(np.mean(np.sum((modes[keys[i]][:m] - modes[keys[j]][:m])**2, axis=1))))
            D[i,j] = D[j,i] = d
    labels = DBSCAN(eps=eps, min_samples=2, metric="precomputed").fit_predict(D)
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    return {"n_modes": n, "n_clusters": n_clusters, "labels": labels.tolist()}


def ligand_sasa_buried_pct(pose_pdbqt: Path, receptor_pdbqt: Path) -> dict:
    """Buried-SASA % = 1 - (complex_lig_SASA / free_lig_SASA)."""
    import freesasa
    # write pose MODEL 1 as standalone PDB
    pose_only_pdb = pose_pdbqt.with_suffix(".pose1.pdb")
    if not pose_only_pdb.exists():
        lines = []; in_m1 = False
        for ln in pose_pdbqt.read_text().splitlines():
            if ln.startswith("MODEL 1"): in_m1 = True; continue
            if ln.startswith("ENDMDL") and in_m1: break
            if in_m1 and (ln.startswith("ATOM") or ln.startswith("HETATM")):
                lines.append(ln[:66])
        pose_only_pdb.write_text("\n".join(lines) + "\nEND\n")
    # receptor as PDB
    rec_pdb = receptor_pdbqt.with_suffix(".pdb")
    if not rec_pdb.exists():
        lines = []
        for ln in receptor_pdbqt.read_text().splitlines():
            if ln.startswith(("ATOM","HETATM")): lines.append(ln[:66])
            elif ln.startswith(("END","TER")): lines.append(ln[:66])
        rec_pdb.write_text("\n".join(lines) + "\nEND\n")
    # combined PDB
    combined = pose_pdbqt.with_suffix(".complex.pdb")
    combined.write_text(rec_pdb.read_text().rstrip("END\n").rstrip() + "\n" + pose_only_pdb.read_text())
    try:
        s_free   = freesasa.calc(freesasa.Structure(str(pose_only_pdb), freesasa.Classifier()))
        s_complex = freesasa.calc(freesasa.Structure(str(combined), freesasa.Classifier()))
        # SASA of the ligand atoms in the complex: take last N atoms of complex structure
        # Simpler: total ligand-only SASA minus total complex SASA from receptor-only baseline
        s_recv   = freesasa.calc(freesasa.Structure(str(rec_pdb), freesasa.Classifier()))
        # buried_lig ≈ free_lig + receptor - complex
        free_lig = float(s_free.totalArea())
        complex_tot = float(s_complex.totalArea())
        rec_alone = float(s_recv.totalArea())
        # change due to ligand-receptor contact: free_lig + rec - complex
        buried_total = (free_lig + rec_alone) - complex_tot
        buried_pct = buried_total / free_lig if free_lig > 0 else None
        return {"free_lig_sasa": free_lig, "complex_total_sasa": complex_tot,
                "receptor_alone_sasa": rec_alone, "buried_sasa": buried_total,
                "buried_pct": buried_pct}
    except Exception as e:
        return {"error": str(e)[:200]}


def main(strategy_dir: str):
    sd = Path(strategy_dir)
    raw_csv = sd / "results_raw.csv"
    if not raw_csv.exists():
        print(f"!!! {raw_csv} not found"); sys.exit(1)
    xtal_1hvy = REPO / "03_structure" / "1hvy.pdb"
    rows = list(csv.DictReader(raw_csv.open()))
    print(f"Analysing {len(rows)} dock rows in {sd.name} …")
    do_waterbridge = ("active_site" in sd.name)
    out_rows = []
    for i, r in enumerate(rows):
        # pose file path differs by strategy schema (compound or peptide or fragment column)
        name = r.get("compound") or r.get("peptide") or r.get("fragment")
        state = r.get("state","") or ""
        seed = r.get("seed","")
        cav = r.get("cavity_id","")
        if cav:
            pose = sd / "docked" / f"cav{cav}_{name}_seed{seed}.pdbqt"
        elif state:
            pose = sd / "docked" / f"{name}_{state}_seed{seed}.pdbqt"
        else:
            pose = sd / "docked" / f"{name}_seed{seed}.pdbqt"
        ar = dict(r)
        if pose.exists() and pose.stat().st_size > 0:
            # cluster
            try: pc = pose_cluster_count(pose)
            except Exception as e: pc = {"error": str(e)[:100]}
            ar["n_pose_clusters"] = pc.get("n_clusters")
            # SASA only for top1 per compound × state (skip multi-seed dup)
            if seed in ("42",):
                try:
                    sa = ligand_sasa_buried_pct(pose, RECEPTOR_APO)
                    ar["sasa_buried_pct"] = sa.get("buried_pct")
                except Exception as e:
                    ar["sasa_buried_pct"] = None; ar["sasa_err"] = str(e)[:100]
            # water bridge
            if do_waterbridge and seed in ("42",):
                try:
                    wb = water_bridge_check(pose, xtal_1hvy)
                    ar["water_bridge_lost"] = wb.get("water_bridge_lost")
                    ar["water_bridge_residues"] = json.dumps(wb.get("water_bridge_residues", []))
                except Exception as e:
                    ar["water_bridge_lost"] = None; ar["wb_err"] = str(e)[:100]
        out_rows.append(ar)
        if (i+1) % 5 == 0: print(f"  [{i+1}/{len(rows)}] {name} {state} s={seed}")
    out_csv = sd / "results_analysed.csv"
    cols = list(out_rows[0].keys())
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(out_rows)
    print(f"  → {out_csv}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "14_inhibitor_design/01_active_site")
