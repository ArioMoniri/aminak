#!/usr/bin/env python3
"""Stage 7 v5: Re-dock 20 mutants HOLO under the v5 (in-place reprotonated)
cofactor receptor. Apo data are reused from v3 verbatim because the apo
receptor has not changed (the v4 cofactor placement bug only affected holo).
The protocol asymmetry (WT-apo and WT-holo: 5-seed exh=96; mutant-apo: single
seed exh=32 from v3; mutant-holo: single seed exh=32) is documented in
methods (sci-off review item 6).
"""
import os, sys, subprocess, json, math, csv, shutil
from datetime import datetime
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common_v5 import (
    parse_vina_pdbqt, parse_pdbqt_models, native_heavy, rmsd_top, split_top,
    restore_atom_names, receptor_max_abs_charge, crystal_dump_centroid,
    make_holo_dimer, prepare_receptor_with_charges,
    OBABEL, VINA, PROJECT,
)

MUT_V5 = os.path.join(PROJECT, "07e_mut_docking_v5")
VIEWER_DIR = os.path.join(MUT_V5, "viewer_files")
WT_V5 = os.path.join(PROJECT, "06e_docking_wt_v5")
MUT_V3 = os.path.join(PROJECT, "07c_mut_docking_v3")
LIG_DIR = os.path.join(PROJECT, "05b_ligand_v2")
STR_V2 = os.path.join(PROJECT, "03b_structure_v2")
STR_V5 = os.path.join(PROJECT, "03e_structure_v5")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v5_07_mutants.log")


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V5] STAGE7: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def dock(rec_pdbqt, lig_pdbqt, out_pdbqt, log_file, centroid):
    cx, cy, cz = centroid
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", "22", "--size_y", "22", "--size_z", "22",
           "--exhaustiveness", "32", "--num_modes", "20", "--seed", "42",
           "--out", out_pdbqt, "--cpu", "4"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_file, "w") as f:
        f.write("STDOUT\n" + proc.stdout + "\nSTDERR\n" + proc.stderr)
    return proc


def build_complex(receptor_pdb, top_pose_pdb, out_complex):
    with open(out_complex, "w") as out:
        with open(receptor_pdb) as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM", "TER")):
                    out.write(line)
        out.write("TER\n")
        with open(top_pose_pdb) as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM")):
                    out.write(line[:21] + "X" + line[22:])
        out.write("TER\nEND\n")


def main():
    os.makedirs(MUT_V5, exist_ok=True)
    os.makedirs(VIEWER_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 7 v5 starting")

    mutants = sorted([d for d in os.listdir(MUT_V3)
                      if os.path.isdir(os.path.join(MUT_V3, d)) and d != "viewer_files"])
    log(f"Found {len(mutants)} mutants from v3")

    cof_a_v5 = os.path.join(STR_V5, "cofactor_chainA_v5.pdb")
    cof_b_v5 = os.path.join(STR_V5, "cofactor_chainB_v5.pdb")
    native_lig = os.path.join(STR_V2, "ligand.pdb")
    ref_ligand_pdb = os.path.join(STR_V2, "ligand_h.pdb")
    centroid = crystal_dump_centroid(native_lig)
    log(f"  centroid: {centroid.tolist()}")

    lig_pdbqt = os.path.join(LIG_DIR, "dump.pdbqt")
    wt_apo = json.load(open(os.path.join(WT_V5, "wt_apo.json")))
    wt_holo = json.load(open(os.path.join(WT_V5, "wt_holo.json")))
    wt_apo_aff = wt_apo["top_affinity"]
    wt_holo_aff = wt_holo["top_affinity"]
    log(f"  WT v5 apo aff = {wt_apo_aff:.2f}, WT v5 holo aff = {wt_holo_aff:.2f}")

    # ---- Inherit category info from v3 ----
    v3_csv = os.path.join(MUT_V3, "mutant_results_v3.csv")
    cat_map = {}
    apo_v3_map = {}
    with open(v3_csv) as f:
        f.readline()
        reader = csv.DictReader(f)
        for row in reader:
            mid = row["mutant"]
            cat_map[mid] = row["category"]
            if row["condition"] == "apo":
                apo_v3_map[mid] = row

    rows = []
    skipped = []
    for mi, mid in enumerate(mutants):
        cat = cat_map.get(mid, "unknown")
        log(f"[{mi+1}/{len(mutants)}] {mid} ({cat})")

        sub = os.path.join(MUT_V5, mid)
        os.makedirs(sub, exist_ok=True)

        v3_mut_h = os.path.join(MUT_V3, mid, f"{mid}_mut_h.pdb")
        if not os.path.exists(v3_mut_h):
            log(f"  v3 mut_h missing for {mid}; skipping")
            skipped.append((mid, "no_v3_mut_h"))
            continue
        v5_mut_h = os.path.join(sub, f"{mid}_mut_h.pdb")
        shutil.copy(v3_mut_h, v5_mut_h)

        # Build holo receptor with v5 cofactor
        holo_pdb = os.path.join(sub, f"{mid}_holo.pdb")
        make_holo_dimer(v5_mut_h, cof_a_v5, cof_b_v5, holo_pdb)
        rec_holo = os.path.join(sub, f"{mid}_holo.pdbqt")
        ok_holo, m_holo = prepare_receptor_with_charges(holo_pdb, rec_holo,
                                                        f"{mid}_holo", log_fn=log)
        if not ok_holo:
            log(f"  holo receptor prep failed; skipping mutant")
            skipped.append((mid, "holo_prep_failed"))
            continue

        # Dock
        out_pdbqt = os.path.join(sub, f"{mid}_holo_dock.pdbqt")
        log_file = os.path.join(sub, f"{mid}_holo.log")
        proc = dock(rec_holo, lig_pdbqt, out_pdbqt, log_file, centroid)
        if proc.returncode != 0:
            log(f"  vina {mid} holo rc={proc.returncode}: {proc.stderr[:200]}")
            rows.append({"mutant": mid, "category": cat, "condition": "holo",
                         "top_affinity": float("nan"), "mean_topk": float("nan"),
                         "n_modes": 0, "rmsd_to_native": float("nan"),
                         "delta_vina_vs_wt": float("nan"),
                         "mis_docked": True, "low_confidence": True,
                         "n_clashes": 0, "error": proc.stderr[:200]})
        else:
            affs = parse_vina_pdbqt(out_pdbqt)
            n_modes = len(affs)
            top = affs[0] if affs else float("nan")
            mean_topk = float(np.mean(affs[:min(3, n_modes)])) if n_modes > 0 else float("nan")

            poses = parse_pdbqt_models(out_pdbqt)
            native = native_heavy(native_lig)
            rmsd = rmsd_top(poses[0], native) if poses else float("nan")
            mis_docked = (not math.isnan(rmsd)) and rmsd > 3.0
            low_confidence = n_modes < 5

            delta = (top - wt_holo_aff) if not math.isnan(top) else float("nan")

            top_pdbqt = os.path.join(sub, f"{mid}_holo_top.pdbqt")
            split_top(out_pdbqt, top_pdbqt)
            top_pose_pdb = os.path.join(VIEWER_DIR, f"{mid}_holo_top_pose.pdb")
            restore_atom_names(top_pdbqt, ref_ligand_pdb, top_pose_pdb,
                               reference_dock_pdbqt=lig_pdbqt)
            complex_pdb = os.path.join(VIEWER_DIR, f"{mid}_holo_complex.pdb")
            build_complex(holo_pdb, top_pose_pdb, complex_pdb)

            rows.append({"mutant": mid, "category": cat, "condition": "holo",
                         "top_affinity": top, "mean_topk": mean_topk,
                         "n_modes": n_modes,
                         "rmsd_to_native": rmsd,
                         "delta_vina_vs_wt": delta,
                         "mis_docked": mis_docked,
                         "low_confidence": low_confidence,
                         "n_clashes": 0,
                         "complex_pdb": complex_pdb, "top_pose_pdb": top_pose_pdb})
            log(f"  {mid} holo: top={top:.2f} delta={delta:+.2f} n={n_modes} "
                f"rmsd={rmsd:.2f} mis_docked={mis_docked} low_conf={low_confidence}")

        # Reuse v3 apo
        if mid in apo_v3_map:
            v3a = apo_v3_map[mid]
            try:
                top_apo = float(v3a["top_affinity"])
            except (ValueError, TypeError):
                top_apo = float("nan")
            try:
                mean_topk_apo = float(v3a["mean_topk"])
            except (ValueError, TypeError):
                mean_topk_apo = float("nan")
            try:
                n_modes_apo = int(v3a["n_modes"])
            except (ValueError, TypeError):
                n_modes_apo = 0
            try:
                rmsd_apo = float(v3a["rmsd_to_native"])
            except (ValueError, TypeError):
                rmsd_apo = float("nan")
            mis_docked_apo = v3a.get("mis_docked", "False").lower() == "true"
            low_confidence_apo = n_modes_apo < 5
            delta_apo = (top_apo - wt_apo_aff) if not math.isnan(top_apo) else float("nan")

            v3_complex = v3a.get("complex_pdb", "")
            v3_top_pose = v3a.get("top_pose_pdb", "")
            v5_top = os.path.join(VIEWER_DIR, f"{mid}_apo_top_pose.pdb")
            v5_cmp = os.path.join(VIEWER_DIR, f"{mid}_apo_complex.pdb")
            if os.path.exists(v3_top_pose) and not os.path.exists(v5_top):
                shutil.copy(v3_top_pose, v5_top)
            if os.path.exists(v3_complex) and not os.path.exists(v5_cmp):
                shutil.copy(v3_complex, v5_cmp)

            rows.append({"mutant": mid, "category": cat, "condition": "apo",
                         "top_affinity": top_apo, "mean_topk": mean_topk_apo,
                         "n_modes": n_modes_apo,
                         "rmsd_to_native": rmsd_apo,
                         "delta_vina_vs_wt": delta_apo,
                         "mis_docked": mis_docked_apo,
                         "low_confidence": low_confidence_apo,
                         "n_clashes": int(v3a.get("n_clashes") or 0),
                         "complex_pdb": v5_cmp, "top_pose_pdb": v5_top,
                         "apo_source": "reused_from_v3"})
            log(f"  {mid} apo (reused v3): top={top_apo:.2f} delta={delta_apo:+.2f} "
                f"n={n_modes_apo} rmsd={rmsd_apo:.2f}")

    fieldnames = ["mutant", "category", "condition", "top_affinity", "mean_topk",
                  "n_modes", "rmsd_to_native", "delta_vina_vs_wt",
                  "mis_docked", "low_confidence", "n_clashes",
                  "complex_pdb", "top_pose_pdb", "apo_source", "error"]
    csv_path = os.path.join(MUT_V5, "mutant_results_v5.csv")
    with open(csv_path, "w", newline="") as f:
        f.write("# Sign convention: delta_vina_vs_wt = top_aff_mut - top_aff_wt_v5; positive = destabilising.\n")
        f.write("# mean_topk = mean(affinities[:min(3, n_modes)]).\n")
        f.write("# low_confidence = (n_modes < 5). Excluded from rankings.\n")
        f.write("# Apo dockings reused from v3 (apo receptor unchanged across v3->v4->v5).\n")
        f.write("# WT apo reused from v4 multi-seed sweep; WT holo re-run in v5 with the in-place reprotonated cofactor.\n")
        wt_apo_row = {"mutant": "WT", "category": "wildtype", "condition": "apo",
                      "top_affinity": wt_apo_aff,
                      "mean_topk": wt_apo.get("mean_topk", float("nan")),
                      "n_modes": wt_apo.get("n_modes", 0),
                      "rmsd_to_native": wt_apo.get("rmsd_top_to_native", float("nan")),
                      "delta_vina_vs_wt": 0.000,
                      "mis_docked": False,
                      "low_confidence": (wt_apo.get("n_modes", 0) < 5),
                      "n_clashes": 0,
                      "complex_pdb": wt_apo.get("complex_pdb", ""),
                      "top_pose_pdb": wt_apo.get("top_pose_pdb", "")}
        wt_holo_row = {"mutant": "WT", "category": "wildtype", "condition": "holo",
                       "top_affinity": wt_holo_aff,
                       "mean_topk": wt_holo.get("mean_topk", float("nan")),
                       "n_modes": wt_holo.get("n_modes", 0),
                       "rmsd_to_native": wt_holo.get("rmsd_top_to_native", float("nan")),
                       "delta_vina_vs_wt": 0.000,
                       "mis_docked": False,
                       "low_confidence": (wt_holo.get("n_modes", 0) < 5),
                       "n_clashes": 0,
                       "complex_pdb": wt_holo.get("complex_pdb", ""),
                       "top_pose_pdb": wt_holo.get("top_pose_pdb", "")}
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerow(wt_apo_row)
        w.writerow(wt_holo_row)
        for r in rows:
            w.writerow(r)
    log(f"wrote {csv_path} ({len(rows)+2} rows including WT reference)")

    json.dump({
        "n_mutants": len(mutants),
        "n_rows": len(rows),
        "skipped": skipped,
        "wt_apo_v5": wt_apo_aff,
        "wt_holo_v5": wt_holo_aff,
        "wt_holo_v5_n_modes": wt_holo.get("n_modes", 0),
        "wt_holo_v5_best_seed": wt_holo.get("best_seed"),
        "wt_holo_v5_rmsd": wt_holo.get("rmsd_top_to_native"),
        "sign_convention": "delta_vina_vs_wt = top_aff_mut - top_aff_wt_v5; positive = destabilising",
        "low_confidence_rule": "n_modes < 5 (holo only); excluded from rankings",
        "apo_source": "reused from v3 (receptor unchanged)",
        "vina_noise_floor_kcal_per_mol": 0.85,
    }, open(os.path.join(MUT_V5, "summary_v5.json"), "w"), indent=2, default=str)
    log(f"Stage 7 v5 DONE: {len(rows)} rows, {len(skipped)} skipped")


if __name__ == "__main__":
    main()
