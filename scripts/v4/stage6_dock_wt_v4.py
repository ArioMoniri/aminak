#!/usr/bin/env python3
"""Stage 6 v4: WT apo+holo dock with multi-seed sweep, affinity-based selection,
atom-name preservation. Holo receptor uses v4 reprotonated cofactor.
"""
import os, sys, subprocess, json, math, shutil
from datetime import datetime
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common_v4 import (
    parse_vina_pdbqt, parse_pdbqt_models, native_heavy, rmsd_top, split_top,
    restore_atom_names, receptor_max_abs_charge, crystal_dump_centroid,
    make_holo_dimer, prepare_receptor_with_charges,
    OBABEL, VINA, PROJECT,
)

DOCK_DIR = os.path.join(PROJECT, "06d_docking_wt_v4")
LIG_DIR = os.path.join(PROJECT, "05b_ligand_v2")
STR_V2 = os.path.join(PROJECT, "03b_structure_v2")
STR_V4 = os.path.join(PROJECT, "03d_structure_v4")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v4_06_dock_wt.log")

VINA_NOISE_FLOOR = 0.85  # Trott&Olson 2010, Forli 2016 -- mean of 0.7-1.0


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V4] STAGE6: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


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


def run_seed_sweep(rec, lig_pdbqt, native_lig, centroid, label, dock_dir, seeds,
                   exh=96, num_modes=32):
    """Run vina with multiple seeds; return dict of seed -> result."""
    seed_results = {}
    for seed in seeds:
        tag = f"{label}_seed{seed}"
        out_pdbqt = os.path.join(dock_dir, f"{tag}.pdbqt")
        log_file = os.path.join(dock_dir, f"{tag}.log")
        proc = dock(rec, lig_pdbqt, out_pdbqt, log_file, centroid, seed, label,
                    exh=exh, num_modes=num_modes)
        if proc.returncode != 0:
            log(f"  vina {label} seed{seed} rc={proc.returncode}: {proc.stderr[:200]}")
            seed_results[seed] = None
            continue
        affs = parse_vina_pdbqt(out_pdbqt)
        poses = parse_pdbqt_models(out_pdbqt)
        native = native_heavy(native_lig)
        rmsd = rmsd_top(poses[0], native) if poses else float("nan")
        top = affs[0] if affs else float("nan")
        mean_topk = float(np.mean(affs[:min(3, len(affs))])) if affs else float("nan")
        log(f"  {label} seed{seed}: top={top:.2f} n={len(affs)} rmsd={rmsd:.2f}")
        seed_results[seed] = {
            "top": top, "mean_topk": mean_topk, "rmsd": rmsd,
            "n_modes": len(affs), "all_affinities": affs,
            "out_pdbqt": out_pdbqt, "log_file": log_file,
        }
    return seed_results


def main():
    os.makedirs(DOCK_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 6 v4 starting (multi-seed sweep, affinity-based selection)")

    prot_h = os.path.join(STR_V2, "protein_dimer_h.pdb")
    cof_a_v4 = os.path.join(STR_V4, "cofactor_chainA_v4.pdb")
    cof_b_v4 = os.path.join(STR_V4, "cofactor_chainB_v4.pdb")
    native_lig = os.path.join(STR_V2, "ligand.pdb")  # original named PDB
    ref_ligand_pdb = os.path.join(STR_V2, "ligand_h.pdb")  # for atom name restore

    # Stage out the v4 cofactor copies into DOCK_DIR for reproducibility
    shutil.copy(cof_a_v4, os.path.join(DOCK_DIR, "cofactor_chainA_v4.pdb"))
    shutil.copy(cof_b_v4, os.path.join(DOCK_DIR, "cofactor_chainB_v4.pdb"))

    centroid = crystal_dump_centroid(native_lig)
    log(f"crystal dUMP centroid: {centroid.tolist()}")

    receptor_methods = {}

    # Apo receptor: same as v3 (no cofactor changed) -- but rebuild for v4 lineage
    rec_apo = os.path.join(DOCK_DIR, "protein_dimer_apo.pdbqt")
    ok_apo, method_apo = prepare_receptor_with_charges(prot_h, rec_apo, "APO_dimer", log_fn=log)
    receptor_methods["apo"] = method_apo
    if not ok_apo:
        log("FATAL: apo receptor prep failed")
        sys.exit(1)
    log(f"  apo receptor OK via {method_apo}, max|q|={receptor_max_abs_charge(rec_apo):.3f}")

    # Holo receptor with v4 reprotonated cofactor
    holo_pdb = os.path.join(DOCK_DIR, "protein_dimer_holo.pdb")
    make_holo_dimer(prot_h, cof_a_v4, cof_b_v4, holo_pdb)
    rec_holo = os.path.join(DOCK_DIR, "protein_dimer_holo.pdbqt")
    ok_holo, method_holo = prepare_receptor_with_charges(holo_pdb, rec_holo, "HOLO_dimer", log_fn=log)
    receptor_methods["holo"] = method_holo
    if not ok_holo:
        log("WARNING: holo receptor prep failed; falling back to apo")
        rec_holo = rec_apo
        receptor_methods["holo"] = "fallback_to_apo"
    else:
        log(f"  holo receptor OK via {method_holo}, max|q|={receptor_max_abs_charge(rec_holo):.3f}")

    lig_pdbqt = os.path.join(LIG_DIR, "dump.pdbqt")
    results = {}

    primary_seeds = [42, 7, 13, 99, 256]
    fallback_seeds = [1, 2025, 31337]
    vina_params_base = {
        "exhaustiveness": 96, "num_modes": 32, "box_size": 22,
        "primary_seeds": primary_seeds,
        "fallback_seeds": fallback_seeds,
        "centroid": centroid.tolist(),
        "centroid_source": "crystal_dump_centroid",
        "receptor_methods": receptor_methods,
        "selection_rule": "lowest top_affinity, tie-break by highest n_modes",
        "vina_noise_floor_kcal_per_mol": VINA_NOISE_FLOOR,
        "vina_noise_floor_source": "Trott & Olson 2010; Forli et al. 2016",
    }

    for cond, rec in [("apo", rec_apo), ("holo", rec_holo)]:
        # Primary sweep
        seed_results = run_seed_sweep(rec, lig_pdbqt, native_lig, centroid,
                                       f"wt_{cond}", DOCK_DIR, primary_seeds,
                                       exh=96, num_modes=32)
        valid = {s: r for s, r in seed_results.items() if r is not None}
        if not valid:
            log(f"  WT {cond}: ALL primary seeds failed; trying fallback")
        max_n_modes = max((r["n_modes"] for r in valid.values()), default=0)

        # Fallback if all seeds give <10 modes
        if max_n_modes < 10:
            log(f"  WT {cond}: max n_modes={max_n_modes} <10, escalating to exh=128 with fallback seeds")
            extra = run_seed_sweep(rec, lig_pdbqt, native_lig, centroid,
                                    f"wt_{cond}", DOCK_DIR, fallback_seeds,
                                    exh=128, num_modes=32)
            for s, r in extra.items():
                if r is not None:
                    valid[s] = r
            max_n_modes = max((r["n_modes"] for r in valid.values()), default=0)
            if max_n_modes < 10:
                log(f"  WT {cond}: STILL n_modes<10 after fallback; documented as undersampled")

        if not valid:
            log(f"  WT {cond}: ALL seeds failed")
            continue

        # Best-by-affinity selection (NOT by RMSD); tie-break = highest n_modes
        best_seed = min(valid.keys(),
                        key=lambda s: (valid[s]["top"] if not math.isnan(valid[s]["top"]) else 1e9,
                                       -valid[s]["n_modes"]))
        best = valid[best_seed]
        log(f"  WT {cond} BEST seed = {best_seed} (top={best['top']:.2f}, "
            f"n_modes={best['n_modes']}, rmsd={best['rmsd']:.2f})")

        # Aggregate affinity stats across seeds
        all_top_affs = [r["top"] for r in valid.values() if not math.isnan(r["top"])]
        all_modes_concat = []
        for r in valid.values():
            all_modes_concat.extend(r["all_affinities"])
        affinity_distribution_width = (max(all_top_affs) - min(all_top_affs)) if all_top_affs else float("nan")

        # Promote best to canonical
        canon_pdbqt = os.path.join(DOCK_DIR, f"wt_{cond}.pdbqt")
        canon_log = os.path.join(DOCK_DIR, f"wt_{cond}.log")
        shutil.copy(best["out_pdbqt"], canon_pdbqt)
        shutil.copy(best["log_file"], canon_log)

        # Top pose extract + restore atom names
        top_pdbqt = os.path.join(DOCK_DIR, f"wt_{cond}_top.pdbqt")
        split_top(canon_pdbqt, top_pdbqt)
        top_pdb = os.path.join(DOCK_DIR, f"wt_{cond}_top_pose.pdb")
        restore_atom_names(top_pdbqt, ref_ligand_pdb, top_pdb,
                           reference_dock_pdbqt=lig_pdbqt)

        # Re-compute RMSD on the named pose to verify reproducibility
        from _common_v4 import native_heavy as _nh
        named_pose = _nh(top_pdb)
        crystal = _nh(native_lig)
        rmsd_named = rmsd_top(named_pose, crystal)
        log(f"  WT {cond} named-pose RMSD vs crystal = {rmsd_named:.3f}")

        # Build complex
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
            "rmsd_top_named": rmsd_named,
            "n_modes": best["n_modes"],
            "best_seed": best_seed,
            "selection_rule": "lowest top_affinity, tie-break highest n_modes",
            "all_affinities": best["all_affinities"],
            "all_seed_results": {str(s): {k: v for k, v in r.items()
                                          if k not in ("out_pdbqt", "log_file")}
                                 for s, r in valid.items()},
            "affinity_distribution_width_kcal": affinity_distribution_width,
            "all_modes_concat_n": len(all_modes_concat),
            "complex_pdb": complex_pdb,
            "top_pose_pdb": top_pdb,
            "vina_params": vina_params_base,
        }
        with open(os.path.join(DOCK_DIR, f"wt_{cond}.json"), "w") as f:
            json.dump(result, f, indent=2, default=str)
        results[cond] = {k: v for k, v in result.items() if k != "all_seed_results"}

    with open(os.path.join(DOCK_DIR, "summary.json"), "w") as f:
        json.dump({"results": results, "vina_params": vina_params_base}, f,
                  indent=2, default=str)
    log("Stage 6 v4 DONE")


if __name__ == "__main__":
    main()
