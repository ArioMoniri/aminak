#!/usr/bin/env python3
"""Stage 6 v5: WT holo dock against the v5 in-place reprotonated cofactor
(crystal heavy-atom coords, deprotonated carboxylates, 0 protein clashes).
WT apo is reused verbatim from v4 (apo receptor unchanged; the cofactor
placement bug only affected the holo receptor).
"""
import os, sys, subprocess, json, math, shutil
from datetime import datetime
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common_v5 import (
    parse_vina_pdbqt, parse_pdbqt_models, native_heavy, rmsd_top, split_top,
    restore_atom_names, receptor_max_abs_charge, crystal_dump_centroid,
    make_holo_dimer, prepare_receptor_with_charges,
    OBABEL, VINA, PROJECT,
)

DOCK_V5 = os.path.join(PROJECT, "06e_docking_wt_v5")
DOCK_V4 = os.path.join(PROJECT, "06d_docking_wt_v4")
LIG_DIR = os.path.join(PROJECT, "05b_ligand_v2")
STR_V2 = os.path.join(PROJECT, "03b_structure_v2")
STR_V5 = os.path.join(PROJECT, "03e_structure_v5")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v5_06_dock_wt.log")

VINA_NOISE_FLOOR = 0.85


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V5] STAGE6: {msg}"
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
    log(f"  vina[{label}] seed={seed} exh={exh}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_file, "w") as f:
        f.write("STDOUT\n" + proc.stdout + "\nSTDERR\n" + proc.stderr)
    return proc


def run_seed_sweep(rec, lig_pdbqt, native_lig, centroid, label, dock_dir, seeds,
                   exh=96, num_modes=32):
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
    os.makedirs(DOCK_V5, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 6 v5 starting")

    prot_h = os.path.join(STR_V2, "protein_dimer_h.pdb")
    cof_a_v5 = os.path.join(STR_V5, "cofactor_chainA_v5.pdb")
    cof_b_v5 = os.path.join(STR_V5, "cofactor_chainB_v5.pdb")
    native_lig = os.path.join(STR_V2, "ligand.pdb")
    ref_ligand_pdb = os.path.join(STR_V2, "ligand_h.pdb")

    shutil.copy(cof_a_v5, os.path.join(DOCK_V5, "cofactor_chainA_v5.pdb"))
    shutil.copy(cof_b_v5, os.path.join(DOCK_V5, "cofactor_chainB_v5.pdb"))

    centroid = crystal_dump_centroid(native_lig)
    log(f"  centroid (crystal dUMP): {centroid.tolist()}")

    receptor_methods = {}

    # ---- Apo: copy from v4 verbatim (apo receptor unchanged) ----
    rec_apo = os.path.join(DOCK_V5, "protein_dimer_apo.pdbqt")
    src_apo = os.path.join(DOCK_V4, "protein_dimer_apo.pdbqt")
    if not os.path.exists(src_apo):
        log("FATAL: v4 apo PDBQT missing")
        sys.exit(1)
    shutil.copy(src_apo, rec_apo)
    receptor_methods["apo"] = "reused_from_v4"
    log(f"  apo PDBQT copied from v4 (max|q|={receptor_max_abs_charge(rec_apo):.3f})")

    # ---- Holo: rebuild with v5 cofactor ----
    holo_pdb = os.path.join(DOCK_V5, "protein_dimer_holo.pdb")
    make_holo_dimer(prot_h, cof_a_v5, cof_b_v5, holo_pdb)
    rec_holo = os.path.join(DOCK_V5, "protein_dimer_holo.pdbqt")
    ok_holo, m_holo = prepare_receptor_with_charges(holo_pdb, rec_holo, "HOLO_dimer_v5", log_fn=log)
    receptor_methods["holo"] = m_holo
    if not ok_holo:
        log("FATAL: holo receptor prep failed")
        sys.exit(1)
    log(f"  holo PDBQT OK via {m_holo}, max|q|={receptor_max_abs_charge(rec_holo):.3f}")

    lig_pdbqt = os.path.join(LIG_DIR, "dump.pdbqt")

    primary_seeds = [42, 7, 13, 99, 256]
    fallback_seeds = [1, 2025, 31337]
    actual_exh = {"apo": None, "holo": 96}  # apo reused; holo runs primary at 96

    vina_params_base = {
        "exhaustiveness": 96, "num_modes": 32, "box_size": 22,
        "primary_seeds": primary_seeds, "fallback_seeds": fallback_seeds,
        "centroid": centroid.tolist(),
        "receptor_methods": receptor_methods,
        "selection_rule": "lowest top_affinity, tie-break highest n_modes",
        "vina_noise_floor_kcal_per_mol": VINA_NOISE_FLOOR,
        "actual_exhaustiveness_used": actual_exh,
    }

    results = {}

    # ---- Apo: reuse v4 result ----
    src_apo_json = os.path.join(DOCK_V4, "wt_apo.json")
    src_apo_pdb_top = os.path.join(DOCK_V4, "wt_apo_top_pose.pdb")
    src_apo_complex = os.path.join(DOCK_V4, "wt_apo_complex.pdb")
    src_apo_pdbqt = os.path.join(DOCK_V4, "wt_apo.pdbqt")
    src_apo_log = os.path.join(DOCK_V4, "wt_apo.log")

    apo = json.load(open(src_apo_json))
    # Update embedded paths to v5 for self-containment (still cite v4 source)
    new_apo_pdbqt = os.path.join(DOCK_V5, "wt_apo.pdbqt")
    new_apo_log = os.path.join(DOCK_V5, "wt_apo.log")
    new_apo_top = os.path.join(DOCK_V5, "wt_apo_top_pose.pdb")
    new_apo_cmp = os.path.join(DOCK_V5, "wt_apo_complex.pdb")
    for src, dst in [(src_apo_pdbqt, new_apo_pdbqt), (src_apo_log, new_apo_log),
                      (src_apo_pdb_top, new_apo_top), (src_apo_complex, new_apo_cmp)]:
        if os.path.exists(src):
            shutil.copy(src, dst)
    apo["complex_pdb"] = new_apo_cmp
    apo["top_pose_pdb"] = new_apo_top
    apo["source"] = "reused_from_v4 (apo receptor unchanged)"
    apo["vina_params"] = vina_params_base
    json.dump(apo, open(os.path.join(DOCK_V5, "wt_apo.json"), "w"), indent=2, default=str)
    results["apo"] = {k: v for k, v in apo.items() if k not in ("all_seed_results",)}
    log(f"  WT apo (reused from v4): top={apo['top_affinity']:.2f} "
        f"rmsd={apo['rmsd_top_to_native']:.2f}")

    # ---- Holo: real seed sweep ----
    cond, rec = "holo", rec_holo
    seed_results = run_seed_sweep(rec, lig_pdbqt, native_lig, centroid,
                                   f"wt_{cond}", DOCK_V5, primary_seeds,
                                   exh=96, num_modes=32)
    valid = {s: r for s, r in seed_results.items() if r is not None}
    if not valid:
        log("FATAL: WT holo: ALL primary seeds failed")
        sys.exit(1)

    max_n_modes = max((r["n_modes"] for r in valid.values()), default=0)
    if max_n_modes < 10:
        log(f"  WT holo: max n_modes={max_n_modes} <10, escalating to exh=128")
        extra = run_seed_sweep(rec, lig_pdbqt, native_lig, centroid,
                                f"wt_{cond}", DOCK_V5, fallback_seeds,
                                exh=128, num_modes=32)
        for s, r in extra.items():
            if r is not None:
                valid[s] = r
        actual_exh["holo"] = 128
        log(f"  WT holo: after fallback, max n_modes={max((r['n_modes'] for r in valid.values()), default=0)}")

    best_seed = min(valid.keys(),
                    key=lambda s: (valid[s]["top"] if not math.isnan(valid[s]["top"]) else 1e9,
                                   -valid[s]["n_modes"]))
    best = valid[best_seed]
    log(f"  WT holo BEST seed = {best_seed} top={best['top']:.2f} n_modes={best['n_modes']} "
        f"rmsd={best['rmsd']:.2f}")

    all_top_affs = [r["top"] for r in valid.values() if not math.isnan(r["top"])]
    aff_width = (max(all_top_affs) - min(all_top_affs)) if all_top_affs else float("nan")

    canon_pdbqt = os.path.join(DOCK_V5, f"wt_{cond}.pdbqt")
    canon_log = os.path.join(DOCK_V5, f"wt_{cond}.log")
    shutil.copy(best["out_pdbqt"], canon_pdbqt)
    shutil.copy(best["log_file"], canon_log)

    top_pdbqt = os.path.join(DOCK_V5, f"wt_{cond}_top.pdbqt")
    split_top(canon_pdbqt, top_pdbqt)
    top_pdb = os.path.join(DOCK_V5, f"wt_{cond}_top_pose.pdb")
    restore_atom_names(top_pdbqt, ref_ligand_pdb, top_pdb,
                       reference_dock_pdbqt=lig_pdbqt)

    named_pose = native_heavy(top_pdb)
    crystal = native_heavy(native_lig)
    rmsd_named = rmsd_top(named_pose, crystal)
    log(f"  WT holo named-pose RMSD vs crystal = {rmsd_named:.3f}")

    complex_pdb = os.path.join(DOCK_V5, f"wt_{cond}_complex.pdb")
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

    holo_result = {
        "top_affinity": best["top"],
        "mean_topk": best["mean_topk"],
        "rmsd_top_to_native": best["rmsd"],
        "rmsd_top_named": rmsd_named,
        "n_modes": best["n_modes"],
        "best_seed": best_seed,
        "all_affinities": best["all_affinities"],
        "all_seed_results": {str(s): {k: v for k, v in r.items()
                                      if k not in ("out_pdbqt", "log_file")}
                             for s, r in valid.items()},
        "affinity_distribution_width_kcal": aff_width,
        "complex_pdb": complex_pdb,
        "top_pose_pdb": top_pdb,
        "vina_params": vina_params_base,
    }
    json.dump(holo_result, open(os.path.join(DOCK_V5, "wt_holo.json"), "w"),
              indent=2, default=str)
    results["holo"] = {k: v for k, v in holo_result.items() if k != "all_seed_results"}

    json.dump({"results": results, "vina_params": vina_params_base},
              open(os.path.join(DOCK_V5, "summary.json"), "w"),
              indent=2, default=str)
    log("Stage 6 v5 DONE")


if __name__ == "__main__":
    main()
