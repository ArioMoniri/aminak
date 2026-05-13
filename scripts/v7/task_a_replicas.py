#!/usr/bin/env python3
"""Task A: Multi-replica Vina docking.

Re-dock WT (apo + holo) and 8 key mutants (holo) 5 times each with
seeds [42, 7, 13, 99, 256] to quantify run-to-run stochastic spread
of the top affinity.

All runs use the canonical box centroid from the v5 mutant pipeline
(-0.137, 4.232, 15.159), size 18 A, exhaustiveness 32, num_modes 20.
"""
from __future__ import annotations
import csv
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

PROJECT = Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
OUT = PROJECT / "12_phase7" / "01_replicas"
LOG = PROJECT / "logs" / "v7_task_a.log"
PIPELOG = PROJECT / "pipeline.log"
VINA = "/opt/homebrew/bin/vina"
LIGAND = PROJECT / "05b_ligand_v2" / "dump.pdbqt"

CENTER = (-0.137, 4.232, 15.159)
SIZE = (18.0, 18.0, 18.0)
EXH = 32
NMODES = 20
SEEDS = [42, 7, 13, 99, 256]

WT_APO = PROJECT / "06e_docking_wt_v5" / "protein_dimer_apo.pdbqt"
WT_HOLO = PROJECT / "06e_docking_wt_v5" / "protein_dimer_holo.pdbqt"

MUTANTS = [
    "R215A_N226A",
    "H196A",
    "R215E",
    "R50A",
    "C195A",
    "R175E_R176E",
    "T170A",
    "Y258F_F225Y",
]

def log(msg: str) -> None:
    line = f"[V7][taskA] {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
    print(line, flush=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")
    with PIPELOG.open("a") as fh:
        fh.write(line + "\n")

def run_vina(receptor: Path, label: str, seed: int) -> dict:
    out_pdbqt = OUT / "raw" / f"{label}_seed{seed}.pdbqt"
    out_log = OUT / "raw" / f"{label}_seed{seed}.log"
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        VINA,
        "--receptor", str(receptor),
        "--ligand", str(LIGAND),
        "--center_x", f"{CENTER[0]}",
        "--center_y", f"{CENTER[1]}",
        "--center_z", f"{CENTER[2]}",
        "--size_x", f"{SIZE[0]}",
        "--size_y", f"{SIZE[1]}",
        "--size_z", f"{SIZE[2]}",
        "--exhaustiveness", str(EXH),
        "--num_modes", str(NMODES),
        "--seed", str(seed),
        "--out", str(out_pdbqt),
    ]
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    dt = time.time() - t0
    out_log.write_text(proc.stdout + "\n--STDERR--\n" + proc.stderr)
    if proc.returncode != 0:
        log(f"FAIL {label} seed {seed}: rc={proc.returncode}")
        return {"top_affinity": None, "n_modes": 0, "elapsed_s": dt, "error": proc.stderr[:200]}
    # Parse top affinity from stdout
    top = None
    n = 0
    for line in proc.stdout.splitlines():
        s = line.strip()
        if not s or not s[0].isdigit():
            continue
        parts = s.split()
        if len(parts) >= 2:
            try:
                rank = int(parts[0])
                aff = float(parts[1])
                if rank == 1:
                    top = aff
                n = max(n, rank)
            except ValueError:
                continue
    return {"top_affinity": top, "n_modes": n, "elapsed_s": round(dt, 2), "error": None}

def receptor_for(label: str) -> Path:
    if label == "WT_apo":
        return WT_APO
    if label == "WT_holo":
        return WT_HOLO
    # mutant
    mut = label.replace("_holo", "")
    return PROJECT / "07e_mut_docking_v5" / mut / f"{mut}_holo.pdbqt"

def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    log("Task A start")
    targets = ["WT_apo", "WT_holo"] + [f"{m}_holo" for m in MUTANTS]
    rows = []
    summary = {"params": {"center": CENTER, "size": SIZE, "exhaustiveness": EXH,
                          "num_modes": NMODES, "seeds": SEEDS},
               "results": {}}
    for label in targets:
        rec = receptor_for(label)
        if not rec.exists():
            log(f"SKIP {label}: receptor missing {rec}")
            continue
        affs = []
        per_seed = {}
        for seed in SEEDS:
            log(f"docking {label} seed={seed}")
            r = run_vina(rec, label, seed)
            per_seed[str(seed)] = r
            if r["top_affinity"] is not None:
                affs.append(r["top_affinity"])
            rows.append({
                "label": label,
                "seed": seed,
                "top_affinity": r["top_affinity"],
                "n_modes": r["n_modes"],
                "elapsed_s": r["elapsed_s"],
                "error": r["error"] or "",
            })
        if affs:
            mean = statistics.mean(affs)
            sd = statistics.stdev(affs) if len(affs) > 1 else 0.0
            mn = min(affs)
            mx = max(affs)
        else:
            mean = sd = mn = mx = None
        summary["results"][label] = {
            "per_seed": per_seed,
            "mean_top_affinity": mean,
            "sd_top_affinity": sd,
            "min_top_affinity": mn,
            "max_top_affinity": mx,
            "spread_kcal_mol": (mx - mn) if mn is not None else None,
            "n_seeds_ok": len(affs),
        }
        log(f"{label}: mean={mean} sd={sd} spread={(mx-mn) if mn is not None else None}")
    # write CSV
    csv_path = OUT / "multi_replica_results.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["label","seed","top_affinity","n_modes","elapsed_s","error"])
        w.writeheader()
        w.writerows(rows)
    # write per-target aggregate CSV
    agg_csv = OUT / "multi_replica_aggregate.csv"
    with agg_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["label","mean_top_affinity","sd_top_affinity","min_top_affinity",
                    "max_top_affinity","spread_kcal_mol","n_seeds_ok",
                    "seed42","seed7","seed13","seed99","seed256"])
        for label, agg in summary["results"].items():
            ps = agg["per_seed"]
            w.writerow([label, agg["mean_top_affinity"], agg["sd_top_affinity"],
                        agg["min_top_affinity"], agg["max_top_affinity"],
                        agg["spread_kcal_mol"], agg["n_seeds_ok"],
                        ps.get("42",{}).get("top_affinity"),
                        ps.get("7",{}).get("top_affinity"),
                        ps.get("13",{}).get("top_affinity"),
                        ps.get("99",{}).get("top_affinity"),
                        ps.get("256",{}).get("top_affinity")])
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    log(f"Task A done: csv={csv_path} agg={agg_csv}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
