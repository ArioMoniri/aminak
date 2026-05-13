#!/usr/bin/env python3
"""Phase 8a: Alternative scoring (Vinardo + AD4) on existing top poses.

For every (mutant, condition) pair in the v5 mutant pipeline plus WT_apo +
WT_holo, re-score the saved Vina top pose with:
  - vina default scoring (sanity check vs the recorded top_affinity)
  - vinardo (Quiroga & Villarreal 2016)
  - ad4 (AutoDock4 force-field) - REQUIRES autogrid4 maps; if unavailable,
    AD4 column is left empty and a note is emitted in the log.

AD4 in Vina 1.2 takes precomputed grid maps via --maps, not a receptor.
Generating those maps requires autogrid4, which is not present on this
Apple Silicon host. The script therefore reports vina + vinardo only,
matching the documented limitation in the README.

Output: 13_phase8/01_alt_scoring/alt_scoring_results.csv with columns
  label, condition, vina_score, vinardo_score, ad4_score,
  delta_vinardo_vs_vina, delta_ad4_vs_vina
"""
from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

PROJECT = Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
OUT = PROJECT / "13_phase8" / "01_alt_scoring"
LOG = PROJECT / "logs" / "v8_phase8.log"
PIPELOG = PROJECT / "pipeline.log"
VINA = "/opt/homebrew/bin/vina"
OBABEL = "/opt/homebrew/bin/obabel"
AUTOGRID = shutil.which("autogrid4")

WT_APO_REC = PROJECT / "06e_docking_wt_v5" / "protein_dimer_apo.pdbqt"
WT_HOLO_REC = PROJECT / "06e_docking_wt_v5" / "protein_dimer_holo.pdbqt"
WT_APO_POSE = PROJECT / "06d_docking_wt_v4" / "wt_apo_top.pdbqt"
WT_HOLO_POSE = PROJECT / "06e_docking_wt_v5" / "wt_holo_top.pdbqt"

# Mutant apo receptor PDBQTs were overwritten by Vina docking output in v3.
# Regenerate them on-the-fly from {mut}_mut_h.pdb (apo source) using obabel
# + Gasteiger charges (same recipe as scripts/v3/stage6_dock_wt_v3.py).
APO_REC_CACHE = PROJECT / "13_phase8" / "01_alt_scoring" / "regen_apo_receptors"


def log(msg: str) -> None:
    line = f"[V8][8a] {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")
    with PIPELOG.open("a") as fh:
        fh.write(line + "\n")


def regen_apo_receptor(label: str) -> Path | None:
    """Regenerate clean mutant apo receptor PDBQT from {mut}_mut_h.pdb using
    obabel + Gasteiger (same recipe as v3 prepare_receptor_with_charges).
    Cached in APO_REC_CACHE.
    """
    APO_REC_CACHE.mkdir(parents=True, exist_ok=True)
    out = APO_REC_CACHE / f"{label}_apo.pdbqt"
    if out.exists() and out.stat().st_size > 0:
        return out
    src = PROJECT / "07e_mut_docking_v5" / label / f"{label}_mut_h.pdb"
    if not src.exists():
        src = PROJECT / "07c_mut_docking_v3" / label / f"{label}_mut_h.pdb"
    if not src.exists():
        return None
    cmd = [OBABEL, str(src), "-O", str(out), "-xr", "-p", "7.4",
           "--partialcharge", "gasteiger"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return None
    if proc.returncode != 0 or not out.exists():
        log(f"obabel apo regen FAIL for {label}: {proc.stderr[:120]}")
        return None
    return out


def receptor_for(label: str, condition: str) -> Path | None:
    """Return the receptor PDBQT used for the original top pose."""
    if label == "WT":
        return WT_APO_REC if condition == "apo" else WT_HOLO_REC
    if condition == "apo":
        # mutant apo receptors were overwritten by Vina output in v3 -> regen
        return regen_apo_receptor(label)
    return PROJECT / "07e_mut_docking_v5" / label / f"{label}_holo.pdbqt"


def pose_for(label: str, condition: str) -> Path:
    if label == "WT":
        return WT_APO_POSE if condition == "apo" else WT_HOLO_POSE
    if condition == "apo":
        return PROJECT / "07c_mut_docking_v3" / label / f"{label}_apo_top.pdbqt"
    return PROJECT / "07e_mut_docking_v5" / label / f"{label}_holo_top.pdbqt"


def strip_model(pose_path: Path, dst_dir: Path) -> Path:
    """Remove MODEL/ENDMDL wrappers - Vina --score_only refuses them."""
    dst_dir.mkdir(parents=True, exist_ok=True)
    out = dst_dir / pose_path.name
    text = pose_path.read_text()
    cleaned_lines = [
        ln for ln in text.splitlines() if not (ln.startswith("MODEL") or ln.startswith("ENDMDL"))
    ]
    out.write_text("\n".join(cleaned_lines) + "\n")
    return out


def run_score(receptor: Path, ligand: Path, scoring: str) -> float | None:
    """Return Estimated Free Energy of Binding in kcal/mol, or None on failure."""
    cmd = [
        VINA,
        "--score_only",
        "--scoring", scoring,
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--autobox",
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    except subprocess.TimeoutExpired:
        return None
    if proc.returncode != 0:
        log(f"vina rc={proc.returncode} scoring={scoring} lig={ligand.name}: {proc.stderr[:120]}")
        return None
    for ln in proc.stdout.splitlines():
        if "Estimated Free Energy of Binding" in ln:
            try:
                # format: "Estimated Free Energy of Binding   : -10.502 (kcal/mol) [...]"
                parts = ln.split(":", 1)[1].strip().split()
                return float(parts[0])
            except (IndexError, ValueError):
                return None
    return None


def collect_targets() -> list[tuple[str, str]]:
    """Read the mutant CSV and return list of (label, condition)."""
    targets: list[tuple[str, str]] = []
    csv_path = PROJECT / "07e_mut_docking_v5" / "mutant_results_v5.csv"
    with csv_path.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("mutant,"):
                continue
            row = line.strip().split(",")
            if len(row) < 3:
                continue
            label, _category, condition = row[0], row[1], row[2]
            targets.append((label, condition))
    return targets


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    tmp_dir = OUT / "stripped_poses"
    log(f"Phase 8a start; AUTOGRID={AUTOGRID}")
    if AUTOGRID is None:
        log("autogrid4 NOT FOUND on host; AD4 column will be empty (documented limitation).")

    targets = collect_targets()
    log(f"Targets: {len(targets)} rows from mutant_results_v5.csv")

    rows = []
    for label, condition in targets:
        rec = receptor_for(label, condition)
        pose = pose_for(label, condition)
        if rec is None or not rec.exists() or not pose.exists():
            log(f"SKIP {label} {condition}: missing receptor={rec is not None and rec.exists()} pose={pose.exists()}")
            rows.append({
                "label": label, "condition": condition,
                "vina_score": "", "vinardo_score": "", "ad4_score": "",
                "delta_vinardo_vs_vina": "", "delta_ad4_vs_vina": "",
                "note": "missing inputs",
            })
            continue

        # Some pose files have MODEL/ENDMDL wrappers (single-model). Strip them.
        clean_pose = strip_model(pose, tmp_dir / f"{label}_{condition}")

        v_vina = run_score(rec, clean_pose, "vina")
        v_vinardo = run_score(rec, clean_pose, "vinardo")
        v_ad4 = None  # AD4 needs precomputed maps - autogrid4 unavailable

        d_vinardo = (v_vinardo - v_vina) if (v_vinardo is not None and v_vina is not None) else None
        d_ad4 = None

        log(f"{label:18s} {condition:4s} vina={v_vina} vinardo={v_vinardo} d_vinardo={d_vinardo}")
        rows.append({
            "label": label,
            "condition": condition,
            "vina_score": "" if v_vina is None else f"{v_vina:.3f}",
            "vinardo_score": "" if v_vinardo is None else f"{v_vinardo:.3f}",
            "ad4_score": "" if v_ad4 is None else f"{v_ad4:.3f}",
            "delta_vinardo_vs_vina": "" if d_vinardo is None else f"{d_vinardo:.3f}",
            "delta_ad4_vs_vina": "" if d_ad4 is None else f"{d_ad4:.3f}",
            "note": "" if (v_vina is not None and v_vinardo is not None) else "score parse fail",
        })

    csv_path = OUT / "alt_scoring_results.csv"
    fields = [
        "label", "condition",
        "vina_score", "vinardo_score", "ad4_score",
        "delta_vinardo_vs_vina", "delta_ad4_vs_vina",
        "note",
    ]
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {csv_path} with {len(rows)} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
