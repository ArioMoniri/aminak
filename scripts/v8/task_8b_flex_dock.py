#!/usr/bin/env python3
"""Phase 8b: Flexible-residue Vina re-dock on the 8 priority mutants.

For each priority mutant:
  1. Regenerate a clean apo PDBQT from {mut}_mut_h.pdb using OpenBabel +
     Gasteiger charges (same recipe as v3 prepare_receptor_with_charges),
     then post-process to restore source residue numbering (obabel renumbers
     1..288; source PDB starts at residue 26, so we add +25 to all resnums).
  2. Append the cofactor PDBQTs (chain A and chain B) as receptor ATOM
     records to make a clean holo PDBQT with chains A/B intact and resnums
     matching the source PDB.
  3. Use scripts/v8/flex_split.py to split into rigid + flex PDBQTs for the
     14 active-site panel residues on chain A (the chain closest to the
     binding box centroid).
  4. Dock with Vina (--exhaustiveness 32, --num_modes 20, --seed 42, same
     box as the v5 mutant pipeline).
  5. Save outputs under 13_phase8/02_flexres/.
  6. Compare flex top affinity to the recorded rigid Vina top affinity from
     07e_mut_docking_v5/mutant_results_v5.csv.

Output:
  13_phase8/02_flexres/<mutant>_flex.pdbqt       Vina docking result
  13_phase8/02_flexres/<mutant>_flex.log         Vina stdout
  13_phase8/02_flexres/<mutant>_rigid.pdbqt      receptor (rigid)
  13_phase8/02_flexres/<mutant>_flexres.pdbqt    flex residues
  13_phase8/02_flexres/flexres_compare.csv       summary table
"""
from __future__ import annotations

import csv
import os
import subprocess
import sys
import time
from pathlib import Path

PROJECT = Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
OUT = PROJECT / "13_phase8" / "02_flexres"
LOG = PROJECT / "logs" / "v8_phase8.log"
PIPELOG = PROJECT / "pipeline.log"

VINA = "/opt/homebrew/bin/vina"
OBABEL = "/opt/homebrew/bin/obabel"
LIGAND = PROJECT / "05b_ligand_v2" / "dump.pdbqt"

CENTER = (-0.137, 4.232, 15.159)
SIZE = (18.0, 18.0, 18.0)
# Flex docking adds ~25-30 rotatable DOF (8 panel residues x 2-4 chi each)
# on top of the ligand's own torsions.  Vina's default exh=8 can take >1h
# per run with this many flex DOF; we reduce to exh=8 (Vina's default) and
# num_modes=10 to keep per-mutant wall time around 5-10 min.
EXH = 8
NMODES = 10
SEED = 42
RESNUM_OFFSET = 25  # obabel renumbers 1..N; source PDB starts at 26

COFACTOR_A = PROJECT / "06f_receptor_fixed" / "cofactor_A.pdbqt"
COFACTOR_B = PROJECT / "06f_receptor_fixed" / "cofactor_B.pdbqt"

# 8 priority mutants (same as scripts/v7/task_a_replicas.py)
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

# Flex residues = the active-site panel on chain A (chain nearest to the
# binding box centroid in this dimer).  We use the 8 highest-conservation
# binding residues that have *rotatable* side chains -- GLY (no side
# chain) and PRO (ring-locked) are excluded.  These 8 residues account
# for ~22 chi DOF, which combined with the ligand's torsions is a
# tractable search at exh=8.
FLEX_PANEL: list[tuple[str, int]] = [
    ("A", 50),    # R50  - dUMP binding
    ("A", 109),   # W109 - dUMP/D16 binding
    ("A", 175),   # R175 - dUMP binding
    ("A", 176),   # R176 - dUMP binding
    ("A", 195),   # C195 - catalytic nucleophile (ALA in C195A mutant)
    ("A", 196),   # H196 - dUMP binding
    ("A", 214),   # Q214 - UMP binding
    ("A", 215),   # R215 - dUMP binding
]

# Ensure flex_split is importable
sys.path.insert(0, str(PROJECT / "scripts" / "v8"))
from flex_split import split_clean_pdbqt  # noqa: E402


def log(msg: str) -> None:
    line = f"[V8][8b] {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")
    with PIPELOG.open("a") as fh:
        fh.write(line + "\n")


def regen_apo_pdbqt(mut: str, out_path: Path) -> bool:
    """Generate clean apo PDBQT (chains A/B preserved, resnums shifted to
    match source PDB)."""
    src = PROJECT / "07e_mut_docking_v5" / mut / f"{mut}_mut_h.pdb"
    if not src.exists():
        src = PROJECT / "07c_mut_docking_v3" / mut / f"{mut}_mut_h.pdb"
    if not src.exists():
        log(f"  no mut_h.pdb for {mut}")
        return False
    tmp = out_path.parent / f".{out_path.stem}.raw.pdbqt"
    cmd = [OBABEL, str(src), "-O", str(tmp), "-xr", "-p", "7.4",
           "--partialcharge", "gasteiger"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        return False
    if proc.returncode != 0 or not tmp.exists():
        log(f"  obabel rc={proc.returncode}: {proc.stderr[:200]}")
        return False

    # Post-process: shift residue numbers by +RESNUM_OFFSET.
    out_lines = []
    for ln in tmp.read_text().splitlines():
        if ln.startswith("ATOM"):
            try:
                n = int(ln[22:26])
                new_n = n + RESNUM_OFFSET
                ln = ln[:22] + f"{new_n:>4d}" + ln[26:]
            except ValueError:
                pass
        out_lines.append(ln)
    out_path.write_text("\n".join(out_lines) + "\n")
    tmp.unlink(missing_ok=True)
    return True


def build_holo_pdbqt(apo_path: Path, out_path: Path) -> bool:
    """Append cofactor A + cofactor B ATOM records to apo to make holo."""
    if not (COFACTOR_A.exists() and COFACTOR_B.exists()):
        return False
    apo_lines = apo_path.read_text().splitlines()
    cof_a = [ln for ln in COFACTOR_A.read_text().splitlines() if ln.startswith("ATOM")]
    cof_b = [ln for ln in COFACTOR_B.read_text().splitlines() if ln.startswith("ATOM")]

    out = []
    maxser = 0
    for ln in apo_lines:
        if ln.startswith(("REMARK", "ATOM")):
            out.append(ln)
            if ln.startswith("ATOM"):
                try:
                    s = int(ln[6:11])
                except ValueError:
                    s = 0
                if s > maxser:
                    maxser = s
    out.append("TER")
    s = maxser + 1
    for ln in cof_a + cof_b:
        out.append(ln[:6] + f"{s:5d}" + ln[11:])
        s += 1
    out.append("TER")
    out.append("END")
    out_path.write_text("\n".join(out) + "\n")
    return True


def run_vina_flex(rigid: Path, flex: Path, out_pdbqt: Path, out_log: Path) -> dict:
    cmd = [
        VINA,
        "--receptor", str(rigid),
        "--flex", str(flex),
        "--ligand", str(LIGAND),
        "--center_x", f"{CENTER[0]}",
        "--center_y", f"{CENTER[1]}",
        "--center_z", f"{CENTER[2]}",
        "--size_x", f"{SIZE[0]}",
        "--size_y", f"{SIZE[1]}",
        "--size_z", f"{SIZE[2]}",
        "--exhaustiveness", str(EXH),
        "--num_modes", str(NMODES),
        "--seed", str(SEED),
        "--out", str(out_pdbqt),
    ]
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    except subprocess.TimeoutExpired:
        return {"top_affinity": None, "n_modes": 0, "elapsed_s": 900, "error": "timeout"}
    dt = time.time() - t0
    out_log.write_text(proc.stdout + "\n--STDERR--\n" + proc.stderr)
    if proc.returncode != 0:
        return {"top_affinity": None, "n_modes": 0, "elapsed_s": dt, "error": proc.stderr[:200]}
    # Parse top affinity from "mode | affinity ..." table
    top = None
    n_modes = 0
    for ln in proc.stdout.splitlines():
        s = ln.strip()
        if not s or not s[0].isdigit():
            continue
        parts = s.split()
        if len(parts) >= 2:
            try:
                rank = int(parts[0])
                aff = float(parts[1])
                if rank == 1:
                    top = aff
                if rank > n_modes:
                    n_modes = rank
            except ValueError:
                continue
    return {"top_affinity": top, "n_modes": n_modes, "elapsed_s": round(dt, 2),
            "error": None}


def load_rigid_top_aff(mut: str) -> float | None:
    csv_path = PROJECT / "07e_mut_docking_v5" / "mutant_results_v5.csv"
    with csv_path.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("mutant,"):
                continue
            row = line.strip().split(",")
            if len(row) < 4:
                continue
            if row[0] == mut and row[2] == "holo":
                try:
                    return float(row[3])
                except ValueError:
                    return None
    return None


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    log(f"Phase 8b start; flex panel = {FLEX_PANEL}")

    rows = []
    for mut in MUTANTS:
        log(f"=== {mut} ===")
        apo_path = OUT / f"{mut}_apo_clean.pdbqt"
        holo_path = OUT / f"{mut}_holo_clean.pdbqt"
        rigid_path = OUT / f"{mut}_rigid.pdbqt"
        flex_path = OUT / f"{mut}_flexres.pdbqt"
        out_dock = OUT / f"{mut}_flex.pdbqt"
        out_log = OUT / f"{mut}_flex.log"

        if not regen_apo_pdbqt(mut, apo_path):
            log(f"  apo regen FAIL")
            rows.append({"label": mut, "rigid_vina_score": "",
                         "flex_vina_score": "", "delta_flex": "",
                         "n_flex_modes": 0, "error": "apo regen fail"})
            continue
        if not build_holo_pdbqt(apo_path, holo_path):
            log(f"  holo build FAIL")
            rows.append({"label": mut, "rigid_vina_score": "",
                         "flex_vina_score": "", "delta_flex": "",
                         "n_flex_modes": 0, "error": "holo build fail"})
            continue

        rigid_str, flex_str, warns = split_clean_pdbqt(
            holo_path.read_text(), FLEX_PANEL,
        )
        for w in warns:
            log(f"  flex_split WARN: {w}")
        rigid_path.write_text(rigid_str)
        flex_path.write_text(flex_str)
        if not flex_str.strip():
            log(f"  no flex blocks generated; skipping")
            rows.append({"label": mut, "rigid_vina_score": "",
                         "flex_vina_score": "", "delta_flex": "",
                         "n_flex_modes": 0, "error": "no flex blocks"})
            continue

        log(f"  docking with --flex (exhaustiveness={EXH}, modes={NMODES})")
        r = run_vina_flex(rigid_path, flex_path, out_dock, out_log)
        log(f"  flex top_aff={r['top_affinity']} n_modes={r['n_modes']} elapsed={r['elapsed_s']}s")

        rigid_aff = load_rigid_top_aff(mut)
        delta = (r["top_affinity"] - rigid_aff) if (r["top_affinity"] is not None and rigid_aff is not None) else None

        rows.append({
            "label": mut,
            "rigid_vina_score": "" if rigid_aff is None else f"{rigid_aff:.3f}",
            "flex_vina_score": "" if r["top_affinity"] is None else f"{r['top_affinity']:.3f}",
            "delta_flex": "" if delta is None else f"{delta:.3f}",
            "n_flex_modes": r["n_modes"],
            "elapsed_s": r["elapsed_s"],
            "error": r["error"] or "",
        })

    csv_path = OUT / "flexres_compare.csv"
    fields = ["label", "rigid_vina_score", "flex_vina_score", "delta_flex",
              "n_flex_modes", "elapsed_s", "error"]
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {csv_path} with {len(rows)} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
