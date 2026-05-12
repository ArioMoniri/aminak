"""Modeller LoopModel refinement of the 93-101 divergent loop.

Uses the best-by-DOPE refined model from 02_refined_models/ as input.
This sub-samples the residue 93-101 loop only (where templates were
uninformative — the gap in 6K7Q_A), with 10 sub-models per template.

Outputs:
  10b_modeller_refined/03_loop_refined/loop_<id>.BL00000001.pdb ... etc.
  10b_modeller_refined/03_loop_refined/best_loop_refined.pdb
  10b_modeller_refined/03_loop_refined/loop_scores.csv
"""
from __future__ import annotations

import csv
import json
import os
import shutil
import sys
from pathlib import Path

PROJECT_DIR = Path(os.environ.get("PROJECT_DIR",
                                  str(Path.home() / "conserved_site_project")))
REFINED_DIR = PROJECT_DIR / "10b_modeller_refined/02_refined_models"
OUT_DIR = PROJECT_DIR / "10b_modeller_refined/03_loop_refined"

LOOP_START = 93
LOOP_END = 101


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Pick best refined model
    best_meta = json.loads((REFINED_DIR / "best_by_dope.json").read_text())
    best_pdb = REFINED_DIR / best_meta["refined_pdb"]
    print(f"Input model: {best_pdb.name}  (DOPE={best_meta['DOPE']})",
          flush=True)

    # Modeller needs a working dir + the input PDB visible there
    work_dir = OUT_DIR / "_modeller_scratch"
    work_dir.mkdir(parents=True, exist_ok=True)
    # Copy the input pdb in (Modeller looks by file basename in atom_files_directory)
    input_basename = "loop_input.pdb"
    shutil.copy(str(best_pdb), str(work_dir / input_basename))

    cwd_before = os.getcwd()
    os.chdir(work_dir)
    try:
        from modeller import Environ, Selection, log
        from modeller.automodel import LoopModel, assess, refine

        class MyLoop(LoopModel):
            def select_loop_atoms(self):
                # Inclusive residue range on chain A
                return Selection(self.residue_range(f"{LOOP_START}:A",
                                                    f"{LOOP_END}:A"))

        log.minimal()
        env = Environ()
        env.io.atom_files_directory = ["."]
        env.io.hetatm = False
        env.io.water = False

        m = MyLoop(
            env,
            inimodel=input_basename,
            sequence="loop_target",
            loop_assess_methods=(assess.DOPE,),
        )
        m.loop.starting_model = 1
        m.loop.ending_model = 10
        m.loop.md_level = refine.very_slow
        m.loop.max_var_iterations = 600
        print(f"Loop refinement starting — residues {LOOP_START}-{LOOP_END}, "
              f"10 sub-models, md_level=refine.very_slow", flush=True)
        m.make()
        print("Loop refinement complete", flush=True)

        rows = []
        for out in m.loop.outputs:
            if out.get("failure"):
                print(f"FAIL: {out.get('name')}: {out.get('failure')}",
                      flush=True)
                continue
            rows.append({
                "loop_pdb": out.get("name"),
                "DOPE": out.get("DOPE score"),
                "molpdf": out.get("molpdf"),
            })
    finally:
        os.chdir(cwd_before)

    if not rows:
        print("No successful loop sub-models", file=sys.stderr)
        return 4

    # Move all loop sub-models out
    moved = []
    for r in rows:
        src = work_dir / r["loop_pdb"]
        if not src.exists():
            print(f"WARN: missing {src}", file=sys.stderr)
            continue
        dst = OUT_DIR / r["loop_pdb"]
        shutil.move(str(src), str(dst))
        r["final_pdb"] = dst.name
        moved.append(r)
        print(f"  -> {dst.name}  DOPE={r['DOPE']}", flush=True)

    with open(OUT_DIR / "loop_scores.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["model_id", "loop_pdb", "molpdf", "DOPE"])
        for i, r in enumerate(moved, start=1):
            w.writerow([i, r["final_pdb"], r["molpdf"], r["DOPE"]])

    def dope_val(rr):
        try:
            return float(rr["DOPE"])
        except Exception:
            return float("inf")

    best = min(moved, key=dope_val)
    shutil.copy(str(OUT_DIR / best["final_pdb"]),
                str(OUT_DIR / "best_loop_refined.pdb"))
    (OUT_DIR / "best_loop_meta.json").write_text(
        json.dumps({"loop_pdb": best["final_pdb"], "DOPE": best["DOPE"],
                    "input_model": best_pdb.name,
                    "loop_residues": [LOOP_START, LOOP_END]}, indent=2))
    print(f"Best loop refinement: {best['final_pdb']}  DOPE={best['DOPE']}",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
