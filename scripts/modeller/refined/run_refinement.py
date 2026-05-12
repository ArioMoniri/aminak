"""Modeller AutoModel refinement at md_level=refine.very_slow.

Re-uses the Phase-6 alignment (10_modeller/03_alignment/alignment.ali)
and template PDBs. Generates 10 refined models with a slower MD
simulated-annealing schedule, which produces better-Ramachandran
geometry than the default 'fast' schedule.

Outputs:
  10b_modeller_refined/02_refined_models/refined_B99990001.pdb ... 10
  10b_modeller_refined/02_refined_models/scores.csv
  10b_modeller_refined/02_refined_models/best_by_dope.json
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
ALIGN_FILE = PROJECT_DIR / "10_modeller/03_alignment/alignment.ali"
KNOWNS_JSON = PROJECT_DIR / "10_modeller/03_alignment/knowns.json"
TEMPLATES_DIR = PROJECT_DIR / "10_modeller/02_blast/templates"
ALIGN_DIR = PROJECT_DIR / "10_modeller/03_alignment"
OUT_DIR = PROJECT_DIR / "10b_modeller_refined/02_refined_models"


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta = json.loads(KNOWNS_JSON.read_text())
    knowns = meta["knowns"]
    sequence_code = meta["sequence_code"]

    # Modeller scratch directory (it dumps lots of intermediate files)
    work_dir = OUT_DIR / "_modeller_scratch"
    work_dir.mkdir(parents=True, exist_ok=True)
    cwd_before = os.getcwd()
    os.chdir(work_dir)
    try:
        from modeller import Environ, log
        from modeller.automodel import AutoModel, assess, refine

        log.minimal()
        env = Environ()
        env.io.atom_files_directory = [
            str(TEMPLATES_DIR.resolve()),
            str(ALIGN_DIR.resolve()),
            ".",
        ]
        env.io.hetatm = False
        env.io.water = False

        a = AutoModel(
            env,
            alnfile=str(ALIGN_FILE.resolve()),
            knowns=knowns,
            sequence=sequence_code,
            assess_methods=(assess.DOPE, assess.GA341),
        )
        # === REFINEMENT KNOBS ===
        a.md_level = refine.very_slow           # slower MD-SA schedule
        a.max_var_iterations = 600              # default 200; allow more optim
        a.repeat_optimization = 2               # extra optim pass
        a.starting_model = 1
        a.ending_model = 10
        print(f"Refinement starting — md_level=refine.very_slow, "
              f"max_var_iterations=600, models=1..10", flush=True)
        a.make()
        print("Refinement complete", flush=True)

        rows = []
        for out in a.outputs:
            if out.get("failure"):
                print(f"FAILED: {out.get('name')}: {out.get('failure')}",
                      flush=True)
                continue
            rows.append({
                "model_pdb": out.get("name"),
                "molpdf": out.get("molpdf"),
                "DOPE": out.get("DOPE score"),
                "GA341": out.get("GA341 score"),
            })
    finally:
        os.chdir(cwd_before)

    if not rows:
        print("No successful models", file=sys.stderr)
        return 3

    # Rename target.B99990001.pdb -> refined_B99990001.pdb and move
    final_rows = []
    for r in rows:
        src = work_dir / r["model_pdb"]
        if not src.exists():
            print(f"WARN: missing {src}", file=sys.stderr)
            continue
        new_name = "refined_" + r["model_pdb"].replace("target.", "")
        dst = OUT_DIR / new_name
        shutil.move(str(src), str(dst))
        r["refined_pdb"] = new_name
        final_rows.append(r)
        print(f"  -> {new_name}  DOPE={r['DOPE']}", flush=True)

    # scores.csv
    with open(OUT_DIR / "scores.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["model_id", "refined_pdb", "molpdf", "DOPE", "GA341"])
        for i, r in enumerate(final_rows, start=1):
            ga = r["GA341"]
            ga_val = ga[0] if isinstance(ga, (list, tuple)) and ga else ga
            w.writerow([i, r["refined_pdb"], r["molpdf"], r["DOPE"], ga_val])

    def dope_val(rr):
        try:
            return float(rr["DOPE"])
        except Exception:
            return float("inf")

    best = min(final_rows, key=dope_val)
    (OUT_DIR / "best_by_dope.json").write_text(
        json.dumps({"refined_pdb": best["refined_pdb"], "DOPE": best["DOPE"]},
                   indent=2))
    # Copy best as a convenience
    shutil.copy(str(OUT_DIR / best["refined_pdb"]),
                str(OUT_DIR / "best_refined.pdb"))
    print(f"Best refined: {best['refined_pdb']}  DOPE={best['DOPE']}",
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
