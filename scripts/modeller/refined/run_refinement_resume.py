"""Resume Modeller refinement for models 9-10 only.

The initial run completed models 1-8 then the parent shell process died.
This script runs models 9-10 with identical settings into the same scratch dir.
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

START = int(os.environ.get("START_MODEL", "9"))
END = int(os.environ.get("END_MODEL", "10"))


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta = json.loads(KNOWNS_JSON.read_text())
    knowns = meta["knowns"]
    sequence_code = meta["sequence_code"]

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
        a.md_level = refine.very_slow
        a.max_var_iterations = 600
        a.repeat_optimization = 2
        a.starting_model = START
        a.ending_model = END
        print(f"Resume: refinement models {START}..{END}, "
              f"md_level=refine.very_slow", flush=True)
        a.make()
        print("Resume refinement complete", flush=True)

        rows = []
        for out in a.outputs:
            if out.get("failure"):
                print(f"FAIL: {out.get('name')}: {out.get('failure')}",
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

    for r in rows:
        src = work_dir / r["model_pdb"]
        if not src.exists():
            print(f"WARN missing {src}", file=sys.stderr)
            continue
        new_name = "refined_" + r["model_pdb"].replace("target.", "")
        dst = OUT_DIR / new_name
        shutil.move(str(src), str(dst))
        r["refined_pdb"] = new_name
        print(f"  -> {new_name}  DOPE={r['DOPE']}", flush=True)

    # Append to scores.csv if exists, else create
    csv_path = OUT_DIR / "scores.csv"
    exists = csv_path.exists()
    with open(csv_path, "a", newline="") as fh:
        w = csv.writer(fh)
        if not exists:
            w.writerow(["model_id", "refined_pdb", "molpdf", "DOPE", "GA341"])
        for i, r in enumerate(rows, start=START):
            ga = r["GA341"]
            ga_val = ga[0] if isinstance(ga, (list, tuple)) and ga else ga
            w.writerow([i, r.get("refined_pdb", r["model_pdb"]),
                        r["molpdf"], r["DOPE"], ga_val])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
