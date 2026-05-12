"""Step 4 — Run Modeller AutoModel to generate 10 models for the target.

Produces target.B9999000{1..10}.pdb (or target.B9999000N.pdb for N>=10).
Parses Modeller log for DOPE / molpdf / GA341 scores and writes scores.csv.
Models are moved to 04_modeller_run/models/.
"""
from __future__ import annotations

import csv
import json
import os
import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    MODELS_DIR, STEP3_DIR, STEP4_DIR, TEMPLATES_DIR, setup_logger,
)

LOG = setup_logger("step4_run_modeller")


def main() -> int:
    STEP4_DIR.mkdir(parents=True, exist_ok=True)
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

    knowns_meta = json.loads((STEP3_DIR / "knowns.json").read_text())
    knowns: list[str] = knowns_meta["knowns"]
    sequence_code: str = knowns_meta["sequence_code"]

    # Modeller needs the alignment file and the template PDBs accessible.
    # We chdir into STEP4_DIR and point atom_files_directory to absolute paths.
    work_dir = STEP4_DIR
    cwd_before = os.getcwd()
    os.chdir(work_dir)
    try:
        # Modeller imports
        from modeller import Environ, log
        from modeller.automodel import AutoModel, assess

        log.minimal()
        env = Environ()
        env.io.atom_files_directory = [
            str(TEMPLATES_DIR.resolve()),
            str(STEP3_DIR.resolve()),
            ".",
        ]
        # Allow heteroatoms / waters off by default; we want clean models
        env.io.hetatm = False
        env.io.water = False

        a = AutoModel(
            env,
            alnfile=str((STEP3_DIR / "alignment.ali").resolve()),
            knowns=knowns,
            sequence=sequence_code,
            assess_methods=(assess.DOPE, assess.GA341),
        )
        a.starting_model = 1
        a.ending_model = 10
        LOG.info("Starting AutoModel.make()  knowns=%s  models=1..10", knowns)
        a.make()
        LOG.info("AutoModel.make() done")

        # `a.outputs` contains a list of dicts per model with 'name', 'failure',
        # 'molpdf', 'DOPE score', 'GA341 score'
        rows = []
        for out in a.outputs:
            if out.get("failure"):
                LOG.warning("Model failed: %s -> %s", out.get("name"), out.get("failure"))
                continue
            row = {
                "model_pdb": out.get("name"),
                "molpdf": out.get("molpdf"),
                "DOPE": out.get("DOPE score"),
                "GA341": out.get("GA341 score"),
            }
            rows.append(row)
            LOG.info("Model %s  molpdf=%s  DOPE=%s  GA341=%s",
                     row["model_pdb"], row["molpdf"], row["DOPE"], row["GA341"])
    finally:
        os.chdir(cwd_before)

    if not rows:
        LOG.error("No successful models — aborting")
        return 4

    # Move models into MODELS_DIR
    moved = []
    for row in rows:
        src = STEP4_DIR / row["model_pdb"]
        dst = MODELS_DIR / row["model_pdb"]
        if src.exists():
            shutil.move(str(src), str(dst))
            moved.append(dst.name)
        else:
            LOG.warning("Expected model file missing: %s", src)
    LOG.info("Moved %d models to %s", len(moved), MODELS_DIR)

    # Write scores.csv
    csv_path = STEP4_DIR / "scores.csv"
    # GA341 might be a list — flatten
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["model_id", "model_pdb", "molpdf", "DOPE", "GA341"])
        for i, row in enumerate(rows, start=1):
            ga = row["GA341"]
            ga_val = ga[0] if isinstance(ga, (list, tuple)) and ga else ga
            w.writerow([i, row["model_pdb"], row["molpdf"], row["DOPE"], ga_val])
    LOG.info("Wrote %s", csv_path)

    # Identify best by DOPE (lowest, more negative = better)
    def dope_val(r):
        try:
            return float(r["DOPE"])
        except Exception:
            return float("inf")

    best = min(rows, key=dope_val)
    (STEP4_DIR / "best_by_dope.json").write_text(
        json.dumps({"model_pdb": best["model_pdb"], "DOPE": best["DOPE"]}, indent=2)
    )
    # Also create a best_model.pdb symlink/copy for the validation step's convenience.
    best_dst = MODELS_DIR / "best_model.pdb"
    src_best = MODELS_DIR / best["model_pdb"]
    if src_best.exists():
        shutil.copy(str(src_best), str(best_dst))
        LOG.info("Best (by DOPE) -> %s (copied to best_model.pdb)", best["model_pdb"])

    LOG.info("STEP 4 OK: %d models", len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
