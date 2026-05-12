"""Collect DOPE/molpdf/seq-id from refined PDB REMARK headers.

After a Modeller AutoModel run, the REMARK 6 lines in each model PDB
carry the MODELLER OBJECTIVE FUNCTION (molpdf) and the % SEQ ID.
DOPE is normally added by assess.DOPE; if absent we recompute it here
by calling assess.DOPE on each finished PDB.
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
from pathlib import Path

PROJECT_DIR = Path(os.environ.get("PROJECT_DIR",
                                  str(Path.home() / "conserved_site_project")))
OUT_DIR = PROJECT_DIR / "10b_modeller_refined/02_refined_models"
ALIGN_DIR = PROJECT_DIR / "10_modeller/03_alignment"
TEMPLATES_DIR = PROJECT_DIR / "10_modeller/02_blast/templates"


def scrape(pdb: Path) -> dict:
    out = {"model_pdb": pdb.name, "molpdf": None, "seq_id": None, "DOPE": None,
           "GA341": None}
    with open(pdb) as fh:
        for line in fh:
            if not line.startswith("REMARK"):
                if line.startswith("ATOM"):
                    break
                continue
            m = re.search(r"OBJECTIVE FUNCTION:\s*([-\d.]+)", line)
            if m:
                out["molpdf"] = float(m.group(1))
            m = re.search(r"SEQ ID:\s*([-\d.]+)", line)
            if m:
                out["seq_id"] = float(m.group(1))
            m = re.search(r"DOPE SCORE:\s*([-\d.]+)", line)
            if m:
                out["DOPE"] = float(m.group(1))
            m = re.search(r"GA341 SCORE:\s*([-\d.]+)", line)
            if m:
                out["GA341"] = float(m.group(1))
    return out


def recompute_dope(pdb: Path) -> float:
    from modeller import Environ, Selection
    from modeller.scripts import complete_pdb
    env = Environ()
    env.io.atom_files_directory = [str(pdb.parent), "."]
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")
    mdl = complete_pdb(env, str(pdb))
    return Selection(mdl).assess_dope()


def main() -> int:
    pdbs = sorted(OUT_DIR.glob("refined_B999*.pdb"))
    if not pdbs:
        print("No refined_B999 PDBs found", file=sys.stderr)
        return 2
    rows = []
    for p in pdbs:
        info = scrape(p)
        # If DOPE missing, compute it
        if info["DOPE"] is None:
            try:
                info["DOPE"] = recompute_dope(p)
            except Exception as e:
                print(f"WARN: DOPE recompute failed for {p.name}: {e}",
                      file=sys.stderr)
        rows.append(info)
        print(f"{p.name}  molpdf={info['molpdf']}  DOPE={info['DOPE']}  "
              f"GA341={info['GA341']}", flush=True)

    csv_path = OUT_DIR / "scores.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["model_id", "refined_pdb", "molpdf", "DOPE", "GA341",
                    "seq_id"])
        for i, r in enumerate(rows, start=1):
            w.writerow([i, r["model_pdb"], r["molpdf"], r["DOPE"],
                        r["GA341"], r["seq_id"]])
    print(f"Wrote {csv_path}", flush=True)

    # Best by DOPE
    def dope_val(rr):
        try:
            return float(rr["DOPE"])
        except Exception:
            return float("inf")
    best = min(rows, key=dope_val)
    (OUT_DIR / "best_by_dope.json").write_text(
        json.dumps({"refined_pdb": best["model_pdb"], "DOPE": best["DOPE"]},
                   indent=2))
    import shutil
    shutil.copy(str(OUT_DIR / best["model_pdb"]),
                str(OUT_DIR / "best_refined.pdb"))
    print(f"Best: {best['model_pdb']} DOPE={best['DOPE']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
