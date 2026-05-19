#!/usr/bin/env python3
"""Phase 14g — HADDOCK3-score on the LR-octapeptide vs TYMS-dimer interface.

Full HADDOCK3 needs CNS which isn't installed. But haddock3-score is a
standalone scoring tool that:
  - reads any PDB complex
  - computes the HADDOCK score (vdW + electrostatic + AIR + desolvation)
  - does NOT need CNS

This lets us at least score the existing Vina LR-octapeptide poses with
HADDOCK's empirical PPI-scoring function, and compare against the scrambled
control. If canonical ranks better than scrambled by HADDOCK score (even
without full flexref refinement), that is a partial answer to the Strategy-3
question.
"""
from __future__ import annotations
import subprocess, csv, re, shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "14_inhibitor_design" / "07_advanced_methods" / "haddock3_run"
OUT.mkdir(parents=True, exist_ok=True)
HADDOCK3_SCORE = Path.home() / "Library/Python/3.14/bin/haddock3-score"

# Phase 14 Strategy 3 docked peptide poses
S3_DIR = REPO / "14_inhibitor_design" / "03_dimer_interface" / "docked"
APO_RECEPTOR_PDB = REPO / "06f_receptor_fixed" / "dimer_noH.pdb"

# Map: peptide label -> docked PDBQT (MODEL 1 only needed)
PEPTIDES = [
    ("LR8_canonical",     S3_DIR / "LR8_LSCQLYQR_seed42.pdbqt",      "canonical Cardinale-2011 octapeptide"),
    ("LR8_scrambled",     S3_DIR / "LR8_scrambled_QLCRQSYL_seed42.pdbqt", "shuffled-seq specificity control"),
    ("LR4_pos1",          S3_DIR / "LR_4mer_pos1_LSCQ_seed42.pdbqt", "overlapping 4-mer LSCQ"),
    ("LR4_pos3",          S3_DIR / "LR_4mer_pos3_CQLY_seed42.pdbqt", "overlapping 4-mer CQLY"),
    ("LR4_pos5",          S3_DIR / "LR_4mer_pos5_LYQR_seed42.pdbqt", "overlapping 4-mer LYQR"),
]


def pdbqt_pose_to_pdb(pose_pdbqt: Path, out_pdb: Path):
    """Extract MODEL 1 of a Vina pose PDBQT to a clean PDB."""
    lines = []; in_m1 = False
    for ln in pose_pdbqt.read_text().splitlines():
        if ln.startswith("MODEL 1"): in_m1 = True; continue
        if ln.startswith("ENDMDL") and in_m1: break
        if in_m1 and (ln.startswith("ATOM") or ln.startswith("HETATM")):
            # rewrite as ATOM, chain Z, resname LIG
            ln_clean = ln[:66]
            # HADDOCK needs each peptide residue to look like a real protein residue —
            # but for haddock3-score, simple ATOM lines work
            lines.append(ln_clean)
    out_pdb.write_text("\n".join(lines) + "\nEND\n")


def build_complex_pdb(receptor_pdb: Path, peptide_pdb: Path, out_pdb: Path):
    """Concatenate receptor + peptide as a single complex for haddock3-score."""
    body = []
    for ln in receptor_pdb.read_text().splitlines():
        if ln.startswith(("ATOM","HETATM","TER")):
            body.append(ln)
    body.append("TER")
    for ln in peptide_pdb.read_text().splitlines():
        if ln.startswith(("ATOM","HETATM")):
            body.append(ln)
    body.append("END")
    out_pdb.write_text("\n".join(body) + "\n")


def haddock3_score(complex_pdb: Path) -> dict:
    """Run haddock3-score on a complex PDB. Returns the parsed score."""
    cmd = [str(HADDOCK3_SCORE), str(complex_pdb), "--run_dir",
           str(complex_pdb.parent / "_run"), "--full"]
    try:
        proc = subprocess.run(cmd, check=False, capture_output=True, timeout=300,
                              cwd=str(complex_pdb.parent))
    except subprocess.TimeoutExpired:
        return {"ok": False, "err": "timeout"}
    out = proc.stdout.decode() + proc.stderr.decode()
    # Look for "HADDOCK-score = X" or similar
    m = re.search(r"HADDOCK\s*score\s*[=:]\s*([-\d.]+)", out, re.I)
    if not m:
        # alternate: "score" on its own line
        for ln in out.splitlines():
            ln = ln.strip()
            if ln.startswith("score") or "haddock-score" in ln.lower():
                m2 = re.search(r"([-\d]+\.\d+)", ln)
                if m2: m = m2; break
    energies = {}
    for term, regex in [("vdw", r"v\s?d\s?w\s*[=:]\s*([-\d.]+)"),
                        ("elec", r"e\s?l\s?e\s?c(?:trostatic)?\s*[=:]\s*([-\d.]+)"),
                        ("air", r"air\s*[=:]\s*([-\d.]+)"),
                        ("desolv", r"desolv\s*[=:]\s*([-\d.]+)"),
                        ("bsa", r"bsa\s*[=:]\s*([-\d.]+)")]:
        em = re.search(regex, out, re.I)
        if em:
            try: energies[term] = float(em.group(1))
            except ValueError: pass
    return {
        "ok": True,
        "haddock_score": float(m.group(1)) if m else None,
        "energies": energies,
        "stdout_tail": out[-400:] if out else "",
    }


def main():
    print("=== Phase 14g — HADDOCK3-score on Strategy-3 peptide complexes ===")
    if not HADDOCK3_SCORE.exists():
        print(f"  ! haddock3-score binary not found at {HADDOCK3_SCORE}")
        return

    rows = []
    for label, pose_pdbqt, comment in PEPTIDES:
        print(f"  --- {label}  ({comment}) ---")
        if not pose_pdbqt.exists():
            print(f"    ! missing pose: {pose_pdbqt}"); continue
        wdir = OUT / label
        wdir.mkdir(parents=True, exist_ok=True)
        pep_pdb = wdir / "peptide.pdb"
        cx_pdb  = wdir / "complex.pdb"
        pdbqt_pose_to_pdb(pose_pdbqt, pep_pdb)
        build_complex_pdb(APO_RECEPTOR_PDB, pep_pdb, cx_pdb)
        r = haddock3_score(cx_pdb)
        if r["ok"] and r["haddock_score"] is not None:
            print(f"    HADDOCK-score = {r['haddock_score']:+.2f}   energies = {r['energies']}")
        else:
            print(f"    ✗ scoring failed; stdout tail: {r.get('stdout_tail','')[:160]}")
        rows.append({"peptide": label, "comment": comment,
                     "haddock_score": r.get("haddock_score"),
                     **{f"E_{k}": v for k, v in r.get("energies", {}).items()}})

    if rows:
        # Coerce field set
        all_keys = set()
        for r in rows: all_keys.update(r.keys())
        all_keys = sorted(all_keys)
        csv_path = OUT / "haddock3_scores.csv"
        with csv_path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=all_keys); w.writeheader(); w.writerows(rows)
        print(f"\n  → {csv_path}")

        # Δ vs scrambled
        scr = next((r for r in rows if r["peptide"] == "LR8_scrambled"), None)
        can = next((r for r in rows if r["peptide"] == "LR8_canonical"), None)
        if scr and can and scr["haddock_score"] and can["haddock_score"]:
            d = can["haddock_score"] - scr["haddock_score"]
            print(f"\n  Δ HADDOCK-score (canonical − scrambled) = {d:+.2f}")
            print(f"    canonical = {can['haddock_score']:+.2f}, scrambled = {scr['haddock_score']:+.2f}")
            print("    → " + ("canonical IS more favourable than scrambled (negative Δ)" if d < -1
                              else "within noise; null result confirmed under HADDOCK scoring"))

if __name__ == "__main__":
    main()
