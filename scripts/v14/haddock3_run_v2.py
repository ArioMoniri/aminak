#!/usr/bin/env python3
"""Phase 14g v2 — HADDOCK3-score with PROPERLY chain-segregated complex PDBs.

v1 issue: I concatenated receptor + peptide both as "chain Z" → haddock3-score
silently treated them as a single chain and returned receptor-only scores.

v2 fix: rebuild each complex PDB with
  - chain A   = receptor (all protein residues, renumbered if needed)
  - chain B   = peptide  (all Vina pose atoms; rewrite resname to standard
                          amino acids so HADDOCK recognises them as a peptide)
  - TER between chains
  - END at file end

Then re-run haddock3-score.
"""
from __future__ import annotations
import subprocess, csv, re, shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "14_inhibitor_design" / "07_advanced_methods" / "haddock3_run_v2"
OUT.mkdir(parents=True, exist_ok=True)
HADDOCK3_SCORE = Path.home() / "Library/Python/3.14/bin/haddock3-score"

S3_DIR = REPO / "14_inhibitor_design" / "03_dimer_interface" / "docked"
APO_RECEPTOR_PDB = REPO / "06f_receptor_fixed" / "dimer_noH.pdb"

PEPTIDES = [
    ("LR8_canonical",     "LSCQLYQR", S3_DIR / "LR8_LSCQLYQR_seed42.pdbqt", "Cardinale-2011 octapeptide"),
    ("LR8_scrambled",     "QLCRQSYL", S3_DIR / "LR8_scrambled_QLCRQSYL_seed42.pdbqt", "shuffled-seq control"),
    ("LR4_pos1",          "LSCQ",     S3_DIR / "LR_4mer_pos1_LSCQ_seed42.pdbqt", "4-mer LSCQ"),
    ("LR4_pos3",          "CQLY",     S3_DIR / "LR_4mer_pos3_CQLY_seed42.pdbqt", "4-mer CQLY"),
    ("LR4_pos5",          "LYQR",     S3_DIR / "LR_4mer_pos5_LYQR_seed42.pdbqt", "4-mer LYQR"),
]

AA_3 = {"L":"LEU","S":"SER","C":"CYS","Q":"GLN","Y":"TYR","R":"ARG",
        "G":"GLY","A":"ALA","V":"VAL","I":"ILE","M":"MET","F":"PHE",
        "W":"TRP","P":"PRO","H":"HIS","K":"LYS","D":"ASP","E":"GLU",
        "N":"ASN","T":"THR"}


def receptor_to_chain_a(pdb_path: Path) -> str:
    """Rewrite receptor PDB lines so col 22 = 'A' (single chain)."""
    out = []
    n_atoms = 0
    for ln in pdb_path.read_text().splitlines():
        if ln.startswith("ATOM"):
            # force chain A
            ln = ln[:21] + "A" + ln[22:]
            out.append(ln)
            n_atoms += 1
        elif ln.startswith("TER"):
            out.append("TER")
    # Determine last receptor resid for peptide-resid offset
    last_resid = 0
    for ln in out:
        if ln.startswith("ATOM"):
            try: last_resid = max(last_resid, int(ln[22:26]))
            except ValueError: pass
    return "\n".join(out), last_resid, n_atoms


def peptide_pose_to_chain_b(pose_pdbqt: Path, sequence: str, start_resid: int) -> str:
    """Rewrite a Vina pose PDBQT so each peptide residue maps to chain B with
    consecutive resids and correct 3-letter resnames.

    This is approximate — Vina pose PDBQTs from RDKit `MolFromSequence` builds
    don't carry residue boundaries reliably. We split the heavy-atom list into
    len(sequence) roughly-equal groups and assign residue ids serially.
    """
    # Extract MODEL 1 ATOM/HETATM lines
    atoms = []
    in_m1 = False
    for ln in pose_pdbqt.read_text().splitlines():
        if ln.startswith("MODEL 1"): in_m1 = True; continue
        if ln.startswith("ENDMDL") and in_m1: break
        if in_m1 and (ln.startswith("ATOM") or ln.startswith("HETATM")):
            # strip AD-type past col 78
            atoms.append(ln[:66])
    if not atoms: return ""
    n_atoms = len(atoms); n_res = len(sequence)
    atoms_per_res = max(1, n_atoms // n_res)
    out = []
    for i, ln in enumerate(atoms):
        ri = min(i // atoms_per_res, n_res - 1)
        aa1 = sequence[ri]
        aa3 = AA_3.get(aa1, "GLY")
        resid = start_resid + 100 + ri  # offset so they don't clash with receptor resids
        # rewrite: ATOM/HETATM, name (col 13-16), resname (17-20), chain (22) = 'B',
        #         resid (23-26)
        head = "ATOM  "
        atom_no = i + 1
        # keep the atom-name as-is (the meeko/RDKit output names are fine for haddock3-score)
        # only force chain + resname + resid
        name_field = ln[12:16]
        new = f"{head}{atom_no:5d} {name_field}{aa3:>3s} B{resid:4d}{ln[26:]}"
        out.append(new)
    out.append("TER")
    return "\n".join(out)


def build_complex(label: str, sequence: str, pose_pdbqt: Path) -> Path:
    wdir = OUT / label; wdir.mkdir(parents=True, exist_ok=True)
    rec_str, last_resid, n_rec = receptor_to_chain_a(APO_RECEPTOR_PDB)
    pep_str = peptide_pose_to_chain_b(pose_pdbqt, sequence, last_resid)
    cx_pdb = wdir / "complex_chain_segregated.pdb"
    cx_pdb.write_text(rec_str + "\nTER\n" + pep_str + "\nEND\n")
    print(f"  complex assembled: {n_rec} chain-A atoms + {len(pep_str.splitlines())-1} chain-B atoms")
    return cx_pdb


def haddock3_score(complex_pdb: Path) -> dict:
    cmd = [str(HADDOCK3_SCORE), str(complex_pdb), "--run_dir",
           str(complex_pdb.parent / "_run"), "--full"]
    try:
        proc = subprocess.run(cmd, check=False, capture_output=True, timeout=600,
                              cwd=str(complex_pdb.parent))
    except subprocess.TimeoutExpired:
        return {"ok": False, "err": "timeout"}
    out = proc.stdout.decode() + proc.stderr.decode()
    score_m = re.search(r"HADDOCK\s*score\s*[=:]\s*([-\d.]+)", out, re.I)
    score = float(score_m.group(1)) if score_m else None
    energies = {}
    for term, regex in [("vdw",   r"v\s?d\s?w\s*[=:]\s*([-\d.]+)"),
                        ("elec",  r"e\s?l\s?e\s?c(?:trostatic)?\s*[=:]\s*([-\d.]+)"),
                        ("desolv",r"desolv\s*[=:]\s*([-\d.]+)"),
                        ("air",   r"air\s*[=:]\s*([-\d.]+)"),
                        ("bsa",   r"bsa\s*[=:]\s*([-\d.]+)")]:
        em = re.search(regex, out, re.I)
        if em:
            try: energies[term] = float(em.group(1))
            except ValueError: pass
    return {"ok": True, "haddock_score": score, "energies": energies,
            "stdout_tail": out[-600:] if out else ""}


def main():
    print("=== Phase 14g v2 — HADDOCK3-score with chain-segregated complex PDBs ===")
    if not HADDOCK3_SCORE.exists():
        print(f"  ! haddock3-score not found"); return
    rows = []
    for label, sequence, pose_pdbqt, comment in PEPTIDES:
        print(f"\n--- {label}  ({comment}) — sequence {sequence} ---")
        if not pose_pdbqt.exists():
            print(f"  ! missing: {pose_pdbqt}"); continue
        cx = build_complex(label, sequence, pose_pdbqt)
        r = haddock3_score(cx)
        if r["ok"] and r.get("haddock_score") is not None:
            e = r["energies"]
            print(f"  HADDOCK-score = {r['haddock_score']:+.2f}    "
                  f"vdw={e.get('vdw','?')}  elec={e.get('elec','?')}  desolv={e.get('desolv','?')}  bsa={e.get('bsa','?')}")
        else:
            print(f"  ✗ score failed; stdout tail:")
            print(r.get("stdout_tail", "")[:500])
        rows.append({"peptide": label, "sequence": sequence, "comment": comment,
                     "haddock_score": r.get("haddock_score"),
                     **{f"E_{k}": v for k, v in r.get("energies", {}).items()}})

    if rows:
        fields = sorted(set(k for r in rows for k in r))
        csv_path = OUT / "haddock3_scores_v2.csv"
        with csv_path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
        print(f"\n  → {csv_path}")

        scr = next((r for r in rows if r["peptide"] == "LR8_scrambled"), None)
        can = next((r for r in rows if r["peptide"] == "LR8_canonical"), None)
        if scr and can and scr["haddock_score"] is not None and can["haddock_score"] is not None:
            d = can["haddock_score"] - scr["haddock_score"]
            print(f"\n  Δ HADDOCK-score (canonical − scrambled) = {d:+.2f}")
            print(f"    canonical = {can['haddock_score']:+.2f}, scrambled = {scr['haddock_score']:+.2f}")
            if abs(d) > 1.0:
                print(f"    ★ above noise — engines DIFFER on canonical vs scrambled now")
            else:
                print(f"    within noise — null result confirmed even with proper chain segregation")


if __name__ == "__main__":
    main()
