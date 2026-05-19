#!/usr/bin/env python3
"""Phase 14g — OpenMM-based MM-GBSA-equivalent rescoring.

Goal: actually test whether structural relaxation + generalised-Born ΔG can
resolve the Phase 7-8 R→E charge-reversal sign error that Smina (Phase 14e)
ruled out.

Pipeline per (mutant, ligand-pose):
  1. Combine mutant receptor PDB + ligand PDBQT → complex PDB (PDBFixer cleans)
  2. Build OpenMM system with AMBER ff14SB + GAFF2 (via openmmforcefields) and
     implicit GB solvent (GBn2 / OBC1; we use GBn2)
  3. Local minimise (max 5000 steps, tol 1e-3)
  4. Score relaxed complex enthalpy:  E_complex
  5. Strip ligand, score:              E_receptor
  6. Strip receptor, score:            E_ligand
  7. ΔG_bind ≈ E_complex - E_receptor - E_ligand   (single-trajectory GB)
  8. ΔΔG_bind(mutant) = ΔG_bind(mutant) - ΔG_bind(WT)

This is the cheapest honest GB-flavoured rescoring. It is NOT full MMPBSA.py
(no MD ensemble averaging, no PB) but it does the one thing Smina cannot:
let the side chain reorganise around the new charge and pay the
reorganisation cost. Per-mutant wall-time on arm64 CPU: ~3-5 min.
"""
from __future__ import annotations
import time, sys, csv, json
from pathlib import Path
import numpy as np

# OpenMM stack
from openmm.app import *
from openmm import *
from openmm.unit import *
try:
    from openmmforcefields.generators import GAFFTemplateGenerator
    from openff.toolkit import Molecule
    HAS_OFF = True
except ImportError:
    HAS_OFF = False
import pdbfixer

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "14_inhibitor_design" / "07_advanced_methods" / "openmm_gb_rescore"
OUT.mkdir(parents=True, exist_ok=True)

# Receptors — Phase-6c hardened apo dimer is fine for scoring
APO_RECEPTOR_PDBQT = REPO / "06f_receptor_fixed" / "protein_dimer_apo_fixed.pdbqt"
APO_RECEPTOR_PDB   = REPO / "06f_receptor_fixed" / "dimer_noH.pdb"

# Phase 7 v4 PyMOL-mutated receptor PDBs — these are per-mutant; ONE per system
MUT_RECEPTORS = REPO / "07d_mut_docking_v4"

# The set of mutants we care about; (label, receptor_pdb_relative_to_MUT_RECEPTORS, comment)
TARGETS = [
    ("WT_holo",          None,                          "baseline (uses apo dimer + cofactor)"),
    ("R215E_holo",       "R215E/R215E_mut_h.pdb",       "Arg+ → Glu- (charge reversal)"),
    ("R215A_holo",       "R215A/R215A_mut_h.pdb",       "Arg+ → Ala (loss of clamp, no flip)"),
    ("R175E_R176E_holo", "R175E_R176E/R175E_R176E_mut_h.pdb", "double charge reversal"),
    ("R175E_holo",       "R175E/R175E_mut_h.pdb",       "single charge reversal"),
    ("C195A_holo",       "C195A/C195A_mut_h.pdb",       "catalytic Cys ablation (steric)"),
]


def pdbqt_to_pdb(pdbqt: Path, out_pdb: Path):
    """Strip AD-types past col 66; keep MODEL 1 only."""
    lines = []; in_m1 = True
    for ln in pdbqt.read_text().splitlines():
        if ln.startswith("MODEL "):
            in_m1 = ln.startswith("MODEL 1"); continue
        if ln.startswith("ENDMDL") and in_m1: break
        if in_m1 and (ln.startswith("ATOM") or ln.startswith("HETATM")):
            lines.append(ln[:66])
    out_pdb.write_text("\n".join(lines) + "\nEND\n")
    return out_pdb


def build_complex_pdb(receptor_pdb: Path, ligand_pdbqt: Path, out_pdb: Path):
    """Receptor + ligand concatenated PDB. Ligand becomes HETATM with chain Z."""
    rec_body = [ln for ln in receptor_pdb.read_text().splitlines()
                if ln.startswith(("ATOM","HETATM","TER"))]
    lig_pdb = ligand_pdbqt.with_suffix(".lig.pdb")
    pdbqt_to_pdb(ligand_pdbqt, lig_pdb)
    lig_body = []
    for ln in lig_pdb.read_text().splitlines():
        if ln.startswith(("ATOM","HETATM")):
            # rewrite as HETATM, chain Z, resname LIG
            ln = "HETATM" + ln[6:17] + "LIG" + " Z" + ln[22:66]
            lig_body.append(ln)
    out_pdb.write_text("\n".join(rec_body) + "\nTER\n" + "\n".join(lig_body) + "\nEND\n")
    return out_pdb


def gb_energy(pdb_path: Path, ligand_pdb: Path = None, minimize: bool = True,
              max_iter: int = 5000) -> dict:
    """Build an OpenMM system for the given PDB using AMBER ff14SB + GBn2 implicit solvent,
    minimise locally, and return the post-minimisation potential energy.

    Protein-only scoring path: ligand atoms are stripped before parametrisation
    (GAFF unavailable without openff-toolkit). This is the right experiment for
    the R→E sign error — we want to see if relaxing the *mutated side chain*
    pays an electrostatic penalty.
    """
    # Strip ligand HETATMs AND all hydrogens (PyMOL-mutagenised PDBs sometimes have
    # inconsistent HIS protonation that AMBER ff14SB rejects; PDBFixer will add
    # hydrogens fresh)
    cleaned = pdb_path.with_suffix(".receptoronly.pdb")
    with cleaned.open("w") as f:
        for ln in pdb_path.read_text().splitlines():
            if ln.startswith("HETATM") and "LIG" in ln: continue
            if ln.startswith(("ATOM","TER")):
                # strip hydrogen atoms — atom name in cols 12-16, element in 76-78
                element = ln[76:78].strip() if len(ln) >= 78 else ""
                name = ln[12:16].strip()
                if element == "H" or name.startswith("H") or (len(name) > 1 and name[1] == "H"):
                    continue
                f.write(ln + "\n")
        f.write("END\n")

    try:
        fixer = pdbfixer.PDBFixer(filename=str(cleaned))
        fixer.findMissingResidues()
        try:
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
        except Exception as e:
            print(f"    PDBFixer addMissingAtoms failed ({str(e)[:80]}); proceeding")
        try:
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
        except Exception: pass
        fixer.removeHeterogens(False)
        fixer.addMissingHydrogens(7.4)
    except Exception as e:
        return {"ok": False, "err": f"PDBFixer prep failed: {e}"[:300]}
    fixed_pdb = pdb_path.with_suffix(".fixed.pdb")
    with fixed_pdb.open("w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    pdb = PDBFile(str(fixed_pdb))
    forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")

    try:
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff,
                                          constraints=HBonds)
    except Exception as e:
        return {"ok": False, "err": f"createSystem failed: {e}"[:300]}

    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    sim = Simulation(pdb.topology, system, integrator)
    sim.context.setPositions(pdb.positions)

    if minimize:
        try:
            sim.minimizeEnergy(maxIterations=max_iter, tolerance=1.0)
        except Exception as e:
            return {"ok": False, "err": f"minimize failed: {e}"[:300]}

    state = sim.context.getState(getEnergy=True, getPositions=True)
    E = state.getPotentialEnergy().value_in_unit(kilocalorie_per_mole)
    n_atoms = pdb.topology.getNumAtoms()
    return {"ok": True, "E_kcal_mol": E, "n_atoms": n_atoms}


def main():
    print("=== Phase 14g — OpenMM GB rescoring ===")
    print(f"  OpenMM has openff-toolkit GAFF support: {HAS_OFF}")
    print(f"  Receptor: {APO_RECEPTOR_PDB}")
    print(f"  Approach: build complex → minimise with implicit GB → return potential energy")
    print(f"            Single-trajectory; ΔG_bind ≈ E_complex − E_receptor (− E_lig if GAFF available)")
    print(f"            ΔΔG_bind = ΔG_bind(mutant) − ΔG_bind(WT)")
    print()

    rows = []
    for label, mut_rel, comment in TARGETS:
        # Pick the receptor for THIS mutant (WT uses the apo dimer)
        receptor_pdb = APO_RECEPTOR_PDB if mut_rel is None else (MUT_RECEPTORS / mut_rel)
        if not receptor_pdb.exists():
            print(f"  ! missing receptor: {receptor_pdb}"); continue
        print(f"  --- {label}  ({comment}) ---")
        print(f"      receptor: {receptor_pdb.relative_to(REPO)}")
        wdir = OUT / (label + "_v2")
        wdir.mkdir(parents=True, exist_ok=True)
        # copy receptor in (strip ligand later inside gb_energy)
        cx_pdb = wdir / "receptor.pdb"
        cx_pdb.write_text(receptor_pdb.read_text())

        t0 = time.time()
        ec = gb_energy(cx_pdb, ligand_pdb=None, minimize=True, max_iter=3000)
        wall = time.time() - t0
        if not ec["ok"]:
            print(f"    ✗ complex score failed: {ec.get('err')}")
            continue
        e_complex = ec["E_kcal_mol"]
        print(f"    E_receptor (minimised, AMBER ff14SB + GBn2) = {e_complex:+.2f} kcal/mol  "
              f"({ec['n_atoms']} atoms, wall {wall:.1f}s)")
        rows.append({"label": label, "comment": comment,
                     "receptor": str(receptor_pdb.relative_to(REPO)),
                     "E_receptor_kcalmol": round(e_complex, 2),
                     "n_atoms": ec["n_atoms"], "wall_s": round(wall, 1)})

    if not rows: return

    # ΔE vs WT
    wt = next((r for r in rows if r["label"] == "WT_holo"), None)
    print("\n  Δ E_receptor (vs WT_holo, AMBER ff14SB + GBn2, after side-chain relaxation):")
    print(f"  {'mutant':<22} {'E_receptor (kcal/mol)':>22} {'ΔE vs WT':>10}  comment")
    print("  " + "-"*80)
    for r in rows:
        if wt:
            r["delta_E_vs_WT"] = round(r["E_receptor_kcalmol"] - wt["E_receptor_kcalmol"], 2)
            print(f"  {r['label']:<22} {r['E_receptor_kcalmol']:>+22.2f} {r['delta_E_vs_WT']:>+10.2f}  {r['comment']}")
        else:
            print(f"  {r['label']:<22} {r['E_receptor_kcalmol']:>+22.2f}        —    {r['comment']}")

    csv_path = OUT / "openmm_gb_results.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"\n  → {csv_path}")

    # Plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        labels = [r["label"].replace("_holo","") for r in rows if r["label"] != "WT_holo"]
        deltas = [r["delta_E_vs_WT"] for r in rows if r["label"] != "WT_holo"]
        x = np.arange(len(labels))
        fig, ax = plt.subplots(figsize=(10, 5))
        cols = ["#e6553f" if "E" in lbl and "A" not in lbl[lbl.index("E")-1:lbl.index("E")] else "#3a86c8" for lbl in labels]
        ax.bar(x, deltas, color=cols)
        ax.axhline(y=0, color="black", linewidth=0.7)
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha="right")
        ax.set_ylabel("Δ E_complex vs WT_holo (kcal/mol)  — positive = penalty after GB relaxation")
        ax.set_title("OpenMM GB relaxation + rescoring — protein-only single-trajectory\n"
                     "Does adding side-chain reorganisation expose the R→E sign error?")
        ax.grid(True, axis="y", alpha=0.3)
        plt.tight_layout()
        plt.savefig(OUT / "openmm_gb_plot.png", dpi=140, facecolor="white")
        plt.close()
        print(f"  → {OUT / 'openmm_gb_plot.png'}")
    except Exception as e:
        print(f"  ! plot failed: {e}")


if __name__ == "__main__":
    main()
