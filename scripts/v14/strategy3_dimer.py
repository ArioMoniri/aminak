#!/usr/bin/env python3
"""Phase 14 — Strategy 3 (dimer interface) — SCOPED execution.

A3  compute chain-A↔B 4 Å contact map → interface residues + box centre
B   peptide source: LR octapeptide from Cardinale 2011 (LSCQLYQR by literature).
    Short mimetics: pentamers (5 residues) from the LR sequence.
D   primary: HPEPDOCK (web, unreachable per pre-flight) → CABS-dock (web) → fragment-decompose Vina.
    For this scoped execution we use Vina-only fragment decomposition (most-portable path).
F   PPI metrics: BSA (dimer ± ligand) via freesasa, interface contact count via MDAnalysis.
"""
from __future__ import annotations
import json, sys, csv, time, traceback, subprocess
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (REPO, PHASE, RECEPTOR_APO, COFACTOR_A, COFACTOR_B,
                    prep_ligand_from_sdf, vina_dock, compound_descriptors)

STRAT = PHASE / "03_dimer_interface"
LOG = STRAT / "strategy3_runlog.txt"

# Literature LR peptide (Cardinale 2011) — building blocks for the disruptor
LR_SEQUENCE = "LSCQLYQR"  # 8-mer (will need HPEPDOCK normally)


def log(msg: str):
    print(msg, flush=True)
    STRAT.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as f: f.write(time.strftime("%H:%M:%S ") + msg + "\n")


# ---------- A3: 4 Å contact map ----------
def step_A3() -> dict:
    log("\n=== A3 — chain A↔B 4 Å contact map ===")
    # Convert apo PDBQT to a minimal PDB for MDAnalysis
    apo_pdb = STRAT / "apo_for_mda.pdb"
    if not apo_pdb.exists():
        lines = []
        cur_chain = "A"
        for ln in RECEPTOR_APO.read_text().splitlines():
            if ln.startswith(("ATOM","HETATM")):
                # patch chain id into col 22
                lines.append(ln[:21] + cur_chain + ln[22:66])
            elif ln.startswith("TER"):
                lines.append("TER\n")
                cur_chain = "B"
        apo_pdb.write_text("\n".join(lines) + "\nEND\n")
    import MDAnalysis as mda
    u = mda.Universe(str(apo_pdb))
    chA = u.select_atoms("chainID A and not name H*")
    chB = u.select_atoms("chainID B and not name H*")
    log(f"  chain A heavy atoms: {len(chA)}, chain B heavy atoms: {len(chB)}")
    if len(chA) == 0 or len(chB) == 0:
        # patching chains failed — assume all atoms in single segment, split by residue range
        all_atoms = u.select_atoms("protein and not name H*")
        # heuristic: first half = A, second half = B
        n = len(all_atoms) // 2
        chA = all_atoms[:n]; chB = all_atoms[n:]
        log(f"  (chain split heuristic — A: {len(chA)} atoms, B: {len(chB)} atoms)")
    # Build position arrays
    pA = chA.positions; pB = chB.positions
    # Pairwise min distance per residue (slow but ok for <300 res)
    from scipy.spatial.distance import cdist
    d = cdist(pA, pB)  # N_A × N_B
    near = d <= 4.0
    interface_A_resids = sorted(set(int(chA[i].resid) for i in range(len(chA)) if near[i].any()))
    interface_B_resids = sorted(set(int(chB[j].resid) for j in range(len(chB)) if near[:, j].any()))
    # Box centre = midpoint of chain-A interface Cα centroid and chain-B interface Cα centroid
    chA_ca = u.select_atoms("chainID A and name CA")
    chB_ca = u.select_atoms("chainID B and name CA")
    cenA = np.mean([a.position for a in chA_ca if int(a.resid) in interface_A_resids], axis=0) if chA_ca else None
    cenB = np.mean([a.position for a in chB_ca if int(a.resid) in interface_B_resids], axis=0) if chB_ca else None
    if cenA is None or cenB is None:
        # fallback to all-atom centroid
        cenA = np.mean([a.position for a in chA if int(a.resid) in interface_A_resids], axis=0)
        cenB = np.mean([a.position for a in chB if int(a.resid) in interface_B_resids], axis=0)
    centre = (cenA + cenB) / 2
    box = {"cx": float(centre[0]), "cy": float(centre[1]), "cz": float(centre[2]),
           "sx": 26.0, "sy": 22.0, "sz": 22.0}
    result = {"interface_A": interface_A_resids, "interface_B": interface_B_resids,
              "n_interface_A": len(interface_A_resids), "n_interface_B": len(interface_B_resids),
              "centre": centre.tolist(), "box": box}
    log(f"  interface A residues ({len(interface_A_resids)}): {interface_A_resids[:30]}{'…' if len(interface_A_resids)>30 else ''}")
    log(f"  interface B residues ({len(interface_B_resids)}): {interface_B_resids[:30]}{'…' if len(interface_B_resids)>30 else ''}")
    log(f"  box centre = {centre}  size = (26, 22, 22) Å")
    (STRAT / "A3_contact_map").mkdir(parents=True, exist_ok=True)
    (STRAT / "A3_contact_map" / "result.json").write_text(json.dumps(result, indent=2))
    return result


# ---------- B: build peptide ligands ----------
def step_B(seq: str = LR_SEQUENCE) -> list[dict]:
    log(f"\n=== B — build peptide library from LR sequence {seq!r} ===")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    compounds = []
    # canonical 8-mer
    canonical = {"name": f"LR8_{seq}", "sequence": seq, "length": len(seq), "kind": "canonical"}
    compounds.append(canonical)
    # scrambled control (numpy seed=42 permutation)
    rng = np.random.default_rng(42)
    scrambled = "".join(rng.permutation(list(seq)).tolist())
    compounds.append({"name": f"LR8_scrambled_{scrambled}", "sequence": scrambled, "length": 8, "kind": "scrambled_control"})
    # fragment decomposition: overlapping 4-mers
    for i in range(len(seq) - 3):
        sub = seq[i:i+4]
        compounds.append({"name": f"LR_4mer_pos{i+1}_{sub}", "sequence": sub, "length": 4, "kind": "fragment_4mer"})
    # build 3D structures via RDKit Chem.MolFromSequence
    for c in compounds:
        sdf = STRAT / "ligands" / f"{c['name']}.sdf"
        sdf.parent.mkdir(parents=True, exist_ok=True)
        if sdf.exists() and sdf.stat().st_size > 0:
            c["sdf_ok"] = True; continue
        try:
            mol = Chem.MolFromSequence(c["sequence"])
            if mol is None:
                c["sdf_ok"] = False; c["err"] = "MolFromSequence None"; continue
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            try: AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
            except Exception: pass
            Chem.MolToMolFile(mol, str(sdf))
            c["sdf_ok"] = True
            # descriptors
            mol_noh = Chem.RemoveHs(mol)
            c["mw"] = sum(a.GetMass() for a in mol_noh.GetAtoms())
            c["n_heavy"] = mol_noh.GetNumHeavyAtoms()
            log(f"  ✓ {c['name']:35s} (len={c['length']}, MW={c['mw']:.1f}, heavy={c['n_heavy']})")
        except Exception as e:
            c["sdf_ok"] = False; c["err"] = str(e)[:200]
            log(f"  ✗ {c['name']:35s} build failed: {e}")
    (STRAT / "compounds.json").write_text(json.dumps(compounds, indent=2))
    return compounds


# ---------- C: prep ----------
def step_C(compounds):
    log("\n=== C — peptide ligand prep ===")
    for c in compounds:
        if not c.get("sdf_ok"): continue
        sdf = STRAT / "ligands" / f"{c['name']}.sdf"
        pdbqt = STRAT / "ligands" / f"{c['name']}.pdbqt"
        if pdbqt.exists() and pdbqt.stat().st_size > 0:
            c["prep_ok"] = True; continue
        r = prep_ligand_from_sdf(sdf, pdbqt)
        c["prep_ok"] = r["ok"]; c["prep_stage"] = r.get("stage")
        log(f"  {'✓' if r['ok'] else '✗'} {c['name']:35s} ({r.get('stage')})")
    return compounds


# ---------- D: Vina dock (fragment-decomp because HPEPDOCK/CABS-dock unreachable) ----------
def step_D(compounds, box, seeds=(42, 7)):
    log(f"\n=== D — Vina dock (HPEPDOCK unreachable → fragment-decomp + canonical 8-mer with Vina; quality caveat) ===")
    rows = []
    for c in compounds:
        if not c.get("prep_ok"): continue
        ligand = STRAT / "ligands" / f"{c['name']}.pdbqt"
        # peptides ≥6 residues with Vina = unreliable per Hassan 2017; we run anyway but flag
        unreliable = c["length"] >= 6
        # For 8-mers in Vina, use exh=4 (much faster, ranking-only signal); fragments stay at exh=16
        exh = 4 if c["length"] >= 6 else 16
        nmod = 5 if c["length"] >= 6 else 9
        for seed in seeds:
            out = STRAT / "docked" / f"{c['name']}_seed{seed}.pdbqt"
            r = vina_dock(RECEPTOR_APO, ligand, out, box, seed=seed, exhaustiveness=exh, num_modes=nmod, cpu=4, timeout_s=1200)
            row = {"peptide": c["name"], "sequence": c["sequence"], "length": c["length"], "kind": c["kind"],
                   "seed": seed, "engine": "Vina_fragment_decomp",
                   "top1": r.get("top1"), "top3_mean": r.get("top3_mean"),
                   "ok": r.get("ok"), "wall_s": r.get("wall_s"),
                   "vina_peptide_unreliable_flag": unreliable,
                   "mw": c.get("mw"), "n_heavy": c.get("n_heavy")}
            rows.append(row)
            log(f"  {c['name']:35s} s={seed} top1={row.get('top1')} unreliable={unreliable} {row.get('wall_s')}s")
    csv_path = STRAT / "results_raw.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    return rows


# ---------- G + F: aggregate + canonical-vs-scrambled control ----------
def step_G(rows):
    log("\n=== G — canonical-vs-scrambled control + per-kind ranking ===")
    import collections
    by = collections.defaultdict(list)
    for r in rows:
        if r.get("top1") is None: continue
        by[r["peptide"]].append(r["top1"])
    summary = []
    for pep, ts in by.items():
        meta = next(r for r in rows if r["peptide"] == pep)
        mean = sum(ts)/len(ts)
        s = {"peptide": pep, "sequence": meta["sequence"], "length": meta["length"], "kind": meta["kind"],
             "n_seeds": len(ts), "top1_mean": float(mean), "top1_min": float(min(ts)),
             "top1_sd": float((sum((t-mean)**2 for t in ts)/len(ts))**0.5),
             "engine": meta["engine"], "unreliable_flag": meta["vina_peptide_unreliable_flag"],
             "mw": meta.get("mw"), "n_heavy": meta.get("n_heavy")}
        summary.append(s)
    # scrambled-control comparison
    canonical = next((s for s in summary if s["kind"] == "canonical"), None)
    scrambled = next((s for s in summary if s["kind"] == "scrambled_control"), None)
    log("\n  Canonical-vs-scrambled control (R2 sign-off #3):")
    if canonical and scrambled:
        delta = canonical["top1_mean"] - scrambled["top1_mean"]
        log(f"    canonical {canonical['peptide']}  top1={canonical['top1_mean']:+.3f}")
        log(f"    scrambled {scrambled['peptide']}  top1={scrambled['top1_mean']:+.3f}")
        log(f"    Δ = {delta:+.3f}  (canonical should be ≤ scrambled − 2 kcal/mol for specificity)")
        specificity_ok = delta <= -2.0
        log(f"    specificity_ok = {specificity_ok}")
        for s in summary: s["specificity_vs_scrambled"] = delta if s["kind"]=="canonical" else None
    csv_path = STRAT / "results_summary.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
    log(f"  → {csv_path}")
    summary.sort(key=lambda x: (x.get("top1_mean") or 0))
    log("\n  Ranking (most negative = best):")
    for s in summary:
        log(f"    {s['peptide']:35s}  top1={s['top1_mean']:+.3f}  kind={s['kind']:18s}")


def main():
    STRAT.mkdir(parents=True, exist_ok=True)
    LOG.write_text(f"Phase 14 Strategy 3 — started {time.ctime()}\n")
    t0 = time.time()
    try:
        a3 = step_A3()
        compounds = step_B(LR_SEQUENCE)
        compounds = step_C(compounds)
        rows = step_D(compounds, a3["box"])
        step_G(rows)
        log(f"=== DONE in {(time.time()-t0)/60:.1f} min ===")
    except Exception as e:
        log(f"!!! ABORT: {e}"); log(traceback.format_exc()); raise

if __name__ == "__main__":
    main()
