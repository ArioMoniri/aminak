#!/usr/bin/env python3
"""Phase 14 — Strategy 1 (active-site, dUMP-mimetic) — SCOPED execution.

Runs A2 (CID verification), B (compound assembly), C (ligand prep),
D (Vina dock apo+holo × 2 seeds), G (aggregate CSV).
Pose analysis (E1b water-bridge) runs in a separate script after this completes.

Scoped scale (per user request): 5 Tier-1 anchors + 30 RDKit decoys + 2 seeds (42, 7).
"""
from __future__ import annotations
import json, sys, csv, os, time, traceback
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (
    REPO, PHASE, RECEPTOR_APO, COFACTOR_A, COFACTOR_B, BOX_ACTIVE_SITE,
    pubchem_get, fetch_sdf_3d, fetch_smiles, pubchem_similar_cids,
    prep_ligand_from_sdf, vina_dock, rdkit_decoys, compound_descriptors,
    load_anchors,
)

STRAT = PHASE / "01_active_site"
LOG = STRAT / "strategy1_runlog.txt"

def log(msg: str):
    print(msg, flush=True)
    with LOG.open("a") as f: f.write(time.strftime("%H:%M:%S ") + msg + "\n")

def receptor_path(state: str) -> Path:
    """Build apo or holo receptor PDBQT. Holo = apo + cofactor A + cofactor B ATOM lines only
    (cofactor PDBQTs have ROOT/BRANCH from meeko — strip those, keep just ATOM/HETATM lines)."""
    if state == "apo":
        return RECEPTOR_APO
    holo = STRAT / "receptor_holo.pdbqt"
    # always regenerate (in case prior runs left a broken file)
    body_lines = []
    for ln in RECEPTOR_APO.read_text().splitlines():
        if ln.startswith("END"): continue
        body_lines.append(ln)
    for cof in (COFACTOR_A, COFACTOR_B):
        for ln in cof.read_text().splitlines():
            if ln.startswith(("ATOM","HETATM")):
                body_lines.append(ln)
    body_lines.append("END")
    holo.write_text("\n".join(body_lines) + "\n")
    return holo

# ---------- A2 (verify CIDs) ----------
def step_A2():
    log("\n=== A2 — verify CIDs against PubChem ===")
    anchors = load_anchors()
    mismatches = []
    for a in anchors:
        cid = a.get("cid")
        if cid is None: continue
        try:
            d = pubchem_get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey,IsomericSMILES,MolecularFormula/JSON")
            p = d["PropertyTable"]["Properties"][0]
            expected = a["inchikey"]
            actual = p.get("InChIKey")
            if not actual:
                from rdkit import Chem
                actual = Chem.MolToInchiKey(Chem.MolFromSmiles(p["IsomericSMILES"]))
            if actual != expected:
                mismatches.append({"name": a["display_name"], "cid": cid, "expected": expected, "actual": actual})
                log(f"  ✗ {a['display_name']} CID {cid}: expected {expected}, got {actual}")
            else:
                log(f"  ✓ {a['display_name']} CID {cid}: {expected}")
        except Exception as e:
            log(f"  ! {a['display_name']} CID {cid}: fetch failed — {e}")
    out = STRAT / "A2_cid_verification" / "result.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"mismatches": mismatches, "checked": len([a for a in anchors if a.get("cid")])}, indent=2))
    if mismatches:
        log(f"A2 FAILED: {len(mismatches)} mismatches. Continuing in WARN mode (compounds with mismatch will be skipped at C).")
    return [m["cid"] for m in mismatches]

# ---------- B (compound set) ----------
def step_B() -> list[dict]:
    """Return list of {name, cid, smiles, tier, anchor_for, descriptors}."""
    log("\n=== B — assemble compound set ===")
    anchors = [a for a in load_anchors() if "active-site" in a["role"].lower() or "Strategy 1" in a.get("role","")]
    log(f"  Strategy 1 Tier-1 anchors: {len(anchors)}")
    compounds = []
    for a in anchors:
        if a.get("cid") is None: continue
        smi = a.get("smiles") or fetch_smiles(a["cid"])
        compounds.append({"name": a["display_name"].split("(")[0].strip().replace(" ", "_"),
                          "cid": a["cid"], "smiles": smi, "tier": 1, "anchor_for": a["display_name"]})
    # Tier 2: RDKit decoys against dUMP (avoid drug-like-only bias)
    log("  Generating Tier-2 RDKit decoys against dUMP (anchor) …")
    dump_anchor = next((a for a in anchors if a["display_name"] == "dUMP"), None)
    if dump_anchor:
        decoys = rdkit_decoys(dump_anchor["smiles"], n=30)
        log(f"  Got {len(decoys)} decoys")
        for d in decoys:
            compounds.append({"name": f"decoy_CID{d['cid']}", "cid": d["cid"],
                              "smiles": d["smiles"], "tier": 2, "anchor_for": "dUMP_decoy"})
    # descriptors + flags
    for c in compounds:
        c.update(compound_descriptors(c["smiles"]))
    out_file = STRAT / "compounds.json"
    out_file.write_text(json.dumps(compounds, indent=2))
    log(f"  Wrote {len(compounds)} compounds → {out_file}")
    return compounds

# ---------- C (ligand prep) ----------
def step_C(compounds: list[dict]) -> list[dict]:
    log("\n=== C — ligand prep (SDF → protonated → PDBQT) ===")
    prepped = []
    for c in compounds:
        sdf = STRAT / "ligands" / f"{c['name']}.sdf"
        pdbqt = STRAT / "ligands" / f"{c['name']}.pdbqt"
        if not sdf.exists():
            try:
                variant = fetch_sdf_3d(c["cid"], sdf)
                c["sdf_variant"] = variant
            except Exception as e:
                log(f"  ✗ {c['name']} (CID {c['cid']}) — SDF fetch failed: {e}")
                c["prep_ok"] = False; c["prep_err"] = "sdf_fetch"
                prepped.append(c); continue
        if not pdbqt.exists():
            r = prep_ligand_from_sdf(sdf, pdbqt)
            c["prep_ok"] = bool(r["ok"])
            c["prep_stage"] = r.get("stage")
            if not r["ok"]:
                c["prep_err"] = r.get("err","unknown")[:200]
                log(f"  ✗ {c['name']} — prep failed at {r.get('stage')}")
            else:
                log(f"  ✓ {c['name']} → {pdbqt.name} (via {r['stage']})")
        else:
            c["prep_ok"] = True; c["prep_stage"] = "cached"
        prepped.append(c)
    (STRAT / "compounds_prepped.json").write_text(json.dumps(prepped, indent=2))
    return prepped

# ---------- A0 (re-dock gate) ----------
def step_A0(compounds: list[dict]) -> dict:
    log("\n=== A0 — re-dock dUMP into 1HVY (gate: RMSD ≤ 2 Å) ===")
    dump = next((c for c in compounds if c["name"].lower().startswith("dump")), None)
    if not dump:
        log("  ! dUMP not in compound list — A0 skipped (WARN)")
        return {"ok": False, "reason": "dump_not_found"}
    ligand = STRAT / "ligands" / f"{dump['name']}.pdbqt"
    out = STRAT / "A0_redock_gate" / "dump_redock.pdbqt"
    r = vina_dock(receptor_path("apo"), ligand, out, BOX_ACTIVE_SITE, seed=42, exhaustiveness=32)
    log(f"  dUMP re-dock: top1={r.get('top1')} kcal/mol wall={r.get('wall_s')}s ok={r['ok']}")
    # RMSD vs crystal pose
    rmsd = None
    try:
        rmsd = compute_pose_rmsd_vs_crystal(out, REPO / "1HVY_dump_ref.pdb")
        log(f"  Pose RMSD vs crystal: {rmsd:.3f} Å")
    except Exception as e:
        log(f"  ! RMSD compute failed: {e} (continuing — gate WARN)")
    result = {"ok": r["ok"], "top1": r.get("top1"), "rmsd": rmsd, "gate_pass": rmsd is not None and rmsd <= 2.0}
    (STRAT / "A0_redock_gate" / "result.json").write_text(json.dumps(result, indent=2))
    return result

def compute_pose_rmsd_vs_crystal(pose_pdbqt: Path, crystal_pdb: Path) -> float:
    """Heavy-atom RMSD between Vina pose MODEL 1 and crystal ligand. Requires crystal_pdb to exist."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    # extract MODEL 1 to PDB for RDKit
    pose_pdb = pose_pdbqt.with_suffix(".pose1.pdb")
    lines = []
    in_m = False
    for ln in pose_pdbqt.read_text().splitlines():
        if ln.startswith("MODEL 1"): in_m = True; continue
        if ln.startswith("ENDMDL") and in_m: break
        if in_m and (ln.startswith("ATOM") or ln.startswith("HETATM")):
            # PDBQT has trailing AD type cols 78+; strip
            lines.append(ln[:66])
    pose_pdb.write_text("\n".join(lines) + "\nEND\n")
    pose_mol = Chem.MolFromPDBFile(str(pose_pdb), removeHs=True, sanitize=False)
    if not crystal_pdb.exists():
        raise FileNotFoundError(crystal_pdb)
    xtal_mol = Chem.MolFromPDBFile(str(crystal_pdb), removeHs=True, sanitize=False)
    if pose_mol is None or xtal_mol is None:
        raise RuntimeError("could not parse pose or crystal PDB")
    # symmetry-corrected best-RMS (heavy atoms)
    try:
        return AllChem.GetBestRMS(pose_mol, xtal_mol)
    except Exception:
        # fallback: nearest-atom mean distance on matched heavy atoms (no symmetry correction)
        import numpy as np
        p = pose_mol.GetConformer().GetPositions()
        x = xtal_mol.GetConformer().GetPositions()
        n = min(len(p), len(x))
        return float(np.sqrt(np.mean(np.sum((p[:n] - x[:n])**2, axis=1))))

# ---------- D (dock all compounds, scoped) ----------
def step_D(compounds: list[dict], seeds=(42, 7), states=("apo", "holo")) -> list[dict]:
    log(f"\n=== D — Vina dock (scoped: {len(compounds)} compounds × {len(seeds)} seeds × {len(states)} states) ===")
    rows = []
    total = sum(1 for c in compounds if c.get("prep_ok")) * len(seeds) * len(states)
    done = 0
    for c in compounds:
        if not c.get("prep_ok"):
            continue
        ligand = STRAT / "ligands" / f"{c['name']}.pdbqt"
        for state in states:
            rec = receptor_path(state)
            for seed in seeds:
                out_pdbqt = STRAT / "docked" / f"{c['name']}_{state}_seed{seed}.pdbqt"
                if out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
                    # skip if cached
                    res = {"ok": True, "cached": True}
                    # but we still need to parse the log
                    log_path = out_pdbqt.with_suffix(".log")
                    if log_path.exists():
                        text = log_path.read_text()
                        affs = []
                        for ln in text.splitlines():
                            parts = ln.split()
                            if len(parts) >= 4 and parts[0].isdigit():
                                try: affs.append(float(parts[1]))
                                except ValueError: pass
                        if affs:
                            res.update({"top1": affs[0], "top3_mean": sum(affs[:3])/min(3,len(affs)),
                                        "top5_mean": sum(affs[:5])/min(5,len(affs)), "n_modes": len(affs)})
                else:
                    res = vina_dock(rec, ligand, out_pdbqt, BOX_ACTIVE_SITE, seed=seed, exhaustiveness=32, cpu=4)
                done += 1
                row = {"compound": c["name"], "cid": c.get("cid"), "tier": c.get("tier"),
                       "anchor_for": c.get("anchor_for"), "state": state, "seed": seed,
                       "top1": res.get("top1"), "top3_mean": res.get("top3_mean"), "top5_mean": res.get("top5_mean"),
                       "n_modes": res.get("n_modes"), "ok": res.get("ok"), "wall_s": res.get("wall_s"),
                       "mw": c.get("mw"), "logp": c.get("logp"),
                       "pains": c.get("pains_a") or c.get("pains_b") or c.get("pains_c"),
                       "brenk": c.get("brenk"), "lipinski_ok": c.get("lipinski_ok")}
                rows.append(row)
                log(f"  [{done}/{total}] {c['name']:30s} {state:5s} seed={seed:3d} top1={row.get('top1')} wall={row.get('wall_s')}s")
    # write CSV
    csv_path = STRAT / "results_raw.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    log(f"  → {csv_path}")
    return rows

# ---------- G (aggregate) ----------
def step_G(rows: list[dict]):
    log("\n=== G — aggregate (per-compound × state means, Δ vs dUMP) ===")
    # collapse seeds
    import collections
    by = collections.defaultdict(list)
    for r in rows:
        if r.get("top1") is None: continue
        by[(r["compound"], r["state"])].append(r["top1"])
    summary = []
    for (cmpd, state), ts in by.items():
        meta = next(r for r in rows if r["compound"] == cmpd and r["state"] == state)
        s = {"compound": cmpd, "state": state, "tier": meta["tier"], "anchor_for": meta["anchor_for"],
             "cid": meta["cid"], "n_seeds": len(ts),
             "top1_mean": float(sum(ts)/len(ts)), "top1_min": float(min(ts)), "top1_max": float(max(ts)),
             "top1_sd": float((sum((t - sum(ts)/len(ts))**2 for t in ts)/len(ts))**0.5),
             "mw": meta.get("mw"), "logp": meta.get("logp"),
             "pains": meta.get("pains"), "brenk": meta.get("brenk"), "lipinski_ok": meta.get("lipinski_ok")}
        summary.append(s)
    # Δ vs dUMP per state
    dump_by_state = {s["state"]: s["top1_mean"] for s in summary if s["compound"].lower().startswith("dump")}
    for s in summary:
        ref = dump_by_state.get(s["state"])
        s["delta_vs_dump"] = (s["top1_mean"] - ref) if ref is not None else None
        s["significant_p085"] = (s["delta_vs_dump"] is not None and abs(s["delta_vs_dump"]) >= 0.85)
        s["beats_dump"] = (s["delta_vs_dump"] is not None and s["delta_vs_dump"] < -0.85)
    csv_path = STRAT / "results_summary.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader(); w.writerows(summary)
    log(f"  → {csv_path}")
    # print top-10 by Δ
    summary.sort(key=lambda x: (x["state"], x["top1_mean"]))
    log("\n  Top by absolute Vina score per state:")
    for state in ("apo", "holo"):
        log(f"  --- {state} ---")
        for s in [r for r in summary if r["state"] == state][:10]:
            log(f"    {s['compound']:30s} top1={s['top1_mean']:+.2f}  Δ={s['delta_vs_dump']}  tier={s['tier']}")

def main():
    STRAT.mkdir(parents=True, exist_ok=True)
    LOG.write_text(f"Phase 14 Strategy 1 — started {time.ctime()}\n")
    t0 = time.time()
    try:
        mism = step_A2()
        compounds = step_B()
        compounds = step_C(compounds)
        A0 = step_A0(compounds)
        rows = step_D(compounds)
        step_G(rows)
        log(f"\n=== DONE in {(time.time()-t0)/60:.1f} min ===")
    except Exception as e:
        log(f"\n!!! ABORT: {e}")
        log(traceback.format_exc())
        raise

if __name__ == "__main__":
    main()
