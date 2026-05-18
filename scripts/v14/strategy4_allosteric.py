#!/usr/bin/env python3
"""Phase 14 — Strategy 4 (allosteric / surface-hotspot) — SCOPED execution.

Pipeline:
  A.0  Run FPocket on chain A of the apo receptor → rank cavities by druggability
  A.1  Exclude cavities within 8 Å of the active-site centroid or cofactor centroid
  B    Pull a small ZINC15 fragment subset (or a built-in fallback list of fragment-like
       PubChem CIDs if ZINC15 is unreachable). MW ≤ 250, logP ≤ 3.5.
  C    Ligand prep
  D    Vina dock against each top allosteric cavity (1 seed = fragment screen scale)
  G    Aggregate; report which fragments land in which cavity
"""
from __future__ import annotations
import json, sys, csv, time, traceback, subprocess, shutil, urllib.request
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (REPO, PHASE, RECEPTOR_APO, COFACTOR_A, BOX_ACTIVE_SITE,
                    fetch_sdf_3d, fetch_smiles, prep_ligand_from_sdf, vina_dock,
                    compound_descriptors, pubchem_get)

STRAT = PHASE / "04_allosteric"
LOG = STRAT / "strategy4_runlog.txt"

FPOCKET = shutil.which("fpocket")
ACTIVE_SITE_CENTRE = np.array([BOX_ACTIVE_SITE["cx"], BOX_ACTIVE_SITE["cy"], BOX_ACTIVE_SITE["cz"]])

# Fallback fragment-like PubChem CIDs (MW < 250, drug-like fragments) if ZINC15 unreachable
def _manual_allosteric_candidates(apo_pdb: Path) -> list[dict]:
    """Fallback for fpocket failure (Strategy 4 S3): pick 3 surface residue centroids
    that are >= 15 A from the active-site centroid and >= 15 A from the cofactor centroid,
    drawn from chain A residues with high SASA (computed via freesasa).
    These are 'manual allosteric candidate boxes' — not druggability-ranked, just spatial.
    """
    import MDAnalysis as mda
    import freesasa
    u = mda.Universe(str(apo_pdb))
    # cofactor centre — read from cofactor A pdbqt
    cof_coords = []
    for ln in COFACTOR_A.read_text().splitlines():
        if ln.startswith(("ATOM","HETATM")):
            try: cof_coords.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
            except ValueError: continue
    cof_centre = np.array(cof_coords).mean(axis=0)
    log(f"  cofactor centre = {cof_centre}")
    # freesasa per-residue
    s = freesasa.calc(freesasa.Structure(str(apo_pdb)))
    res_areas = s.residueAreas()
    chains = list(res_areas.keys())
    chain_a_key = "A" if "A" in chains else (chains[0] if chains else None)
    if not chain_a_key:
        log("  ! freesasa produced no chains; aborting strategy 4")
        return []
    # candidate residues: highest SASA in chain A, distance gates passed
    chA_ca = u.select_atoms("chainID A and name CA")
    if len(chA_ca) == 0:
        chA_ca = u.select_atoms("name CA")[:len(u.select_atoms('name CA'))//2]
    candidates = []
    for resid_str, area_obj in res_areas[chain_a_key].items():
        resid = int(resid_str)
        ca = chA_ca.select_atoms(f"resid {resid}")
        if len(ca) == 0: continue
        pos = ca[0].position
        d_site = float(np.linalg.norm(pos - ACTIVE_SITE_CENTRE))
        d_cof  = float(np.linalg.norm(pos - cof_centre))
        if d_site < 15.0 or d_cof < 15.0: continue
        candidates.append({"resid": resid, "pos": pos, "sasa": float(area_obj.total),
                           "d_site": d_site, "d_cof": d_cof})
    # pick top 3 by SASA, but enforce mutual distance >= 10 A so they don't all cluster on one face
    candidates.sort(key=lambda x: -x["sasa"])
    picked = []
    for c in candidates:
        if all(np.linalg.norm(c["pos"] - p["pos"]) >= 10.0 for p in picked):
            picked.append(c)
        if len(picked) >= 3: break
    out = []
    for i, p in enumerate(picked):
        out.append({"cavity_id": f"manual_chainA_res{p['resid']}",
                    "centre": p["pos"].tolist(),
                    "Druggability Score": None,    # n/a for manual
                    "d_active_site": p["d_site"],
                    "d_cofactor": p["d_cof"],
                    "source_residue": p["resid"],
                    "residue_sasa": p["sasa"],
                    "fallback": True})
        log(f"    cavity {out[-1]['cavity_id']}: centre={p['pos'].tolist()}, d_site={p['d_site']:.1f}, d_cof={p['d_cof']:.1f}, sasa={p['sasa']:.1f}")
    return out


FALLBACK_FRAGMENTS = [
    1983, 2244, 2519, 2662, 3672, 4017, 4173, 5564, 6202, 6253,
    7032, 8094, 8654, 9087, 10257, 13558, 14888, 22311, 25245, 31703,
    35814, 39042, 49846, 60656, 65329, 67272, 67451, 69658, 92395, 9433,
]


def log(msg: str):
    print(msg, flush=True)
    STRAT.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as f: f.write(time.strftime("%H:%M:%S ") + msg + "\n")


def step_A_fpocket() -> list[dict]:
    """Run FPocket; return list of pockets ranked by druggability, excluding active/cofactor shells.
    Falls back to manual surface-residue centroids if fpocket fails (arm64-darwin Qhull bug)."""
    log("\n=== A — FPocket on apo receptor ===")
    apo_pdb = STRAT / "apo_for_fpocket.pdb"
    if not apo_pdb.exists():
        lines = []
        for ln in RECEPTOR_APO.read_text().splitlines():
            if ln.startswith(("ATOM","HETATM")): lines.append(ln[:66])
            elif ln.startswith(("END","TER")): lines.append(ln[:66])
        apo_pdb.write_text("\n".join(lines) + "\n")
    fp_out = STRAT / "apo_for_fpocket_out"
    if fp_out.exists(): shutil.rmtree(fp_out)
    if FPOCKET:
        proc = subprocess.run([FPOCKET, "-f", str(apo_pdb)], capture_output=True, timeout=300)
        log(f"  fpocket rc={proc.returncode}")
    info_file = STRAT / "apo_for_fpocket_out" / "apo_for_fpocket_info.txt"
    if not info_file.exists():
        log("  ! fpocket failed (arm64 Qhull/Voronoi bug per roadmap S3) — falling back to manual chain-A surface centroids")
        return _manual_allosteric_candidates(apo_pdb)
    pockets = []
    cur = None
    for ln in info_file.read_text().splitlines():
        if ln.startswith("Pocket "):
            if cur: pockets.append(cur)
            cur = {"cavity_id": ln.split()[1].rstrip(":")}
        elif cur and ":" in ln:
            k, v = ln.split(":", 1)
            cur[k.strip()] = v.strip()
    if cur: pockets.append(cur)
    log(f"  found {len(pockets)} pockets")
    # for each pocket, parse the pocket atoms PDB to get centroid
    for p in pockets:
        pdb = STRAT / "apo_for_fpocket_out" / "pockets" / f"pocket{p['cavity_id']}_atm.pdb"
        if not pdb.exists(): continue
        coords = []
        for ln in pdb.read_text().splitlines():
            if ln.startswith(("ATOM","HETATM")):
                coords.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
        if coords:
            c = np.array(coords).mean(axis=0)
            p["centre"] = c.tolist()
            p["d_active_site"] = float(np.linalg.norm(c - ACTIVE_SITE_CENTRE))
    # Rank by druggability score, filter cavities within 8 Å of active site
    def druggability(p):
        try: return float(p.get("Druggability Score", 0))
        except Exception: return 0.0
    ranked = [p for p in pockets if p.get("d_active_site", 99) > 8.0]
    ranked.sort(key=druggability, reverse=True)
    keep = ranked[:5]
    log(f"  kept {len(keep)} allosteric cavities (≥8 Å from active site, top 5 by druggability)")
    for p in keep:
        log(f"    {p['cavity_id']:>4s}  drug={druggability(p):.3f}  d_site={p.get('d_active_site','?')}  centre={p.get('centre')}")
    (STRAT / "allosteric_cavities.json").write_text(json.dumps(keep, indent=2, default=str))
    return keep


def step_B_fragments(target_n: int = 25) -> list[dict]:
    log("\n=== B — assemble fragment library ===")
    # try ZINC15 first
    fragments = []
    try:
        url = "https://zinc15.docking.org/substances/subsets/fragment/?count=50&format=smi"
        req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=10) as r:
            data = r.read().decode()
        for line in data.strip().splitlines()[:target_n]:
            parts = line.split()
            if len(parts) >= 2:
                fragments.append({"name": parts[1], "cid": None, "smiles": parts[0], "source": "ZINC15"})
        log(f"  fetched {len(fragments)} ZINC15 fragments")
    except Exception as e:
        log(f"  ZINC15 unreachable: {e} — falling back to PubChem fragment list")
    # fallback: PubChem fragment-like
    if len(fragments) < target_n:
        for cid in FALLBACK_FRAGMENTS:
            if len(fragments) >= target_n: break
            try:
                smi = fetch_smiles(cid)
                fragments.append({"name": f"frag_CID{cid}", "cid": cid, "smiles": smi, "source": "PubChem_fallback"})
            except Exception:
                continue
            time.sleep(0.2)
    # descriptors
    for f in fragments: f.update(compound_descriptors(f["smiles"]))
    # filter: drop non-fragment compounds (MW > 350 or n_heavy > 25)
    fragments = [f for f in fragments if f.get("mw", 999) <= 350 and f.get("n_heavy", 99) <= 25]
    log(f"  → {len(fragments)} fragments after MW/heavy-atom filter")
    (STRAT / "fragments.json").write_text(json.dumps(fragments, indent=2))
    return fragments


def step_C_prep(fragments: list[dict]):
    log("\n=== C — fragment prep ===")
    for f in fragments:
        sdf = STRAT / "ligands" / f"{f['name']}.sdf"
        pdbqt = STRAT / "ligands" / f"{f['name']}.pdbqt"
        if not sdf.exists():
            if f.get("cid"):
                try: fetch_sdf_3d(f["cid"], sdf)
                except Exception as e:
                    log(f"  ✗ {f['name']} SDF fetch failed: {e}"); f["prep_ok"]=False; continue
            else:
                # ZINC15 entries — build 3D via RDKit from SMILES
                from rdkit import Chem
                from rdkit.Chem import AllChem
                try:
                    m = Chem.MolFromSmiles(f["smiles"])
                    m = Chem.AddHs(m)
                    AllChem.EmbedMolecule(m, randomSeed=42)
                    AllChem.MMFFOptimizeMolecule(m, maxIters=500)
                    Chem.MolToMolFile(m, str(sdf))
                except Exception as e:
                    log(f"  ✗ {f['name']} RDKit embed failed: {e}"); f["prep_ok"]=False; continue
        if not pdbqt.exists():
            r = prep_ligand_from_sdf(sdf, pdbqt)
            f["prep_ok"] = r["ok"]; f["prep_stage"] = r.get("stage")
            log(f"  {'✓' if r['ok'] else '✗'} {f['name']:25s} → {pdbqt.name} ({r.get('stage')})")
        else:
            f["prep_ok"] = True
    return fragments


def step_D_dock(fragments, cavities, seeds=(42,)):
    log(f"\n=== D — dock fragments into top allosteric cavities ===")
    rows = []
    total = sum(1 for f in fragments if f.get("prep_ok")) * len(cavities) * len(seeds)
    done = 0
    for cav in cavities:
        if not cav.get("centre"): continue
        box = {"cx": cav["centre"][0], "cy": cav["centre"][1], "cz": cav["centre"][2],
               "sx": 20.0, "sy": 20.0, "sz": 20.0}
        log(f"  --- cavity {cav['cavity_id']} centre={cav['centre']} ---")
        for f in fragments:
            if not f.get("prep_ok"): continue
            ligand = STRAT / "ligands" / f"{f['name']}.pdbqt"
            for seed in seeds:
                out = STRAT / "docked" / f"cav{cav['cavity_id']}_{f['name']}_seed{seed}.pdbqt"
                r = vina_dock(RECEPTOR_APO, ligand, out, box, seed=seed, exhaustiveness=16, num_modes=9, cpu=4)
                done += 1
                row = {"fragment": f["name"], "cid": f.get("cid"), "cavity_id": cav["cavity_id"],
                       "cavity_centre_x": cav["centre"][0], "cavity_centre_y": cav["centre"][1], "cavity_centre_z": cav["centre"][2],
                       "cavity_druggability": cav.get("Druggability Score"),
                       "d_active_site": cav.get("d_active_site"),
                       "seed": seed, "top1": r.get("top1"), "top3_mean": r.get("top3_mean"),
                       "ok": r.get("ok"), "wall_s": r.get("wall_s"),
                       "mw": f.get("mw"), "logp": f.get("logp")}
                rows.append(row)
                if done <= 3 or done % 10 == 0:
                    log(f"    [{done}/{total}] cav{cav['cavity_id']} {f['name']:22s} top1={row.get('top1')} {row.get('wall_s')}s")
    csv_path = STRAT / "results_raw.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    log(f"  → {csv_path} ({len(rows)} rows)")
    return rows


def step_G(rows):
    log("\n=== G — aggregate (absolute Vina ranking; no Δ ref for allosteric) ===")
    rows.sort(key=lambda r: (r["cavity_id"], r.get("top1") or 0))
    summary = []
    # per (cavity, fragment) best
    import collections
    by = collections.defaultdict(list)
    for r in rows:
        if r.get("top1") is None: continue
        by[(r["cavity_id"], r["fragment"])].append(r)
    for (cav, frag), rs in by.items():
        best = min(rs, key=lambda x: x["top1"])
        summary.append(best)
    summary.sort(key=lambda r: r.get("top1") or 0)
    csv_path = STRAT / "results_summary.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
    log(f"  Top 10 by absolute Vina score across all cavities:")
    for s in summary[:10]:
        log(f"    cav{s['cavity_id']:>3s}  {s['fragment']:25s}  top1={s['top1']:+.2f}  d_site={s.get('d_active_site')}")


def main():
    STRAT.mkdir(parents=True, exist_ok=True)
    LOG.write_text(f"Phase 14 Strategy 4 — started {time.ctime()}\n")
    t0 = time.time()
    try:
        cavities = step_A_fpocket()
        if not cavities:
            log("!!! No allosteric cavities found — strategy 4 NULL RESULT per Stop Condition S1"); return
        fragments = step_B_fragments()
        fragments = step_C_prep(fragments)
        rows = step_D_dock(fragments, cavities)
        step_G(rows)
        log(f"=== DONE in {(time.time()-t0)/60:.1f} min ===")
    except Exception as e:
        log(f"!!! ABORT: {e}"); log(traceback.format_exc()); raise

if __name__ == "__main__":
    main()
