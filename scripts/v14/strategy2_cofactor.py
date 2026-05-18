#!/usr/bin/env python3
"""Phase 14 — Strategy 2 (cofactor-site, antifolates) — SCOPED execution.

Box centre computed from holo cofactor A coords (D16).
Anchors: methotrexate, raltitrexed, pemetrexed, nolatrexed, plevitrexed, ibuprofen (neg control).
"""
from __future__ import annotations
import json, sys, csv, time, traceback
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (
    REPO, PHASE, RECEPTOR_APO, COFACTOR_A, COFACTOR_B,
    pubchem_get, fetch_sdf_3d, fetch_smiles, prep_ligand_from_sdf,
    vina_dock, rdkit_decoys, compound_descriptors, load_anchors,
)

STRAT = PHASE / "02_cofactor_site"
LOG = STRAT / "strategy2_runlog.txt"


def log(msg: str):
    print(msg, flush=True)
    STRAT.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as f: f.write(time.strftime("%H:%M:%S ") + msg + "\n")


def compute_cofactor_box() -> dict:
    """Centroid of cofactor A heavy atoms from the holo PDBQT."""
    coords = []
    for line in COFACTOR_A.read_text().splitlines():
        if line.startswith(("ATOM","HETATM")):
            elem = line[76:78].strip() if len(line) >= 78 else ""
            if elem.startswith("H"): continue
            try:
                x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append((x,y,z))
            except ValueError: continue
    if not coords: raise RuntimeError("no cofactor coords parsed")
    c = np.array(coords).mean(axis=0)
    box = {"cx": float(c[0]), "cy": float(c[1]), "cz": float(c[2]),
           "sx": 22.0, "sy": 22.0, "sz": 22.0}
    log(f"  cofactor box centre (from holo D16:A heavy-atom centroid) = ({c[0]:+.3f}, {c[1]:+.3f}, {c[2]:+.3f})")
    return box


def receptor_path(state: str) -> Path:
    if state == "apo": return RECEPTOR_APO
    holo = STRAT / "receptor_holo.pdbqt"
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


def step_B(box: dict) -> list[dict]:
    log("\n=== B — assemble cofactor-site compound set ===")
    cofactor_names = ["Methotrexate", "Raltitrexed", "Pemetrexed", "Nolatrexed", "Plevitrexed", "Ibuprofen"]
    anchors = [a for a in load_anchors()
               if any(a["display_name"].startswith(n) for n in cofactor_names)]
    log(f"  Strategy 2 Tier-1 anchors: {len(anchors)}")
    compounds = []
    for a in anchors:
        smi = a.get("smiles") or fetch_smiles(a["cid"])
        compounds.append({"name": a["display_name"].split("(")[0].strip().split("/")[0].strip().replace(" ", "_"),
                          "cid": a["cid"], "smiles": smi, "tier": 1, "anchor_for": a["display_name"]})
    # Tier 2: RDKit decoys against raltitrexed (the cofactor-site reference)
    ralt = next((a for a in anchors if a["display_name"].startswith("Raltitrexed")), None)
    if ralt:
        log("  Generating RDKit decoys against raltitrexed …")
        ralt_smi = ralt.get("smiles") or fetch_smiles(ralt["cid"])
        decoys = rdkit_decoys(ralt_smi, n=25)
        log(f"  Got {len(decoys)} decoys")
        for d in decoys:
            compounds.append({"name": f"decoy_CID{d['cid']}", "cid": d["cid"],
                              "smiles": d["smiles"], "tier": 2, "anchor_for": "raltitrexed_decoy"})
    for c in compounds:
        c.update(compound_descriptors(c["smiles"]))
    (STRAT / "compounds.json").write_text(json.dumps(compounds, indent=2))
    return compounds


def step_C(compounds: list[dict]) -> list[dict]:
    log("\n=== C — ligand prep ===")
    for c in compounds:
        sdf = STRAT / "ligands" / f"{c['name']}.sdf"
        pdbqt = STRAT / "ligands" / f"{c['name']}.pdbqt"
        if not sdf.exists():
            try: fetch_sdf_3d(c["cid"], sdf)
            except Exception as e:
                log(f"  ✗ {c['name']} SDF fetch failed: {e}"); c["prep_ok"]=False; continue
        if not pdbqt.exists():
            r = prep_ligand_from_sdf(sdf, pdbqt)
            c["prep_ok"] = r["ok"]; c["prep_stage"] = r.get("stage")
            log(f"  {'✓' if r['ok'] else '✗'} {c['name']} → {pdbqt.name} ({r.get('stage')})")
        else:
            c["prep_ok"] = True; c["prep_stage"] = "cached"
    (STRAT / "compounds_prepped.json").write_text(json.dumps(compounds, indent=2))
    return compounds


def step_A0(compounds: list[dict], box: dict):
    log("\n=== A0 — re-dock raltitrexed into the holo cofactor box (sanity) ===")
    ralt = next((c for c in compounds if c["name"].lower().startswith("raltit")), None)
    if not ralt or not ralt.get("prep_ok"):
        log("  ! raltitrexed prep failed; A0 skipped"); return {"ok": False}
    ligand = STRAT / "ligands" / f"{ralt['name']}.pdbqt"
    out = STRAT / "A0_redock_gate" / "raltitrexed_redock.pdbqt"
    r = vina_dock(receptor_path("apo"), ligand, out, box, seed=42, exhaustiveness=32)
    log(f"  raltitrexed top1={r.get('top1')} (apo cofactor box)")
    (STRAT / "A0_redock_gate").mkdir(parents=True, exist_ok=True)
    (STRAT / "A0_redock_gate" / "result.json").write_text(json.dumps(r, indent=2))
    return r


def step_D(compounds, box, seeds=(42,7), states=("apo","holo")) -> list[dict]:
    log(f"\n=== D — Vina dock ({len(compounds)} × {len(seeds)} × {len(states)} = {len(compounds)*len(seeds)*len(states)} runs) ===")
    rows = []; done=0
    total = sum(1 for c in compounds if c.get("prep_ok")) * len(seeds) * len(states)
    for c in compounds:
        if not c.get("prep_ok"): continue
        ligand = STRAT / "ligands" / f"{c['name']}.pdbqt"
        for state in states:
            rec = receptor_path(state)
            for seed in seeds:
                out_pdbqt = STRAT / "docked" / f"{c['name']}_{state}_seed{seed}.pdbqt"
                if out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
                    log_path = out_pdbqt.with_suffix(".log")
                    affs = []
                    if log_path.exists():
                        for ln in log_path.read_text().splitlines():
                            p = ln.split()
                            if len(p)>=4 and p[0].isdigit():
                                try: affs.append(float(p[1]))
                                except ValueError: pass
                    res = {"ok": True, "cached": True,
                           "top1": affs[0] if affs else None,
                           "top3_mean": (sum(affs[:3])/min(3,len(affs))) if affs else None,
                           "top5_mean": (sum(affs[:5])/min(5,len(affs))) if affs else None,
                           "n_modes": len(affs)}
                else:
                    res = vina_dock(rec, ligand, out_pdbqt, box, seed=seed, exhaustiveness=16, cpu=4, num_modes=10)
                done += 1
                row = {"compound": c["name"], "cid": c.get("cid"), "tier": c.get("tier"),
                       "anchor_for": c.get("anchor_for"), "state": state, "seed": seed,
                       "top1": res.get("top1"), "top3_mean": res.get("top3_mean"), "top5_mean": res.get("top5_mean"),
                       "n_modes": res.get("n_modes"), "ok": res.get("ok"), "wall_s": res.get("wall_s"),
                       "mw": c.get("mw"), "logp": c.get("logp"),
                       "pains": c.get("pains_a") or c.get("pains_b") or c.get("pains_c"),
                       "brenk": c.get("brenk"), "lipinski_ok": c.get("lipinski_ok")}
                rows.append(row)
                log(f"  [{done}/{total}] {c['name']:30s} {state:5s} s={seed} top1={row.get('top1')} {row.get('wall_s')}s")
    csv_path = STRAT / "results_raw.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    return rows


def step_G(rows):
    log("\n=== G — aggregate (Δ vs raltitrexed) ===")
    import collections
    by = collections.defaultdict(list)
    for r in rows:
        if r.get("top1") is None: continue
        by[(r["compound"], r["state"])].append(r["top1"])
    summary = []
    for (cmpd, state), ts in by.items():
        meta = next(r for r in rows if r["compound"] == cmpd and r["state"] == state)
        mean = sum(ts)/len(ts)
        s = {"compound": cmpd, "state": state, "tier": meta["tier"], "anchor_for": meta["anchor_for"],
             "cid": meta["cid"], "n_seeds": len(ts), "top1_mean": float(mean),
             "top1_min": float(min(ts)), "top1_max": float(max(ts)),
             "top1_sd": float((sum((t-mean)**2 for t in ts)/len(ts))**0.5),
             "mw": meta.get("mw"), "logp": meta.get("logp"),
             "pains": meta.get("pains"), "brenk": meta.get("brenk"), "lipinski_ok": meta.get("lipinski_ok")}
        summary.append(s)
    ralt_by_state = {s["state"]: s["top1_mean"] for s in summary if s["compound"].lower().startswith("raltit")}
    for s in summary:
        ref = ralt_by_state.get(s["state"])
        s["delta_vs_raltitrexed"] = (s["top1_mean"] - ref) if ref is not None else None
        s["beats_raltitrexed"] = (s["delta_vs_raltitrexed"] is not None and s["delta_vs_raltitrexed"] < -0.85)
    csv_path = STRAT / "results_summary.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
    log(f"  → {csv_path}")
    summary.sort(key=lambda x: (x["state"], x["top1_mean"]))
    for state in ("apo","holo"):
        log(f"  --- {state} ---")
        for s in [r for r in summary if r["state"] == state][:10]:
            log(f"    {s['compound']:25s} top1={s['top1_mean']:+.2f}  Δralt={s['delta_vs_raltitrexed']}  tier={s['tier']}")


def main():
    STRAT.mkdir(parents=True, exist_ok=True)
    LOG.write_text(f"Phase 14 Strategy 2 — started {time.ctime()}\n")
    t0 = time.time()
    try:
        box = compute_cofactor_box()
        (STRAT / "box.json").write_text(json.dumps(box, indent=2))
        compounds = step_B(box)
        compounds = step_C(compounds)
        step_A0(compounds, box)
        rows = step_D(compounds, box)
        step_G(rows)
        log(f"=== DONE in {(time.time()-t0)/60:.1f} min ===")
    except Exception as e:
        log(f"!!! ABORT: {e}"); log(traceback.format_exc()); raise

if __name__ == "__main__":
    main()
