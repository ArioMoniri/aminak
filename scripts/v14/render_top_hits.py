#!/usr/bin/env python3
"""Phase 14 — render PyMOL images + interaction analysis for the Strategy-4 top hits.

For each top hit:
  1. Compute heavy-atom contacts within 4 Å between ligand and protein (residue-level).
  2. Classify each contact (hydrophobic / H-bond donor-acceptor / aromatic π-stacking / salt-bridge)
     using a simple element + distance + functional-group filter — no PROLIF dependency needed.
  3. Render a PyMOL ray-traced PNG: cartoon protein (chain A grey, chain B lightgrey) + cavity
     residues as yellow sticks + ligand as cyan sticks + dashed yellow lines for polar contacts.

Outputs:
  14_inhibitor_design/04_allosteric/poses/<hit_id>.png
  14_inhibitor_design/04_allosteric/poses/<hit_id>_interactions.json
  14_inhibitor_design/04_allosteric/poses/all_interactions.csv  (one row per ligand-residue contact)
"""
from __future__ import annotations
import json, csv, subprocess, shutil, tempfile
from pathlib import Path
from collections import defaultdict
import numpy as np

REPO = Path(__file__).resolve().parents[2]
STRAT = REPO / "14_inhibitor_design" / "04_allosteric"
OUT = STRAT / "poses"; OUT.mkdir(parents=True, exist_ok=True)
APO_PDB = STRAT / "apo_for_fpocket.pdb"
PYMOL = shutil.which("pymol")

# Top 5 hits from results_summary.csv
TOP_HITS = [
    {"id": "cav18_CID7032", "compound": "1H-indazole", "cid": 7032, "cavity": 18, "top1": -7.52, "drug": 0.994,
     "use_chain": "B",
     "summary": "Kinase-inhibitor privileged scaffold (e.g. axitinib, niraparib, pazopanib all carry indazole). A small bicyclic 5-6 heteroaromatic; can act as H-bond donor (N1-H) and acceptor (N2)."},
    {"id": "cav18_CID3672", "compound": "ibuprofen", "cid": 3672, "cavity": 18, "top1": -7.28, "drug": 0.994,
     "use_chain": "B",
     "summary": "Non-steroidal anti-inflammatory drug; COX1/2 inhibitor. Carboxylate + isobutyl-phenyl combination is famously promiscuous, binding many off-targets (HSA pocket I, FABP, CRBN)."},
    {"id": "cav2_CID5564",  "compound": "tolnaftate",  "cid": 5564, "cavity": 2,  "top1": -6.88, "drug": 0.009,
     "use_chain": "B",
     "summary": "Topical antifungal (thiocarbamate); not a TYMS literature ligand. Strongly lipophilic (logP ~5.1)."},
    {"id": "cav2_CID7032",  "compound": "1H-indazole", "cid": 7032, "cavity": 2,  "top1": -6.86, "drug": 0.009,
     "use_chain": "B",
     "summary": "Same scaffold as above but docked in the low-druggability cavity 2 — comparison shows the −7.5 kcal/mol signal at cavity 18 tracks the *pocket*, not the *library*."},
    {"id": "cav12_CID35814","compound": "flurbiprofen","cid": 35814,"cavity": 12, "top1": -6.52, "drug": 0.010,
     "use_chain": "B",
     "summary": "NSAID (COX inhibitor), related to ibuprofen but with a fluorobiphenyl. Surface binding only."},
]

# Cavity centroids from results_summary
CAV_CENTRE = {
    18: (4.564,  -12.706, -14.884),
    17: (-4.564, +12.706, +14.884),   # C2 mirror estimate
     2: (12.200, -14.477,  -9.151),
    12: (17.956,   0.959,  -1.468),
}


def parse_atoms(path: Path, model_idx: int = 1):
    """Return list of dicts: {idx, name, resname, chain, resid, x, y, z, element}.
    PDBQT MODEL 1 only."""
    atoms = []
    in_model = (model_idx == 0)
    cur_model = 0
    for ln in path.read_text().splitlines():
        if ln.startswith("MODEL "):
            cur_model = int(ln.split()[1])
            in_model = (cur_model == model_idx) or (model_idx == 0)
            continue
        if ln.startswith("ENDMDL") and in_model and model_idx != 0:
            break
        if not in_model: continue
        if ln.startswith(("ATOM","HETATM")):
            name    = ln[12:16].strip()
            resname = ln[17:20].strip()
            chain   = ln[21].strip() or "A"
            try: resid = int(ln[22:26])
            except ValueError: resid = -1
            try:
                x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            except ValueError: continue
            # element: cols 76-78 if present, else first letter of name
            el = ln[76:78].strip() if len(ln) >= 78 else ""
            if not el:
                el = "".join(c for c in name if c.isalpha())[:1]
            atoms.append({"name": name, "resname": resname, "chain": chain, "resid": resid,
                          "x": x, "y": y, "z": z, "element": el.upper()[:1]})
    return atoms


def contacts(ligand_atoms, protein_atoms, cutoff: float = 4.0):
    """Per-residue interaction summary."""
    by_res = defaultdict(lambda: {"heavy_count": 0, "min_d": 999.0, "hbond": False,
                                  "salt_bridge": False, "hydrophobic": False, "aromatic": False,
                                  "ligand_partners": set(), "atoms": set()})
    for la in ligand_atoms:
        if la["element"] == "H": continue
        for pa in protein_atoms:
            if pa["element"] == "H": continue
            d = ((la["x"]-pa["x"])**2 + (la["y"]-pa["y"])**2 + (la["z"]-pa["z"])**2)**0.5
            if d > cutoff: continue
            key = (pa["chain"], pa["resid"], pa["resname"])
            r = by_res[key]
            r["heavy_count"] += 1
            r["min_d"] = min(r["min_d"], d)
            r["atoms"].add(pa["name"])
            r["ligand_partners"].add(la["name"])
            # very crude interaction classifier:
            # H-bond candidate: O/N on protein <=3.5 Å from O/N on ligand
            if la["element"] in ("O","N") and pa["element"] in ("O","N") and d <= 3.5:
                r["hbond"] = True
            # salt-bridge: anionic side chain (ASP/GLU) <=4.5 Å from a ligand N (or vice versa)
            if pa["resname"] in ("ASP","GLU") and pa["name"] in ("OD1","OD2","OE1","OE2") and la["element"]=="N" and d <= 4.5:
                r["salt_bridge"] = True
            if pa["resname"] in ("LYS","ARG","HIS") and pa["name"] in ("NZ","NH1","NH2","NE2","ND1") and la["element"]=="O" and d <= 4.5:
                r["salt_bridge"] = True
            # aromatic: ring residues (PHE/TYR/TRP/HIS) ring atom <= 5 Å from ligand C
            if pa["resname"] in ("PHE","TYR","TRP","HIS") and pa["name"].startswith(("CD","CE","CG","CZ","NE","CH")) and la["element"]=="C" and d <= 5.0:
                r["aromatic"] = True
            # hydrophobic contact: C-C within 4.5 Å
            if la["element"] == "C" and pa["element"] == "C" and d <= 4.5:
                r["hydrophobic"] = True
    # finalize
    out = []
    for (chain, resid, resname), r in sorted(by_res.items(), key=lambda kv: kv[1]["min_d"]):
        r["chain"] = chain; r["resid"] = resid; r["resname"] = resname
        r["atoms"] = sorted(list(r["atoms"]))
        r["ligand_partners"] = sorted(list(r["ligand_partners"]))
        out.append(r)
    return out


def render_pymol(hit: dict, ligand_pose_pdb: Path, png_out: Path):
    """Render: cartoon protein + cavity-residues sticks + ligand sticks + polar contacts."""
    cx, cy, cz = CAV_CENTRE[hit["cavity"]]
    pml = f"""
reinitialize
load {APO_PDB}, recep
load {ligand_pose_pdb}, lig
# styling
bg_color white
hide everything
show cartoon, recep
color grey80, recep and chain A
color grey50, recep and chain B
set cartoon_transparency, 0.20
# cavity residues — within 6 Å of the ligand
select cav_res, byres (recep within 6.0 of lig)
show sticks, cav_res and not name H*
util.cnc cav_res
color yelloworange, cav_res and elem C
# ligand
show sticks, lig
color cyan, lig and elem C
color red, lig and elem O
color blue, lig and elem N
color yellow, lig and elem S
color magenta, lig and elem F
color forest, lig and elem Cl
# polar contacts (H-bond < 3.5 Å between N/O atoms)
distance hbond, (recep and (elem N+O) and not elem H), (lig and (elem N+O)), 3.5
hide labels, hbond
color yellow, hbond
set dash_color, yellow
set dash_width, 2.5
# label key residues
label cav_res and name CA, "%s/%s%d" % (chain, resn, resi)
set label_size, 11
set label_color, black
# camera
orient lig
zoom lig, 12
ray 1200, 1200
png {png_out}, dpi=140
quit
"""
    pml_path = png_out.with_suffix(".pml")
    pml_path.write_text(pml)
    subprocess.run([PYMOL, "-cq", str(pml_path)], capture_output=True, timeout=120)
    return png_out.exists() and png_out.stat().st_size > 0


def main():
    if not PYMOL:
        print("PyMOL not found"); return
    # Build a chained PDB for the apo (no PDBQT-isms) — keep chains A and B
    if not APO_PDB.exists():
        print(f"!!! {APO_PDB} missing"); return
    protein_atoms = parse_atoms(APO_PDB, model_idx=0)
    # Assign chain id by atom-count split (Phase-6c apo PDB doesn't carry chain ids in cols 22)
    if all(a["chain"] in ("","A","_") for a in protein_atoms):
        n = len(protein_atoms) // 2
        for i, a in enumerate(protein_atoms):
            a["chain"] = "A" if i < n else "B"
    print(f"Protein: {len(protein_atoms)} heavy + H atoms, chains {set(a['chain'] for a in protein_atoms)}")

    all_rows = []
    for hit in TOP_HITS:
        pose_pdbqt = STRAT / "docked" / f"cav{hit['cavity']}_frag_CID{hit['cid']}_seed42.pdbqt"
        if not pose_pdbqt.exists():
            print(f"!!! {pose_pdbqt} missing"); continue
        # extract MODEL 1 to a PDB
        pose_pdb = OUT / f"{hit['id']}_pose1.pdb"
        lines = []
        in_m1 = False
        for ln in pose_pdbqt.read_text().splitlines():
            if ln.startswith("MODEL 1"): in_m1 = True; continue
            if ln.startswith("ENDMDL") and in_m1: break
            if in_m1 and (ln.startswith("ATOM") or ln.startswith("HETATM")):
                lines.append(ln[:66])
        pose_pdb.write_text("\n".join(lines) + "\nEND\n")
        ligand_atoms = parse_atoms(pose_pdb, model_idx=0)
        # filter to model 1 only (already done)
        ints = contacts(ligand_atoms, protein_atoms, cutoff=4.0)
        ints_filt = [r for r in ints if r["heavy_count"] >= 1 and r["chain"] == hit["use_chain"]]
        # save interactions JSON
        json_out = OUT / f"{hit['id']}_interactions.json"
        json_out.write_text(json.dumps({"hit": hit, "n_residue_contacts": len(ints_filt),
                                       "contacts": ints_filt}, indent=2,
                                       default=lambda o: list(o) if isinstance(o, set) else o))
        # accumulate cross-strategy rows
        for r in ints_filt:
            all_rows.append({
                "hit_id": hit["id"], "compound": hit["compound"], "cid": hit["cid"],
                "cavity": hit["cavity"], "top1": hit["top1"], "druggability": hit["drug"],
                "chain": r["chain"], "resid": r["resid"], "resname": r["resname"],
                "min_d": round(r["min_d"], 3), "heavy_count": r["heavy_count"],
                "hbond": r["hbond"], "salt_bridge": r["salt_bridge"],
                "hydrophobic": r["hydrophobic"], "aromatic": r["aromatic"],
                "lig_partners": ";".join(r["ligand_partners"]),
                "prot_atoms": ";".join(r["atoms"]),
            })
        # render
        png_out = OUT / f"{hit['id']}.png"
        ok = render_pymol(hit, pose_pdb, png_out)
        print(f"  {'✓' if ok else '✗'} {hit['id']:18s} top1={hit['top1']:+.2f}  drug={hit['drug']:.3f}  "
              f"n_contacts={len(ints_filt)}  render={png_out.name}")
    csv_out = OUT / "all_interactions.csv"
    if all_rows:
        with csv_out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            w.writeheader(); w.writerows(all_rows)
        print(f"  → {csv_out} ({len(all_rows)} rows)")

if __name__ == "__main__":
    main()
