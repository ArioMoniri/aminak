#!/usr/bin/env python3
"""Stage 2 v2: Active-site selection from UniProt+PDBe ∩ top-conserved.
No force-augment of catalytic residues — they should now naturally rank in top set.
"""
import os, sys, json
from datetime import datetime
import requests
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

PROJECT = os.path.expanduser("~/conserved_site_project")
AS_DIR = os.path.join(PROJECT, "02b_active_site_v2")
MSA_DIR = os.path.join(PROJECT, "01b_msa_v2")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v2_02_active_site.log")
REF = "P04818"
PDB = "1hvy"

AA3to1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
          "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
          "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V2] STAGE2: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def fetch_uniprot_features(acc):
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=30)
    r.raise_for_status()
    return r.json()


def fetch_pdbe_ligand_monomers(pdb):
    r = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb}", timeout=30)
    r.raise_for_status()
    return r.json()


def fetch_pdbe_bm_interactions(pdb, bm_id):
    r = requests.get(f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/bound_molecule_interactions/{pdb}/{bm_id}", timeout=30)
    r.raise_for_status()
    return r.json()


def main():
    os.makedirs(AS_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 2 v2 starting")

    # UniProt features
    uj = fetch_uniprot_features(REF)
    seq = uj.get("sequence", {}).get("value", "")
    uniprot_residues = {}
    for feat in uj.get("features", []):
        ftype = feat.get("type", "")
        if ftype not in ("Active site", "Binding site", "Site"):
            continue
        loc = feat.get("location", {})
        s = loc.get("start", {}).get("value")
        e = loc.get("end", {}).get("value")
        desc = feat.get("description", "")
        ligand = feat.get("ligand", {}).get("name", "") if feat.get("ligand") else ""
        if s is None:
            continue
        for pos in range(s, (e or s) + 1):
            res = seq[pos-1] if 0 < pos <= len(seq) else "?"
            note = f"{ftype}: {desc}" + (f" [{ligand}]" if ligand else "")
            uniprot_residues.setdefault(pos, (res, []))[1].append(note)
    log(f"UniProt: {len(uniprot_residues)} annotated positions")

    # PDBe binding via bound_molecule_interactions for both UMP and D16, both chains A and B
    monomers = fetch_pdbe_ligand_monomers(PDB).get(PDB, [])
    target_ligands = {"UMP", "D16"}
    pdbe_residues = {}
    bm_targets = []
    for m in monomers:
        if m.get("chain_id") in ("A", "B") and m.get("chem_comp_id") in target_ligands:
            bm_targets.append((m["bm_id"], m["chem_comp_id"], m["chain_id"]))
    log(f"chain-A/B target ligands: {bm_targets}")
    for bm_id, lig, chid in bm_targets:
        try:
            data = fetch_pdbe_bm_interactions(PDB, bm_id).get(PDB, [])
        except Exception as e:
            log(f"  failed bm {bm_id}: {e}")
            continue
        for entry in data:
            for inter in entry.get("interactions", []):
                for end in (inter.get("begin"), inter.get("end")):
                    if not end:
                        continue
                    rn = end.get("chem_comp_id", "")
                    if rn in AA3to1 and end.get("chain_id") in ("A", "B"):
                        ri = end.get("author_residue_number")
                        if ri is None:
                            continue
                        aa = AA3to1[rn]
                        pdbe_residues.setdefault(ri, (aa, []))[1].append(f"PDBe binding [{lig}/{chid}]")
    log(f"PDBe: {len(pdbe_residues)} unique residues binding UMP/D16")

    all_pos = sorted(set(uniprot_residues.keys()) | set(pdbe_residues.keys()))
    rows = []
    for p in all_pos:
        if p in uniprot_residues and p in pdbe_residues:
            src = "both"
            res = uniprot_residues[p][0]
            notes = "; ".join(uniprot_residues[p][1] + pdbe_residues[p][1])
        elif p in uniprot_residues:
            src = "uniprot"
            res = uniprot_residues[p][0]
            notes = "; ".join(uniprot_residues[p][1])
        else:
            src = "pdbe"
            res = pdbe_residues[p][0]
            notes = "; ".join(pdbe_residues[p][1])
        rows.append({"ref_position": p, "residue": res, "source": src, "notes": notes})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(AS_DIR, "active_site_residues.csv"), index=False)
    log(f"Wrote active_site_residues.csv ({len(df)} rows)")

    # Cross with v2 conservation
    cons = pd.read_csv(os.path.join(MSA_DIR, "conservation_scores.csv"))
    merged = df.merge(cons, on="ref_position", how="left", suffixes=("", "_cons"))
    merged["js_score"] = merged["js_score"].fillna(0)
    merged["percentile"] = merged["percentile"].fillna(0)
    merged.to_csv(os.path.join(AS_DIR, "overlap_table.csv"), index=False)
    log("Wrote overlap_table.csv")

    db_set = set(df["ref_position"].tolist())
    top25 = set(cons[cons["percentile"] >= 75]["ref_position"].dropna().tolist())
    top10 = set(cons[cons["percentile"] >= 90]["ref_position"].dropna().tolist())
    log(f"DB-annotated set: {len(db_set)}; top25%: {len(top25)}; top10%: {len(top10)}")
    log(f"DB ∩ top25%: {len(db_set & top25)}; DB ∩ top10%: {len(db_set & top10)}")

    fig, ax = plt.subplots(figsize=(6, 6))
    venn2([db_set, top25], set_labels=("DB-annotated\n(UniProt+PDBe)", "Top-25%\nconserved"))
    ax.set_title("v2: Active site DB annotations vs conservation")
    plt.tight_layout()
    plt.savefig(os.path.join(AS_DIR, "overlap_figure.png"), dpi=130)
    plt.close()

    selected = sorted(db_set & top25)
    threshold_used = 25
    if len(selected) < 6:
        selected = sorted(db_set & set(cons[cons["percentile"] >= 60]["ref_position"].dropna().tolist()))
        threshold_used = 40
        log(f"Loosened to top-40% (got {len(selected)})")

    # Cap at 10 for downstream mutation panel coverage
    if len(selected) > 10:
        cons_lookup = {p: float(cons[cons.ref_position == p].js_score.iloc[0]) for p in selected}
        selected = sorted(selected, key=lambda x: -cons_lookup[x])[:10]
        selected = sorted(selected)
        log(f"Capped to 10 highest-cons: {selected}")
    log(f"Final selected ({len(selected)}, threshold top-{threshold_used}%): {selected}")

    sel_df = merged[merged["ref_position"].isin(selected)].copy()
    sel_df["selected"] = True
    sel_df.to_csv(os.path.join(AS_DIR, "selected_residues.csv"), index=False)

    with open(os.path.join(AS_DIR, "selected_meta.json"), "w") as f:
        json.dump({"selected": selected, "threshold_percentile": threshold_used,
                   "n_db": len(db_set), "n_topcons_25": len(top25), "n_topcons_10": len(top10),
                   "force_augmented": []}, f, indent=2)
    log("Stage 2 v2 DONE")


if __name__ == "__main__":
    main()
