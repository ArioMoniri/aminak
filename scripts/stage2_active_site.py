#!/usr/bin/env python3
"""Stage 2: Merge UniProt + PDBe active-site/binding annotations, cross with conservation."""
import os, sys, json, time
from datetime import datetime
import requests
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

PROJECT = os.path.expanduser("~/conserved_site_project")
AS_DIR = os.path.join(PROJECT, "02_active_site")
MSA_DIR = os.path.join(PROJECT, "01_msa")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
REF = "P04818"
PDB = "1hvy"

def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] STAGE2: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

def fetch_uniprot_features(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.json()

def fetch_pdbe_ligand_monomers(pdb):
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb}"
    r = requests.get(url, timeout=30); r.raise_for_status(); return r.json()

def fetch_pdbe_bm_interactions(pdb, bm_id):
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/pdb/bound_molecule_interactions/{pdb}/{bm_id}"
    r = requests.get(url, timeout=30); r.raise_for_status(); return r.json()

def main():
    os.makedirs(AS_DIR, exist_ok=True)
    log("Stage 2 starting")

    # UniProt features
    log("Fetching UniProt P04818 features")
    uj = fetch_uniprot_features(REF)
    uniprot_residues = {}  # pos -> (residue_letter, [notes])
    seq = uj.get("sequence", {}).get("value", "")
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
    log(f"UniProt: {len(uniprot_residues)} unique annotated positions")

    # PDBe binding via bound_molecule_interactions
    log("Fetching PDBe ligand monomers for 1HVY")
    monomers = fetch_pdbe_ligand_monomers(PDB).get(PDB, [])
    target_ligands = {"UMP", "D16"}
    AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
           "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
           "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
    pdbe_residues = {}
    # Find chain-A bm_ids for UMP/D16
    bm_targets = []
    for m in monomers:
        if m.get("chain_id") == "A" and m.get("chem_comp_id") in target_ligands:
            bm_targets.append((m["bm_id"], m["chem_comp_id"]))
    log(f"Chain-A target ligands: {bm_targets}")
    for bm_id, lig in bm_targets:
        try:
            data = fetch_pdbe_bm_interactions(PDB, bm_id).get(PDB, [])
        except Exception as e:
            log(f"failed bm {bm_id}: {e}")
            continue
        for entry in data:
            for inter in entry.get("interactions", []):
                # we want the partner that is a protein residue (not the ligand itself)
                for end in (inter.get("begin"), inter.get("end")):
                    if not end:
                        continue
                    rn = end.get("chem_comp_id", "")
                    if rn in AA3 and end.get("chain_id") == "A":
                        ri = end.get("author_residue_number")
                        if ri is None:
                            continue
                        aa = AA3[rn]
                        pdbe_residues.setdefault(ri, (aa, []))[1].append(f"PDBe binding [{lig}]")
    log(f"PDBe: {len(pdbe_residues)} unique residues binding UMP/D16")

    # Note: PDB 1HVY chain A residue numbering follows the human TS sequence
    # (no offset for human TS — confirmed by literature: Cys195, His196, etc.)
    # We treat PDBe positions as already aligned to UniProt numbering.

    # Combine
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
    log(f"Wrote active_site_residues.csv with {len(df)} residues")

    # Cross with conservation
    cons = pd.read_csv(os.path.join(MSA_DIR, "conservation_scores.csv"))
    merged = df.merge(cons, on="ref_position", how="left", suffixes=("","_cons"))
    merged["js_score"] = merged["js_score"].fillna(0)
    merged["percentile"] = merged["percentile"].fillna(0)
    merged.to_csv(os.path.join(AS_DIR, "overlap_table.csv"), index=False)
    log(f"Wrote overlap_table.csv")

    # Venn figure: DB-annotated set vs top-25% conserved
    db_set = set(df["ref_position"].tolist())
    top25 = set(cons[cons["percentile"] >= 75]["ref_position"].tolist())
    fig, ax = plt.subplots(figsize=(6, 6))
    venn2([db_set, top25], set_labels=("DB-annotated\n(UniProt+PDBe)", "Top-25%\nconserved"))
    ax.set_title("Active site DB annotations vs conservation")
    plt.tight_layout()
    plt.savefig(os.path.join(AS_DIR, "overlap_figure.png"), dpi=120)
    plt.close()
    log("Wrote overlap_figure.png")

    # Final selected set: DB-annotated AND top-25% conserved
    selected = sorted(db_set & top25)
    threshold_used = 25
    if len(selected) < 5:
        top40 = set(cons[cons["percentile"] >= 60]["ref_position"].tolist())
        selected = sorted(db_set & top40)
        threshold_used = 40
        log(f"Loosened to top-40% conserved (got {len(selected)})")
    # Augment with literature-known catalytic residues to ensure they're tested
    # even if global conservation pseudocount masks them (TYMS is so conserved
    # that most positions cluster tightly).
    catalytic_must_include = [195, 196]  # Cys195 nucleophile, His196
    for p in catalytic_must_include:
        if p in db_set and p not in selected:
            selected.append(p)
            log(f"Augmented with catalytic residue {p}")
    selected = sorted(selected)
    # Cap at 8 for runtime budget
    if len(selected) > 8:
        # keep highest conservation 6 + the 2 catalytic
        cons_lookup = {p: float(cons[cons.ref_position==p].js_score.iloc[0]) for p in selected}
        keep = set(catalytic_must_include) & set(selected)
        remaining = sorted([p for p in selected if p not in keep],
                           key=lambda x: -cons_lookup[x])[:8 - len(keep)]
        selected = sorted(keep | set(remaining))
        log(f"Capped to 8: {selected}")
    log(f"Final selected residues ({len(selected)}, threshold top-{threshold_used}% +catalytic): {selected}")

    sel_df = merged[merged["ref_position"].isin(selected)].copy()
    sel_df["selected"] = True
    sel_df.to_csv(os.path.join(AS_DIR, "selected_residues.csv"), index=False)

    with open(os.path.join(AS_DIR, "selected_meta.json"), "w") as f:
        json.dump({"selected": selected, "threshold_percentile": threshold_used,
                   "n_db": len(db_set), "n_topcons": len(top25)}, f, indent=2)

    log("Stage 2 DONE")

if __name__ == "__main__":
    main()
