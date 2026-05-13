#!/usr/bin/env python3
"""Task B: enumerate all single mutations at the 14 active-site residues
and all double mutations restricted to those residues.

Outputs:
  12_phase7/02_enum/all_singles.csv
  12_phase7/02_enum/all_doubles_sample.csv
  12_phase7/02_enum/feasibility_note.md
  12_phase7/02_enum/all_singles_chemistry_map.html
"""
from __future__ import annotations
import csv
import itertools
import json
import os
import sys
from pathlib import Path

import plotly.express as px
import pandas as pd

PROJECT = Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
OUT = PROJECT / "12_phase7" / "02_enum"

# 14 active-site residues. Use the panel that lines the binding pocket.
# Drawn from 02b_active_site_v2/active_site_residues.csv "both" rows + the
# critical PDBe-only residues at distance < 5 A from dUMP.
PANEL = [
    (50, "R"),
    (109, "W"),
    (170, "T"),    # included since user references T170A in mutant panel
    (175, "R"),
    (176, "R"),
    (195, "C"),
    (196, "H"),
    (214, "Q"),
    (215, "R"),
    (216, "S"),
    (225, "F"),
    (226, "N"),
    (256, "H"),
    (258, "Y"),
]
ALL_AA = list("ACDEFGHIKLMNPQRSTVWY")

# Kyte-Doolittle hydropathy
KD = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}
# Side-chain volume (A^3, Zamyatnin)
VOL = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5,
    "E": 138.4, "Q": 143.8, "G": 60.1, "H": 153.2, "I": 166.7,
    "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7,
    "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0,
}

CLASS = {
    "A": "hydrophobic", "V": "hydrophobic", "L": "hydrophobic",
    "I": "hydrophobic", "M": "hydrophobic", "F": "aromatic",
    "W": "aromatic", "Y": "aromatic", "P": "special",
    "G": "special", "C": "polar_thiol",
    "S": "polar", "T": "polar", "N": "polar", "Q": "polar",
    "K": "positive", "R": "positive", "H": "positive",
    "D": "negative", "E": "negative",
}

def classify(wt: str, new: str) -> str:
    cw, cn = CLASS[wt], CLASS[new]
    if cw == cn:
        return "conservative"
    if (cw, cn) in {("positive","negative"), ("negative","positive")}:
        return "charge_reversal"
    if cw in {"positive","negative"} and cn not in {"positive","negative"}:
        return "loss_of_charge"
    if cn in {"positive","negative"} and cw not in {"positive","negative"}:
        return "gain_of_charge"
    if cn == "hydrophobic" and cw in {"polar","positive","negative","polar_thiol"}:
        return "polar_to_hydrophobic"
    if cw == "hydrophobic" and cn in {"polar","positive","negative","polar_thiol"}:
        return "hydrophobic_to_polar"
    if cw == "aromatic" and cn != "aromatic":
        return "loss_of_aromatic"
    if cn == "aromatic" and cw != "aromatic":
        return "gain_of_aromatic"
    return "other"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    # Singles
    singles = []
    for pos, wt in PANEL:
        for new in ALL_AA:
            if new == wt:
                continue
            singles.append({
                "residue_position": pos,
                "wt_aa": wt,
                "new_aa": new,
                "mutation_id": f"{wt}{pos}{new}",
                "hydropathy_change": round(KD[new] - KD[wt], 3),
                "volume_change": round(VOL[new] - VOL[wt], 2),
                "functional_class": classify(wt, new),
            })
    spath = OUT / "all_singles.csv"
    with spath.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(singles[0].keys()))
        w.writeheader()
        w.writerows(singles)
    print(f"singles written: {len(singles)} -> {spath}")

    # Doubles restricted to panel pairs
    doubles = []
    for (p1, wt1), (p2, wt2) in itertools.combinations(PANEL, 2):
        for n1 in ALL_AA:
            if n1 == wt1: continue
            for n2 in ALL_AA:
                if n2 == wt2: continue
                doubles.append({
                    "pos1": p1, "wt1": wt1, "new1": n1,
                    "pos2": p2, "wt2": wt2, "new2": n2,
                    "mutation_id": f"{wt1}{p1}{n1}_{wt2}{p2}{n2}",
                    "sum_hydropathy_change": round((KD[n1]-KD[wt1]) + (KD[n2]-KD[wt2]), 3),
                    "sum_volume_change": round((VOL[n1]-VOL[wt1]) + (VOL[n2]-VOL[wt2]), 2),
                    "class1": classify(wt1, n1),
                    "class2": classify(wt2, n2),
                })
    dpath = OUT / "all_doubles_sample.csv"
    with dpath.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(doubles[0].keys()))
        w.writeheader()
        w.writerows(doubles)
    print(f"doubles written: {len(doubles)} -> {dpath}")

    # Feasibility note
    n_pos = len(PANEL)
    n_singles = n_pos * (len(ALL_AA) - 1)
    n_doubles = (n_pos * (n_pos - 1) // 2) * ((len(ALL_AA) - 1) ** 2)
    sec_per = 5
    note = f"""# Feasibility note: enumeration of active-site mutations

We restrict the enumeration to the **{n_pos} residues** that line the dUMP /
folate active site of human TYMS (positions: {", ".join(str(p) for p,_ in PANEL)}).

## Singles
- Count: {n_pos} positions x 19 alternative AAs = **{n_singles} singles**.
- Per-mutation cost: build mutant model (PyMOL mutagenesis), prep PDBQT
  (PDB2PQR + Meeko), Vina dock 5x replicates @ exh 32 num_modes 20.
- Realistic wall-time on this 4-CPU laptop: ~5-10 s per single dock.
  Single-seed sweep: {n_singles} x 5 s ~= {n_singles*sec_per} s ({n_singles*sec_per//60} min).
  5-replicate sweep: {n_singles*sec_per*5} s (~{n_singles*sec_per*5//60} min).
  Plus PyMOL mutagenesis + protonation pre-step: ~10 s/mutant -> ~{n_singles*10//60} min.
  Realistic: 1-2 hours for the full singles sweep with replicates.

## Doubles
- Count restricted to *pairs from this panel only*:
  C({n_pos},2) x 19^2 = {n_pos*(n_pos-1)//2} x 361 = **{n_doubles} doubles**.
- At 5 s per single dock that is ~{n_doubles*sec_per/3600:.1f} h for one seed each.
  With 5 replicates that becomes ~{n_doubles*sec_per*5/3600:.1f} h (~{n_doubles*sec_per*5/86400:.1f} days).
- **Infeasible** without GPU acceleration (e.g. Uni-Dock, Vina-GPU 2.1) or
  cluster batch.  Even GPU Vina at ~0.3 s/dock would be ~{n_doubles*0.3*5/3600:.1f} h.

## Why we still ship the IDs
Even without docking, the enumerated lists are useful for:
- prioritising biologically interesting subsets (charge-reversal, gain-of-aromatic,
  loss-of-thiol-nucleophile);
- cross-referencing against ClinVar / COSMIC / gnomAD population variants;
- guiding a smaller targeted sub-sweep (e.g. only the {n_pos*4} singles whose
  functional_class is "charge_reversal" or "loss_of_charge").

## What we actually docked
Phase 5/7 dock 20 hand-picked mutants 1x and 8 priority mutants 5x.  The
multi-replica per-(target, box, ligand) numerical SD measured here was
**0.01-0.05 kcal/mol** at exhaustiveness 32 and box 18 A
(see `12_phase7/01_replicas/multi_replica_aggregate.csv`).  Note this is
*within-seed search reproducibility for one tuple*, not the published
Vina absolute-affinity noise floor (~0.85 kcal/mol; Trott & Olson 2010),
which is the right number to quote when comparing Vina deltaG against a
measured Kd.
"""
    (OUT / "feasibility_note.md").write_text(note)
    print(f"feasibility note written")

    # Plotly chemistry map (wider canvas + clearly separated legend / colorbar)
    df = pd.DataFrame(singles)
    fig = px.scatter(
        df,
        x="residue_position",
        y="hydropathy_change",
        color="volume_change",
        color_continuous_scale="RdBu_r",
        symbol="functional_class",
        hover_data=["mutation_id","wt_aa","new_aa","functional_class","volume_change","hydropathy_change"],
        labels={"residue_position":"Residue position",
                "hydropathy_change":"Δ Kyte-Doolittle (new − wt)",
                "volume_change":"Δ side-chain volume (Å³)"},
        height=720,
        width=1400,
    )
    fig.update_traces(marker=dict(size=11, line=dict(width=0.5, color="DarkSlateGrey")))
    fig.update_layout(
        template="plotly_white",
        # Title sits above the plot with breathing room — no axis collision
        title=dict(
            text=(f"<b>Active-site singles ({len(df)} mutations) at {n_pos} positions</b>"
                  "<br><sup>shape = functional class · colour = ΔV<sub>side-chain</sub> (Å³) · y = ΔKyte-Doolittle</sup>"),
            x=0.02, xanchor="left",
            y=0.97, yanchor="top",
            font=dict(size=15),
        ),
        # Generous right margin so legend AND colorbar fit side-by-side
        margin=dict(l=70, r=260, t=110, b=70),
        # Symbol legend (functional_class) -> top-right column
        legend=dict(
            title=dict(text="<b>Functional class</b>"),
            orientation="v",
            x=1.02, xanchor="left",
            y=1.0, yanchor="top",
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor="lightgrey", borderwidth=1,
        ),
        # Continuous colour-bar (Δvolume) -> outside the legend column, lower
        coloraxis_colorbar=dict(
            title=dict(text="ΔV<br>(Å³)", side="top"),
            x=1.18, xanchor="left",
            y=0.40, yanchor="middle",
            len=0.55, thickness=14,
        ),
        xaxis=dict(dtick=10),
    )
    html_path = OUT / "all_singles_chemistry_map.html"
    fig.write_html(html_path, include_plotlyjs="cdn")
    # Static PNG fallback (kaleido); skip silently if not available
    try:
        fig.write_image(str(OUT / "all_singles_chemistry_map.png"), width=1400, height=720, scale=2)
    except Exception as e:
        print(f"[warn] PNG fallback unavailable: {e}")
    print(f"plot written -> {html_path}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
