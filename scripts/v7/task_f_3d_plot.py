#!/usr/bin/env python3
"""Task F: master 3D dynamic plot.

x = residue position (mutated site)
y = chemistry index (sum hydropathy change for the mutation)
z = delta Vina (kcal/mol)
color = functional class
size = abs(delta Vina)
"""
from __future__ import annotations
import json
import re
import sys
import time
from pathlib import Path

import pandas as pd
import plotly.express as px

import os
PROJECT = Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
OUT = PROJECT / "12_phase7" / "06_3d_plot"
LOG = PROJECT / "logs" / "v7_task_f.log"
PIPELOG = PROJECT / "pipeline.log"

KD = {"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,
      "K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,
      "T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}

def parse_mut(name: str) -> list[tuple[str,int,str]]:
    out = []
    # patterns: A123B  optionally joined by _
    for tok in name.split("_"):
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", tok)
        if m:
            out.append((m.group(1), int(m.group(2)), m.group(3)))
    return out


def log(msg: str) -> None:
    line = f"[V7][taskF] {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
    print(line, flush=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a") as fh: fh.write(line+"\n")
    with PIPELOG.open("a") as fh: fh.write(line+"\n")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    log("Task F start")
    df = pd.read_csv(PROJECT / "07e_mut_docking_v5" / "mutant_results_v5.csv", comment="#")
    df = df[(df.condition == "holo") & (df.mutant != "WT")].copy()
    rows = []
    for _, r in df.iterrows():
        mut = parse_mut(r["mutant"])
        if not mut:
            continue
        # Use the position of the FIRST mutated residue as x
        for wt, pos, new in mut:
            sum_dkd = sum(KD[n] - KD[w] for w,_,n in mut)
            rows.append({
                "mutant_id": r["mutant"],
                "category": r["category"],
                "residue_position": pos,
                "wt_aa": wt,
                "new_aa": new,
                "sum_hydropathy_change": round(sum_dkd, 3),
                "delta_vina": float(r["delta_vina_vs_wt"]),
                "abs_delta_vina": abs(float(r["delta_vina_vs_wt"])),
                "top_affinity": float(r["top_affinity"]),
                "n_mutated_residues": len(mut),
            })
    plot_df = pd.DataFrame(rows)
    log(f"plotting {len(plot_df)} points")
    fig = px.scatter_3d(
        plot_df, x="residue_position", y="sum_hydropathy_change", z="delta_vina",
        color="category", size="abs_delta_vina",
        hover_data=["mutant_id","wt_aa","new_aa","top_affinity","n_mutated_residues"],
        title="Phase 5 mutants in (position, chemistry, dVina) space",
        labels={"residue_position":"Residue position",
                "sum_hydropathy_change":"Sum delta Kyte-Doolittle",
                "delta_vina":"Delta Vina vs WT (kcal/mol)"},
        height=750,
    )
    fig.update_traces(marker=dict(line=dict(color="black", width=0.5), opacity=0.85))
    fig.update_layout(scene=dict(zaxis=dict(title="Delta Vina (kcal/mol)")),
                      template="plotly_white")
    html = OUT / "mutation_3d.html"
    fig.write_html(html, include_plotlyjs="cdn")
    log(f"plotly 3D -> {html}")
    # PNG snapshot via kaleido if available
    png = OUT / "mutation_3d.png"
    try:
        import kaleido  # noqa: F401
        fig.write_image(png, scale=2, width=1400, height=900)
        log(f"png -> {png}")
    except Exception as e:
        (OUT / "PNG_FALLBACK.txt").write_text(
            f"PNG snapshot skipped: {e}\nUse the HTML viewer or run kaleido manually.")
        log(f"png skipped: {e}")
    plot_df.to_csv(OUT / "mutation_3d_data.csv", index=False)
    log("Task F done")
    return 0

if __name__ == "__main__":
    sys.exit(main())
