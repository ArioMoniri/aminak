#!/usr/bin/env python3
"""
Build "Ramachandran-style" 2D mutation-effect plots for v2 results:
- x-axis: ΔVina score vs WT (apo)
- y-axis: pose RMSD vs crystal dUMP
- one panel for apo, one for holo
- residues labelled, type-coloured
- quadrants annotated

Plus a paired-condition apo/holo concordance scatter.

Reads:  07b_mut_docking_v2/results_full.csv
Writes: 08b_analysis_v2/mutation_effect_plot.png (panel)
        08b_analysis_v2/apo_vs_holo_concordance.png

Run from project root:
    python scripts/v2/build_mutation_plot.py
"""
from __future__ import annotations
import os, sys, pathlib, csv
import matplotlib.pyplot as plt
import matplotlib.patches as mp

ROOT = pathlib.Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
CSV  = ROOT / "07b_mut_docking_v2" / "mutant_results_v2.csv"
OUT  = ROOT / "08b_analysis_v2"
OUT.mkdir(exist_ok=True)

if not CSV.exists():
    sys.exit(f"missing {CSV} — run Stage 7 v2 first")

rows = []
with open(CSV) as f:
    for r in csv.DictReader(f):
        rows.append(r)

def num(r, k, default=0.0):
    try:
        return float(r.get(k) or default)
    except (TypeError, ValueError):
        return default

# Type → color/marker
TYPE_STYLE = {
    "wildtype":         ("#000000", "*", 200),
    "single_ala":       ("#1f77b4", "o", 80),
    "single_opposite":  ("#ff7f0e", "s", 80),
    "arg_clamp":        ("#2ca02c", "D", 80),
    "double_catalytic": ("#9467bd", "^", 100),
    "double":           ("#d62728", "P", 100),
    "control_surface":  ("#7f7f7f", "X", 100),
}
def style_for(typ):
    for k, v in TYPE_STYLE.items():
        if k in typ:
            return v
    return ("#999999", ".", 60)


def panel(ax, condition: str):
    title = f"Δ-Score vs Pose-RMSD — {condition.upper()}"
    plotted_types = set()
    for r in rows:
        cond = r.get("condition") or "apo"
        if cond != condition:
            continue
        typ = r.get("category", "")
        col, mk, sz = style_for(typ)
        x = num(r, "ddG_vs_wt")
        y = num(r, "rmsd_to_native")
        ax.scatter(x, y, c=col, marker=mk, s=sz,
                   edgecolors="white", linewidths=0.7, alpha=0.9)
        # Label
        mid = r.get("mutant", "")
        ax.annotate(mid, (x, y), fontsize=7.0, color="#222",
                    xytext=(4, 4), textcoords="offset points")
        plotted_types.add(typ)

    ax.axhline(2.0, color="#888", lw=0.8, ls="--")
    ax.axvline(0.0, color="#888", lw=0.8, ls="--")
    ax.set_xlabel("ΔVina score vs WT (kcal/mol)")
    ax.set_ylabel("Top-pose RMSD vs crystal dUMP (Å)")
    ax.set_title(title, fontsize=11, fontweight="bold")
    # Quadrant annotations
    ax.text(0.99, 0.98, "destabilising +\npose-displacing",
            transform=ax.transAxes, ha="right", va="top",
            color="#a33", fontsize=8.5, alpha=0.8)
    ax.text(0.01, 0.98, "stabilising score\n(but pose-displaced)",
            transform=ax.transAxes, ha="left", va="top",
            color="#888", fontsize=8.5, alpha=0.8)
    ax.text(0.01, 0.02, "stabilising +\npose-preserving",
            transform=ax.transAxes, ha="left", va="bottom",
            color="#393", fontsize=8.5, alpha=0.8)
    ax.text(0.99, 0.02, "destabilising score\n(but pose-preserved)",
            transform=ax.transAxes, ha="right", va="bottom",
            color="#888", fontsize=8.5, alpha=0.8)


fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
panel(axes[0], "apo")
panel(axes[1], "holo")

# Build a shared legend from observed types
seen = []
for r in rows:
    typ = r.get("category", "")
    if typ and typ not in [t for t, _ in seen]:
        seen.append((typ, style_for(typ)))
handles = [plt.Line2D([], [], marker=mk, color=col, linestyle="",
                      markersize=8, markeredgecolor="white",
                      label=typ.replace("_", " "))
           for typ, (col, mk, sz) in seen]
fig.legend(handles=handles, loc="upper center",
           bbox_to_anchor=(0.5, -0.01), ncol=min(7, len(handles)), frameon=False, fontsize=9)
fig.suptitle("Mutation-effect map — Vina score change vs pose RMSD",
             fontsize=13, fontweight="bold", color="#1f3b5e")
fig.tight_layout(rect=[0, 0.03, 1, 0.97])
out_panel = OUT / "mutation_effect_plot.png"
fig.savefig(out_panel, dpi=160, bbox_inches="tight", facecolor="white")
print("wrote", out_panel)

# --- apo vs holo concordance scatter ---
apo_by_id = {r["mutant"]: num(r, "ddG_vs_wt") for r in rows if (r.get("condition") or "apo") == "apo"}
hol_by_id = {r["mutant"]: num(r, "ddG_vs_wt") for r in rows if (r.get("condition") or "apo") == "holo"}
common = sorted(set(apo_by_id) & set(hol_by_id))
fig2, ax = plt.subplots(figsize=(7.5, 6))
xs, ys = [], []
for mid in common:
    typ = next((r.get("category", "") for r in rows if r["mutant"] == mid), "")
    col, mk, sz = style_for(typ)
    x = apo_by_id[mid]; y = hol_by_id[mid]
    xs.append(x); ys.append(y)
    ax.scatter(x, y, c=col, marker=mk, s=sz,
               edgecolors="white", linewidths=0.7, alpha=0.9)
    ax.annotate(mid, (x, y), fontsize=7, color="#222",
                xytext=(4, 4), textcoords="offset points")
lo = min(xs + ys) - 0.2 if xs else -1
hi = max(xs + ys) + 0.2 if xs else 1
ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.6, label="y = x (perfect concordance)")
ax.axhline(0, color="#bbb", lw=0.5); ax.axvline(0, color="#bbb", lw=0.5)
ax.set_xlabel("Δ score (apo) — kcal/mol")
ax.set_ylabel("Δ score (holo) — kcal/mol")
ax.set_title("Apo vs holo — does the cofactor change the mutant ranking?",
             fontsize=11, fontweight="bold", color="#1f3b5e")
ax.legend(loc="lower right", fontsize=9)
fig2.tight_layout()
out_conc = OUT / "apo_vs_holo_concordance.png"
fig2.savefig(out_conc, dpi=160, bbox_inches="tight", facecolor="white")
print("wrote", out_conc)
