#!/usr/bin/env python3
"""Build mutation-effect plots for v3 — FIX 10: legend keys match actual categories."""
from __future__ import annotations
import os, sys, pathlib, csv
import matplotlib.pyplot as plt

ROOT = pathlib.Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
CSV  = ROOT / "07c_mut_docking_v3" / "mutant_results_v3.csv"
OUT  = ROOT / "08c_analysis_v3"
OUT.mkdir(exist_ok=True)

if not CSV.exists():
    sys.exit(f"missing {CSV} — run Stage 7 v3 first")

rows = []
with open(CSV) as f:
    for line in f:
        if line.startswith("#"):
            continue
        # rewind: use DictReader from this content
        rest = line + f.read()
        import io
        for r in csv.DictReader(io.StringIO(rest)):
            rows.append(r)
        break

def num(r, k, default=0.0):
    try:
        return float(r.get(k) or default)
    except (TypeError, ValueError):
        return default

# FIX 10: TYPE_STYLE keys match actual categories from stage 7 v3
TYPE_STYLE = {
    "wildtype":                    ("#000000", "*", 200),
    "ala_scan":                    ("#1f77b4", "o", 80),
    "opposite":                    ("#ff7f0e", "s", 80),
    "arg_clamp":                   ("#2ca02c", "D", 80),
    "double_dyad":                 ("#9467bd", "^", 100),
    "double_phosclamp":            ("#d62728", "P", 100),
    "double_polar_neutral":        ("#8c564b", "v", 100),
    "double_substrate_orient":     ("#e377c2", "X", 100),
    "double_aromatic_swap":        ("#17becf", "<", 100),
    "control_surface":             ("#7f7f7f", ">", 100),
}

def style_for(cat):
    return TYPE_STYLE.get(cat, ("#999999", ".", 60))


def panel(ax, condition: str):
    title = f"Δ Vina score vs Pose-RMSD — {condition.upper()}"
    for r in rows:
        cond = r.get("condition") or "apo"
        if cond != condition:
            continue
        cat = r.get("category", "")
        col, mk, sz = style_for(cat)
        x = num(r, "delta_vina_vs_wt")
        y = num(r, "rmsd_to_native")
        # mis-docked: faded marker
        mis = (r.get("mis_docked", "False").lower() == "true")
        alpha = 0.35 if mis else 0.9
        ax.scatter(x, y, c=col, marker=mk, s=sz,
                   edgecolors="white", linewidths=0.7, alpha=alpha)
        mid = r.get("mutant", "")
        ax.annotate(mid, (x, y), fontsize=7.0, color="#222",
                    xytext=(4, 4), textcoords="offset points")

    ax.axhline(3.0, color="#888", lw=0.8, ls="--")
    ax.axhline(2.0, color="#ccc", lw=0.6, ls=":")
    ax.axvline(0.0, color="#888", lw=0.8, ls="--")
    ax.set_xlabel("Δ Vina score vs WT (kcal/mol; +ve = destabilising)")
    ax.set_ylabel("Top-pose RMSD vs crystal dUMP (Å)")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.text(0.99, 0.98, "destabilising +\npose-displacing",
            transform=ax.transAxes, ha="right", va="top",
            color="#a33", fontsize=8.5, alpha=0.8)
    ax.text(0.01, 0.02, "stabilising +\npose-preserving",
            transform=ax.transAxes, ha="left", va="bottom",
            color="#393", fontsize=8.5, alpha=0.8)


fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
panel(axes[0], "apo")
panel(axes[1], "holo")

# Build legend from observed categories (in TYPE_STYLE)
seen = []
for r in rows:
    cat = r.get("category", "")
    if cat and cat not in [t for t, _ in seen]:
        seen.append((cat, style_for(cat)))
handles = [plt.Line2D([], [], marker=mk, color=col, linestyle="",
                      markersize=8, markeredgecolor="white",
                      label=cat.replace("_", " "))
           for cat, (col, mk, sz) in seen]
fig.legend(handles=handles, loc="upper center",
           bbox_to_anchor=(0.5, -0.01), ncol=min(5, len(handles)), frameon=False, fontsize=9)
fig.suptitle("v3 Mutation-effect map — Δ Vina vs pose RMSD (faded = mis-docked, RMSD > 3 Å)",
             fontsize=13, fontweight="bold", color="#1f3b5e")
fig.tight_layout(rect=[0, 0.03, 1, 0.96])
out = OUT / "mutation_effect_plot.png"
fig.savefig(out, dpi=160, bbox_inches="tight", facecolor="white")
print("wrote", out, "- categories:", [c for c, _ in seen])

# Concordance scatter
apo_by = {r["mutant"]: num(r, "delta_vina_vs_wt") for r in rows if (r.get("condition") or "apo") == "apo"}
hol_by = {r["mutant"]: num(r, "delta_vina_vs_wt") for r in rows if (r.get("condition") or "apo") == "holo"}
common = sorted(set(apo_by) & set(hol_by))
fig2, ax = plt.subplots(figsize=(7.5, 6))
xs, ys = [], []
for mid in common:
    cat = next((r.get("category", "") for r in rows if r["mutant"] == mid), "")
    col, mk, sz = style_for(cat)
    x = apo_by[mid]; y = hol_by[mid]
    xs.append(x); ys.append(y)
    ax.scatter(x, y, c=col, marker=mk, s=sz, edgecolors="white", linewidths=0.7, alpha=0.9)
    ax.annotate(mid, (x, y), fontsize=7, color="#222",
                xytext=(4, 4), textcoords="offset points")
if xs:
    lo = min(xs + ys) - 0.2
    hi = max(xs + ys) + 0.2
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, alpha=0.6, label="y = x")
ax.axhline(0, color="#bbb", lw=0.5); ax.axvline(0, color="#bbb", lw=0.5)
ax.set_xlabel("Δ Vina (apo) — kcal/mol")
ax.set_ylabel("Δ Vina (holo) — kcal/mol")
ax.set_title("v3 Apo vs Holo concordance", fontsize=11, fontweight="bold", color="#1f3b5e")
ax.legend(loc="lower right", fontsize=9)
fig2.tight_layout()
out2 = OUT / "apo_vs_holo_concordance.png"
fig2.savefig(out2, dpi=160, bbox_inches="tight", facecolor="white")
print("wrote", out2)
