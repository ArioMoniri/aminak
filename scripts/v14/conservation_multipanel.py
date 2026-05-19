#!/usr/bin/env python3
"""Build a single multi-panel conservation figure matching the user's example layout:

    [   Top panel — JSD conservation × 313 residues, top-10% in red,    ]
    [   catalytic markers, p90 dashed line                              ]
    [                                                                    ]
    [  Phylogeny tree (NJ)        |   Per-residue mutation table        ]

Outputs:
  14_inhibitor_design/presentation/figures/conservation_multipanel.png
"""
from pathlib import Path
import csv, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle
import numpy as np

REPO = Path(__file__).resolve().parents[2]
CONS = REPO / "01b_msa_v2" / "conservation_scores.csv"
MSA = REPO / "01b_msa_v2" / "aligned.fa"
TREE_NWK = REPO / "12_phase7" / "05_phylogeny" / "tymS_tree.nwk"
OUT = REPO / "14_inhibitor_design" / "presentation" / "figures" / "conservation_multipanel.png"
OUT.parent.mkdir(parents=True, exist_ok=True)

# Load conservation
positions, js = [], []
with CONS.open() as f:
    for row in csv.DictReader(f):
        try:
            positions.append(int(row["ref_position"]))
            js.append(float(row["js_score"]))
        except ValueError: continue
js = np.array(js); positions = np.array(positions)
p90 = np.percentile(js, 90)

# Catalytic + clamp residues (mark with triangles)
catalytic = [(175, "R175"), (176, "R176"), (195, "C195"), (196, "H196"),
             (215, "R215"), (226, "N226"), (258, "Y258")]

# Phase 7 mutated-residue panel (a few representative ones to show in the table)
table_rows = [
    ("R50",  "R", "→Ala, →Glu",  "R = 10/10 (100 %, invariant)",
     "none", "Phosphate clamp"),
    ("F80",  "F", "→Ala, →Asp",  "F = 7/10 (70 %), H = 1/10, P = 1/10, A = 1/10",
     "E. coli = H · L. casei = P · phage T4 = A (all three bacterial / phage lineages)",
     "Pocket scaffold"),
    ("W109", "W", "→Ala",        "W = 10/10 (100 %, invariant)",
     "none", "Pocket scaffold"),
    ("T170", "T", "→Ala",        "T = 5/10 (50 %), N = 4/10, K = 1/10 — variable, exactly as expected for the distant-surface control",
     "E. coli, D. melanogaster, A. thaliana, P. falciparum all = N · phage T4 = K (the 5 mammalian / yeast / human-aligned lineages keep T)",
     "Distant control"),
    ("C195", "C", "→Ala",        "C = 10/10 (100 %, invariant) — catalytic Michael nucleophile",
     "none", "Catalytic dyad"),
    ("H196", "H", "→Ala",        "H = 10/10 (100 %, invariant)",
     "none", "Catalytic dyad"),
    ("R215", "R", "→Ala, →Glu",  "R = 10/10 (100 %, invariant)",
     "none", "Phosphate clamp"),
]

# Try to load the phylogeny
try:
    from Bio import Phylo
    tree = Phylo.read(str(TREE_NWK), "newick")
    have_tree = True
except Exception:
    have_tree = False

# ---- figure ----
fig = plt.figure(figsize=(17, 9.5))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1.4], width_ratios=[1, 2.4], hspace=0.45, wspace=0.18)

# Top: spans both columns — JSD plot
ax_top = fig.add_subplot(gs[0, :])
top10_mask = js >= np.percentile(js, 90)
colors = ["#c0392b" if t else "#3a86c8" for t in top10_mask]
ax_top.bar(positions, js, color=colors, width=1.0, linewidth=0)
ax_top.axhline(y=p90, ls="--", color="#888", linewidth=1.0, label=f"p90 = {p90:.3f}")
# catalytic markers (triangles at top)
ymax = js.max()
for pos, lbl in catalytic:
    if pos in positions:
        i = list(positions).index(pos)
        ax_top.plot(pos, js[i] + 0.012, marker="v", color="black", markersize=6)
        ax_top.annotate(lbl, xy=(pos, js[i] + 0.022), ha="center", fontsize=8, color="black")
ax_top.set_xlim(0, max(positions) + 2)
ax_top.set_ylim(0, ymax * 1.18)
ax_top.set_xlabel("Reference position (P04818)", fontsize=11)
ax_top.set_ylabel("JSD conservation (Capra-Singh, weighted window)", fontsize=10)
ax_top.set_title("v2 TYMS conservation across 10 orthologs (top 10% red, catalytic ▼)", fontsize=12)
ax_top.legend(loc="upper right", fontsize=9)
ax_top.grid(False)
for s in ("top", "right"): ax_top.spines[s].set_visible(False)

# Bottom-left: phylogeny
ax_tree = fig.add_subplot(gs[1, 0])
if have_tree:
    # ladderise + use coloured labels
    colour_map = {
        "Homo_sapiens": "#e74c3c", "Mus_musculus": "#e67e22", "Rattus_norvegicus": "#f1c40f",
        "Drosophila_melanogaster": "#9b59b6", "Saccharomyces_cerevisiae": "#8e44ad",
        "Escherichia_coli": "#2c3e50", "Lactobacillus_casei": "#7f8c8d",
        "Bacteriophage_T4": "#34495e", "Arabidopsis_thaliana": "#27ae60",
        "Plasmodium_falciparum": "#16a085",
    }
    # custom label_func to strip pipe
    def label_func(c):
        return (c.name or "").split("|")[0]
    # draw without confidence; can't set per-leaf colour easily with Phylo.draw, but ladderise first
    tree.ladderize()
    Phylo.draw(tree, axes=ax_tree, label_func=label_func, do_show=False, show_confidence=False)
    ax_tree.set_title("TYMS NJ tree (BLOSUM62 distance, 10 orthologs)", fontsize=11)
    # recolour terminal labels (matplotlib texts) to match the species
    for t in ax_tree.texts:
        s = t.get_text().strip()
        for k, c in colour_map.items():
            if s == k:
                t.set_color(c); t.set_fontsize(9); break
    for s in ("top","right"): ax_tree.spines[s].set_visible(False)
else:
    ax_tree.text(0.5, 0.5, "Phylogeny not available", ha="center", va="center")
    ax_tree.axis("off")

# Bottom-right: residue mutation table (custom matplotlib table)
ax_tbl = fig.add_subplot(gs[1, 1]); ax_tbl.axis("off")
header = ["Position", "WT", "Mutated to",
          "Conservation across\nthe 10 orthologs",
          "Which species deviate from WT?", "Class"]
cell_text = []
for pos, wt, mut, cons, dev, cls in table_rows:
    cell_text.append([pos, wt, mut, cons, dev, cls])

tbl = ax_tbl.table(cellText=cell_text, colLabels=header,
                    cellLoc="left", loc="upper center",
                    colWidths=[0.07, 0.05, 0.13, 0.27, 0.34, 0.14])
tbl.auto_set_font_size(False); tbl.set_fontsize(8.5)
tbl.scale(1.0, 2.6)
# style: header dark, alternating zebra
for (r, c), cell in tbl.get_celld().items():
    cell.set_linewidth(0.4)
    if r == 0:
        cell.set_facecolor("#21293A")
        cell.set_text_props(color="white", weight="bold", fontsize=9)
    else:
        cell.set_facecolor("#f7f4ee" if r % 2 == 1 else "#fdfbf7")
        if c == 0:
            cell.set_text_props(weight="bold", color="#21293A")

plt.savefig(OUT, dpi=140, bbox_inches="tight", facecolor="white")
print(f"→ {OUT}  ({OUT.stat().st_size/1024:.0f} KB)")
