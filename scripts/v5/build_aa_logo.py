#!/usr/bin/env python3
"""
Sequence-logo plot of the active-site residues in TYMS, computed from the
v2 verified-orthologs MSA (01b_msa_v2/aligned.fa). Letter size = frequency
× information content (bits). Uses the same conventions as WebLogo /
MolProbity / logomaker:
  - Stack height at column j = log2(20) − Shannon_entropy(j) bits
                              = sum over aa of p(aa,j) * log2(20*p(aa,j))
  - Letter height within stack = p(aa,j) × stack height
  - Most conserved → tallest single letter; variable → multiple short letters.

Highlights residues we mutated in the docking panel with a coloured band
beneath each column, and annotates the WT residue + the substitution(s)
that we tried.

Output:
  11_enhanced/aa_logo_active_site.png
"""
from __future__ import annotations
import os, pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import logomaker as lm
from Bio import SeqIO

ROOT = pathlib.Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
MSA  = ROOT / "01b_msa_v2" / "aligned.fa"
OUT  = ROOT / "11_enhanced" / "aa_logo_active_site.png"
OUT2 = ROOT / "11_enhanced" / "aa_logo_full_chain.png"

# Active-site / mutation-panel residues (numbering = ungapped P04818)
PANEL = {
    50:  ("R", ["A", "E"], "Phosphate clamp"),
    80:  ("F", ["A", "D"], "Pocket scaffold"),
    109: ("W", ["A"],      "Pocket scaffold"),
    170: ("T", ["A"],      "Distant control"),
    175: ("R", ["A", "E"], "Phosphate clamp"),
    176: ("R", ["A", "E"], "Phosphate clamp"),
    195: ("C", ["A", "S"], "Catalytic nucleophile"),
    196: ("H", ["A", "F"], "Catalytic dyad"),
    214: ("Q", ["A"],      "Pocket scaffold"),
    215: ("R", ["A", "E"], "Phosphate clamp"),
    218: ("D", ["A", "K"], "Pocket scaffold"),
    225: ("F", ["A", "D"], "Pocket scaffold"),
    226: ("N", ["A", "D"], "Substrate orientation"),
    258: ("Y", ["A", "F"], "Substrate orientation"),
}
COLOR_OF = {
    "Catalytic nucleophile":  "#d62728",
    "Catalytic dyad":         "#d62728",
    "Phosphate clamp":        "#1f77b4",
    "Substrate orientation":  "#2ca02c",
    "Pocket scaffold":        "#ff7f0e",
    "Distant control":        "#7f7f7f",
}

# 1) Parse MSA → list of aligned sequences
seqs = list(SeqIO.parse(str(MSA), "fasta"))
print(f"  loaded {len(seqs)} aligned orthologs")
ref = next(s for s in seqs if "P04818" in s.id)
ref_aln = str(ref.seq)

# 2) Build column → ungapped reference position map
col_to_ref = {}      # alignment column index → ungapped 1-based residue number
ungapped_seq = []    # sequence with ref residue at each kept position
ref_pos = 0
for col_idx, ch in enumerate(ref_aln):
    if ch == "-": continue
    ref_pos += 1
    col_to_ref[col_idx] = ref_pos
    ungapped_seq.append(ch)
print(f"  reference (P04818): {len(ungapped_seq)} ungapped residues")

# 3) Build a per-column count matrix in alignment columns where the reference is non-gap
aas = "ACDEFGHIKLMNPQRSTVWY"
rows = []
for col_idx in sorted(col_to_ref):
    counts = {a: 0 for a in aas}
    for s in seqs:
        c = str(s.seq)[col_idx].upper()
        if c in counts: counts[c] += 1
    rows.append(counts)
df = pd.DataFrame(rows, index=[col_to_ref[c] for c in sorted(col_to_ref)])
df.index.name = "ref_position"

# 4) Convert counts → probability → information content
prob = df.div(df.sum(axis=1).replace(0, 1), axis=0)
# Shannon entropy per column
entropy = -(prob * np.where(prob > 0, np.log2(prob), 0)).sum(axis=1)
ic = np.log2(20) - entropy
# Letter heights
heights = prob.mul(ic, axis=0)

# ===================== PLOT 1: focus on active-site residues =====================
# Re-index by COMPACT positions (1, 2, 3, ...) so columns don't crowd
positions = sorted(PANEL.keys())
compact_idx = list(range(len(positions)))
sub = heights.loc[positions].reset_index(drop=True)

# Two-stack figure: top = logo, bottom = function-class band, with extra
# vertical room for the "WT=X / →A,B" labels above each column.
fig, (ax_top, ax) = plt.subplots(2, 1, figsize=(16, 7),
                                  gridspec_kw={"height_ratios": [1.4, 6]},
                                  sharex=True)
ax_top.axis("off")
logo = lm.Logo(sub, ax=ax, color_scheme="chemistry", font_name="DejaVu Sans Bold")

# Use COMPACT x positions; relabel ticks with the actual residue numbers
ax.set_xticks(compact_idx)
ax.set_xticklabels([f"{wt}{p}" for p, (wt,_,_) in zip(positions, [PANEL[p] for p in positions])],
                   fontsize=10, fontweight="bold")
ax.set_xlabel("Active-site / mutated residues  (numbering = ungapped P04818)",
              fontsize=11, labelpad=10)
ax.set_ylabel("Information content (bits)\nlog₂(20) − Shannon entropy",
              fontsize=11)
ax.set_ylim(0, np.log2(20) + 0.2)
ax.axhline(np.log2(20), ls="--", color="#bbb", lw=0.6)
ax.text(compact_idx[-1] + 0.4, np.log2(20) - 0.05, "fully conserved (4.32 bits)",
        ha="left", va="top", fontsize=8, color="#666")

# Functional-class colour band IN A SEPARATE STRIP just below the logo
band_y = -0.50
band_h = 0.40
for x, p in zip(compact_idx, positions):
    cls = PANEL[p][2]
    ax.add_patch(mpatches.Rectangle((x - 0.45, band_y), 0.90, band_h,
                                    facecolor=COLOR_OF[cls], alpha=0.85,
                                    edgecolor="white", linewidth=0.5,
                                    transform=ax.transData, clip_on=False))

# Mutation panel labels in the TOP axis (no overlap with logo)
ax_top.set_xlim(ax.get_xlim())
for x, p in zip(compact_idx, positions):
    wt, muts, _ = PANEL[p]
    ax_top.text(x, 0.65, f"WT={wt}", ha="center", va="bottom",
                fontsize=9, fontweight="bold", color="#1f3b5e")
    ax_top.text(x, 0.20, f"→{', '.join(muts)}", ha="center", va="bottom",
                fontsize=8, color="#b8593c")
ax_top.set_ylim(0, 1)

# Class legend (below x label)
handles = [mpatches.Patch(color=col, label=name) for name, col in COLOR_OF.items()]
ax.legend(handles=handles, loc="upper center", ncol=3,
          bbox_to_anchor=(0.5, -0.30), fontsize=9, frameon=False,
          title="Functional class (band colour)", title_fontsize=9)

fig.suptitle(
    "Amino-acid sequence logo at TYMS active-site & mutated residues\n"
    "Letter height = frequency × information content (10 verified TYMS orthologs)",
    fontsize=13, fontweight="bold", color="#1f3b5e", y=0.98,
)
plt.subplots_adjust(left=0.07, right=0.98, bottom=0.18, top=0.88, hspace=0.05)
plt.savefig(OUT, dpi=180, bbox_inches="tight", facecolor="white")
print(f"  wrote {OUT} ({OUT.stat().st_size//1024} KB)")

# ===================== PLOT 2: full chain logo (small, for reference) =====================
fig2, ax2 = plt.subplots(figsize=(20, 2.5))
lm.Logo(heights, ax=ax2, color_scheme="chemistry", font_name="DejaVu Sans Bold")
ax2.set_xlabel("Reference position (P04818)", fontsize=10)
ax2.set_ylabel("bits", fontsize=10)
# Mark active-site columns
for p in positions:
    ax2.axvspan(p - 0.5, p + 0.5, alpha=0.15, color=COLOR_OF[PANEL[p][2]], zorder=0)
ax2.set_xlim(0, len(heights) + 1)
ax2.set_ylim(0, np.log2(20) + 0.2)
ax2.set_title(
    "Full-chain TYMS sequence logo (P04818 numbering); active-site columns highlighted",
    fontsize=11, fontweight="bold", color="#1f3b5e",
)
plt.subplots_adjust(bottom=0.20, top=0.85)
plt.savefig(OUT2, dpi=160, bbox_inches="tight", facecolor="white")
print(f"  wrote {OUT2} ({OUT2.stat().st_size//1024} KB)")

# ===================== Per-position frequency table for the README =====================
table_csv = ROOT / "11_enhanced" / "aa_logo_active_site.csv"
rows = []
for p in positions:
    counts = df.loc[p].to_dict()
    total = sum(counts.values())
    top3 = sorted(counts.items(), key=lambda kv: -kv[1])[:3]
    pct = ", ".join(f"{aa}={cnt}/{total} ({100*cnt/total:.0f}%)" for aa, cnt in top3 if cnt > 0)
    wt, muts, cls = PANEL[p]
    is_invariant = (counts.get(wt, 0) == total)
    rows.append(dict(
        position=p, wt=wt, ic_bits=round(float(ic.loc[p]), 3),
        invariant=is_invariant,
        top3_observed=pct,
        functional_class=cls,
        mutated_to=",".join(muts),
    ))
pd.DataFrame(rows).to_csv(table_csv, index=False)
print(f"  wrote {table_csv}")
