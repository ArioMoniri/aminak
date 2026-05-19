#!/usr/bin/env python3
"""Re-render the Modeller-vs-AlphaFold comparison correctly — the v1 script
silently couldn't parse the column names and fell through to placeholder."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "14_inhibitor_design" / "07_advanced_methods" / "modeller_vs_alphafold.csv"
OUT = REPO / "14_inhibitor_design" / "presentation" / "figures" / "modeller_vs_alphafold.png"

rows = list(csv.DictReader(SRC.open()))
sources = [r["source"] for r in rows]
rmsd = [float(r["rmsd_vs_1HVY_align_A"]) for r in rows]
pct_fav = [float(r["pct_favoured"]) for r in rows]

colours = ["#e6553f" if "AlphaFold" in s else "#3a86c8" for s in sources]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Left: RMSD vs crystal
ax1.barh(range(len(sources)), rmsd, color=colours)
ax1.set_yticks(range(len(sources))); ax1.set_yticklabels(sources, fontsize=10)
ax1.set_xlabel("Cα RMSD vs 1HVY chain A (Å)  —  lower = closer to crystal")
ax1.set_title("Backbone fidelity vs the experimental crystal", fontsize=11)
ax1.grid(True, axis="x", alpha=0.3)
ax1.set_xlim(0, max(rmsd) * 1.30)
for i, v in enumerate(rmsd):
    ax1.annotate(f"{v:.3f} Å", xy=(v, i), xytext=(6, 0), textcoords="offset points",
                 va="center", fontsize=10, color="black")
for s in ("top","right"): ax1.spines[s].set_visible(False)

# Right: Lovell %favoured
ax2.barh(range(len(sources)), pct_fav, color=colours)
ax2.set_yticks(range(len(sources))); ax2.set_yticklabels(sources, fontsize=10)
ax2.set_xlabel("Ramachandran % favoured  (Lovell 4-map)")
ax2.set_title("Local geometry quality", fontsize=11)
ax2.axvline(x=92.2, ls="--", color="#888", alpha=0.7, label="1HVY crystal (92.2%)")
ax2.set_xlim(80, 100)
ax2.grid(True, axis="x", alpha=0.3); ax2.legend(loc="lower right", fontsize=9)
for i, v in enumerate(pct_fav):
    ax2.annotate(f"{v:.2f} %", xy=(v, i), xytext=(6, 0), textcoords="offset points",
                 va="center", fontsize=10, color="black")
for s in ("top","right"): ax2.spines[s].set_visible(False)

plt.suptitle("Modeller vs AlphaFold — all three beat the 1HVY crystal on Lovell favoured", fontsize=13)
plt.tight_layout(rect=(0, 0, 1, 0.95))
plt.savefig(OUT, dpi=140, facecolor="white")
plt.close()
print(f"→ {OUT}")
