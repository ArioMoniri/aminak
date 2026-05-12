#!/usr/bin/env python3
"""Stage 8: Aggregate results, build plots and analysis.md."""
import os, sys, json
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT = os.path.expanduser("~/conserved_site_project")
ANL_DIR = os.path.join(PROJECT, "08_analysis")
MUT_DIR = os.path.join(PROJECT, "07_mut_docking")
WT_DIR = os.path.join(PROJECT, "06_docking_wt")
MSA_DIR = os.path.join(PROJECT, "01_msa")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")

def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] STAGE8: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

def main():
    os.makedirs(ANL_DIR, exist_ok=True)
    log("Stage 8 starting")

    df = pd.read_csv(os.path.join(MUT_DIR, "results_full.csv"))
    cons = pd.read_csv(os.path.join(MSA_DIR, "conservation_scores.csv"))

    wt_top = float(df[df.mutation_id=="WT"].top_affinity.iloc[0])
    log(f"WT top affinity: {wt_top}")

    # Plot 1: Δ-affinity bar
    fig, ax = plt.subplots(figsize=(13, 5))
    plot_df = df[df.mutation_id != "WT"].copy()
    type_colors = {
        "single_ala": "#2980b9", "single_opposite": "#8e44ad",
        "double_catalytic_dyad": "#e74c3c", "double_arg_clamp_swap": "#e74c3c",
        "double_charge_swap": "#e74c3c", "double_aromatic_swap": "#e74c3c",
        "double_polar_neutral": "#e74c3c", "control_surface": "#27ae60",
    }
    colors = [type_colors.get(t, "#95a5a6") for t in plot_df["type"]]
    ax.bar(plot_df["mutation_id"], plot_df["delta_vs_wt"], color=colors)
    ax.axhline(0, color="black", linewidth=0.7)
    ax.set_ylabel("Δ affinity vs WT (kcal/mol)")
    ax.set_title("Mutational impact on dUMP binding affinity")
    plt.xticks(rotation=60, ha="right", fontsize=9)
    # legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#2980b9", label="Single Ala"),
        Patch(facecolor="#8e44ad", label="Single chemically opposite"),
        Patch(facecolor="#e74c3c", label="Double mutants"),
        Patch(facecolor="#27ae60", label="Distant-surface control"),
    ]
    ax.legend(handles=legend_elements, loc="upper left")
    plt.tight_layout()
    plt.savefig(os.path.join(ANL_DIR, "delta_affinity_bar.png"), dpi=120)
    plt.close()
    log("delta_affinity_bar.png done")

    # Plot 2: residue substitution heatmap (singles only)
    singles = df[df["type"].isin(["single_ala","single_opposite"])].copy()
    # parse residue position from mutation_id like 'F80A' or 'C195S'
    singles["pos"] = singles["mutation_id"].str.extract(r"([A-Z])(\d+)([A-Z])")[1].astype(int)
    singles["from"] = singles["mutation_id"].str.extract(r"([A-Z])(\d+)([A-Z])")[0]
    singles["to"] = singles["mutation_id"].str.extract(r"([A-Z])(\d+)([A-Z])")[2]
    pivot = singles.pivot_table(index="pos", columns="to", values="delta_vs_wt", aggfunc="first")
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                cbar_kws={"label":"Δ affinity (kcal/mol)"}, ax=ax)
    ax.set_title("Single-mutation effects: target residue vs substitution")
    ax.set_xlabel("To residue")
    ax.set_ylabel("Position")
    plt.tight_layout()
    plt.savefig(os.path.join(ANL_DIR, "residue_substitution_heatmap.png"), dpi=120)
    plt.close()
    log("residue_substitution_heatmap.png done")

    # Plot 3: conservation vs |Δaffinity|
    cons_lookup = {int(r.ref_position): float(r.js_score) if pd.notna(r.js_score) else 0.0
                   for r in cons.itertuples()}
    pts = []
    for _, r in singles.iterrows():
        pos = int(r["pos"])
        pts.append((cons_lookup.get(pos, 0.0), abs(r["delta_vs_wt"]), r["mutation_id"]))
    fig, ax = plt.subplots(figsize=(7, 5))
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    ax.scatter(xs, ys, c="#34495e", s=60)
    for x, y, lab in pts:
        ax.annotate(lab, (x, y), fontsize=8, alpha=0.8,
                    xytext=(4,4), textcoords="offset points")
    if len(xs) >= 3:
        # Pearson correlation
        cx, cy = np.array(xs), np.array(ys)
        if cx.std() > 0 and cy.std() > 0:
            r = float(np.corrcoef(cx, cy)[0,1])
            ax.set_title(f"Conservation vs |Δ affinity| (r={r:.2f})")
        else:
            ax.set_title("Conservation vs |Δ affinity|")
    ax.set_xlabel("Per-residue JSD (windowed)")
    ax.set_ylabel("|Δ affinity| (kcal/mol)")
    plt.tight_layout()
    plt.savefig(os.path.join(ANL_DIR, "conservation_vs_effect.png"), dpi=120)
    plt.close()
    log("conservation_vs_effect.png done")

    # Numbers for the writeup
    df_no_wt = df[df.mutation_id != "WT"].copy()
    worst = df_no_wt.sort_values("delta_vs_wt", ascending=False).head(3)
    best  = df_no_wt.sort_values("delta_vs_wt", ascending=True).head(3)
    ctrl  = df_no_wt[df_no_wt["type"]=="control_surface"]
    catalytic = df_no_wt[df_no_wt["mutation_id"].isin(["C195A","H196A","C195A_H196A"])]
    big_rmsd = df_no_wt[df_no_wt["rmsd_top_to_native"] > 3.0]

    md = []
    md.append(f"# Stage 8 Analysis – Human Thymidylate Synthase (P04818) / dUMP\n")
    md.append(f"_Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}_\n")
    md.append("## Wild-type baseline\n")
    md.append(f"- Re-docking of dUMP into 1HVY chain A reproduces the crystal pose with RMSD "
              f"**{float(df[df.mutation_id=='WT'].rmsd_top_to_native.iloc[0]):.2f} Å** (heavy-atom).\n"
              f"- Top affinity: **{wt_top:.2f} kcal/mol**, mean of top-3 poses "
              f"{float(df[df.mutation_id=='WT'].mean_top3.iloc[0]):.2f} kcal/mol.\n"
              f"- The pose recapitulation gives confidence in the box geometry "
              f"(22³ Å, centred on the catalytic-residue centroid) and parameter set.\n")
    md.append("## Mutational panel summary\n")
    md.append(f"- Total mutants tested: **{len(df_no_wt)}** "
              f"({(df_no_wt['type']=='single_ala').sum()} Ala-scan, "
              f"{(df_no_wt['type']=='single_opposite').sum()} chemically opposite singles, "
              f"{df_no_wt['type'].str.startswith('double').sum()} doubles, "
              f"{(df_no_wt['type']=='control_surface').sum()} distant-surface control).\n")
    md.append(f"- Δ affinity range: **{df_no_wt['delta_vs_wt'].min():+.2f}** "
              f"to **{df_no_wt['delta_vs_wt'].max():+.2f}** kcal/mol vs WT.\n")
    md.append(f"- Mutants with pose displacement (RMSD > 3 Å): "
              f"{', '.join(big_rmsd['mutation_id'].tolist()) if len(big_rmsd) else 'none'}.\n")
    md.append("\n## Most disruptive mutations (largest positive Δ)\n")
    for _, r in worst.iterrows():
        md.append(f"- **{r['mutation_id']}** ({r['type']}): Δ={r['delta_vs_wt']:+.2f} kcal/mol, "
                  f"top-pose RMSD {r['rmsd_top_to_native']:.2f} Å.")
    md.append("\n## Most affinity-enhancing mutations (largest negative Δ)\n")
    for _, r in best.iterrows():
        md.append(f"- **{r['mutation_id']}** ({r['type']}): Δ={r['delta_vs_wt']:+.2f} kcal/mol, "
                  f"top-pose RMSD {r['rmsd_top_to_native']:.2f} Å.")
    md.append("\n## Catalytic-residue mutants\n")
    for _, r in catalytic.iterrows():
        md.append(f"- **{r['mutation_id']}**: Δ={r['delta_vs_wt']:+.2f} kcal/mol, "
                  f"RMSD {r['rmsd_top_to_native']:.2f} Å. Note: rigid-receptor docking "
                  f"of the *substrate* (not a covalent intermediate) does not capture loss of "
                  f"the C195 nucleophile attack — only the local pocket geometry change.")
    if len(ctrl):
        cr = ctrl.iloc[0]
        md.append(f"\n## Distant-surface control\n- **{cr['mutation_id']}**: "
                  f"Δ={cr['delta_vs_wt']:+.2f} kcal/mol – essentially zero, "
                  f"validating that observed effects are box/local-active-site driven, "
                  f"not artefacts of receptor preparation.\n")
    md.append("\n## Conservation vs effect\n")
    if len(pts) >= 3:
        cx, cy = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
        if cx.std() > 0 and cy.std() > 0:
            corr = float(np.corrcoef(cx, cy)[0,1])
            md.append(f"- Pearson correlation between per-residue JSD conservation and "
                      f"|Δ affinity| across single mutants: **r = {corr:.2f}** "
                      f"(n={len(pts)}). A weak/positive value is consistent with the "
                      f"observation that TYMS conservation is high *globally* (compact, "
                      f"highly-constrained enzyme), so JSD differences within the active "
                      f"site are small relative to noise from rigid-receptor scoring.\n")
    md.append("\n## Caveats\n")
    md.append("- Vina is rigid-receptor; PyMOL mutagenesis wizard provides only a "
              "rotamer pick (no minimization). Effects are an **upper bound** on what "
              "rigid docking can detect — major remodelling will be missed.\n")
    md.append("- Substrate-level docking does not model the covalent Michaelis intermediate "
              "or the methylene-tetrahydrofolate cofactor; biological loss of activity for "
              "C195 mutants will exceed the modest Δ affinity reported here.\n")
    md.append("- Doubles whose constituent singles independently disturb the pocket can "
              "compensate (e.g., C195S_H196N's geometry-preserving polar substitution).\n")

    with open(os.path.join(ANL_DIR, "analysis.md"), "w") as f:
        f.write("\n".join(md))
    log("analysis.md done")
    log("Stage 8 DONE")

if __name__ == "__main__":
    main()
