#!/usr/bin/env python3
"""Stage 8 v3: Analysis with corrected sign convention, mean_topk, mis-dock annotation.

Fixes:
- FIX 3: positive delta_vina = destabilising (already in CSV)
- FIX 5: mean_topk column already in CSV
- FIX 8: mis_docked column annotated; excluded from top destabilisers ranking
- FIX 9: column name and prose use 'delta_vina' / 'Delta Vina score'
"""
import os, json, math
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = os.path.expanduser("~/conserved_site_project")
ANA_DIR = os.path.join(PROJECT, "08c_analysis_v3")
MUT_DIR = os.path.join(PROJECT, "07c_mut_docking_v3")
WT_DIR = os.path.join(PROJECT, "06c_docking_wt_v3")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v3_08_analysis.log")


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V3] STAGE8: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def main():
    os.makedirs(ANA_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 8 v3 starting")

    df = pd.read_csv(os.path.join(MUT_DIR, "mutant_results_v3.csv"), comment="#")
    log(f"loaded {len(df)} rows")

    # Pivot
    pv = df.pivot_table(index=["mutant", "category"], columns="condition",
                        values=["top_affinity", "delta_vina_vs_wt", "rmsd_to_native",
                                "mis_docked", "n_modes"],
                        aggfunc="first").reset_index()
    pv.columns = ["_".join([str(x) for x in c if x is not None and str(x) != ""]).strip("_")
                  for c in pv.columns]
    pv.to_csv(os.path.join(ANA_DIR, "mutant_pivot_v3.csv"), index=False)
    log(f"wrote pivot: {pv.columns.tolist()[:8]}...")

    # Ranked tables (sorted by delta_vina_vs_wt descending; positive=destabilising)
    df_apo = df[df.condition == "apo"].sort_values("delta_vina_vs_wt", ascending=False).reset_index(drop=True)
    df_holo = df[df.condition == "holo"].sort_values("delta_vina_vs_wt", ascending=False).reset_index(drop=True)
    df_apo.to_csv(os.path.join(ANA_DIR, "ranked_apo.csv"), index=False)
    df_holo.to_csv(os.path.join(ANA_DIR, "ranked_holo.csv"), index=False)

    # Top destabilisers — exclude mis_docked rows
    df_apo_clean = df_apo[~df_apo.mis_docked.fillna(True)]
    df_holo_clean = df_holo[~df_holo.mis_docked.fillna(True)]
    log(f"top apo destab (well-docked): "
        f"{df_apo_clean[['mutant','delta_vina_vs_wt','rmsd_to_native']].head(5).values.tolist()}")
    log(f"top holo destab (well-docked): "
        f"{df_holo_clean[['mutant','delta_vina_vs_wt','rmsd_to_native']].head(5).values.tolist()}")

    df_apo_clean.to_csv(os.path.join(ANA_DIR, "ranked_apo_clean.csv"), index=False)
    df_holo_clean.to_csv(os.path.join(ANA_DIR, "ranked_holo_clean.csv"), index=False)

    # Plot 1: bar chart of delta_vina apo vs holo, sorted by apo delta
    mutants_ordered = df_apo.mutant.tolist()
    apo_d, holo_d, apo_mis, holo_mis = [], [], [], []
    for m in mutants_ordered:
        a = df[(df.mutant == m) & (df.condition == "apo")]
        h = df[(df.mutant == m) & (df.condition == "holo")]
        apo_d.append(float(a.delta_vina_vs_wt.iloc[0]) if len(a) else float("nan"))
        holo_d.append(float(h.delta_vina_vs_wt.iloc[0]) if len(h) else float("nan"))
        apo_mis.append(bool(a.mis_docked.iloc[0]) if len(a) else True)
        holo_mis.append(bool(h.mis_docked.iloc[0]) if len(h) else True)

    fig, ax = plt.subplots(figsize=(14, 5))
    x = np.arange(len(mutants_ordered))
    w = 0.4
    apo_colors = ["#3498db" if not m else "#bdc3c7" for m in apo_mis]
    holo_colors = ["#e67e22" if not m else "#bdc3c7" for m in holo_mis]
    ax.bar(x - w/2, apo_d, w, color=apo_colors, label="apo (well-docked)")
    ax.bar(x + w/2, holo_d, w, color=holo_colors, label="holo (well-docked)")
    ax.bar([], [], color="#bdc3c7", label="mis-docked (RMSD>3)")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xticks(x); ax.set_xticklabels(mutants_ordered, rotation=70, ha="right", fontsize=8)
    ax.set_ylabel("Δ Vina score vs WT (kcal/mol; positive = destabilising)")
    ax.set_title("v3 mutant docking: Δ Vina score vs WT (apo vs holo)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(ANA_DIR, "delta_vina_apo_holo.png"), dpi=130)
    plt.close()

    # Plot 2: scatter apo vs holo
    fig, ax = plt.subplots(figsize=(7, 7))
    apo_arr = np.array(apo_d, dtype=float)
    holo_arr = np.array(holo_d, dtype=float)
    valid_mask = (~np.isnan(apo_arr)) & (~np.isnan(holo_arr))
    ax.scatter(apo_arr[valid_mask], holo_arr[valid_mask], s=50, alpha=0.7)
    for i, m in enumerate(mutants_ordered):
        if valid_mask[i]:
            ax.annotate(m, (apo_arr[i], holo_arr[i]), fontsize=7, alpha=0.8)
    if valid_mask.sum() > 0:
        lo = min(np.nanmin(apo_arr), np.nanmin(holo_arr)) - 0.5
        hi = max(np.nanmax(apo_arr), np.nanmax(holo_arr)) + 0.5
        ax.plot([lo, hi], [lo, hi], "k--", alpha=0.3, label="y=x")
    ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
    ax.set_xlabel("Δ Vina apo (kcal/mol)")
    ax.set_ylabel("Δ Vina holo (kcal/mol)")
    ax.set_title("Apo vs Holo Δ Vina score correlation")
    if valid_mask.sum() > 2:
        r = np.corrcoef(apo_arr[valid_mask], holo_arr[valid_mask])[0, 1]
        ax.text(0.05, 0.95, f"r = {r:.2f} (n={valid_mask.sum()})",
                transform=ax.transAxes, fontsize=10, va="top")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(ANA_DIR, "delta_vina_apo_vs_holo.png"), dpi=130)
    plt.close()

    # Plot 3: category boxplot (using clean data)
    df_clean = df[~df.mis_docked.fillna(True)]
    fig, ax = plt.subplots(figsize=(11, 5))
    cats = sorted(df_clean.category.unique())
    box_data, box_labels = [], []
    for c in cats:
        sub = df_clean[df_clean.category == c]
        box_data.append(sub.delta_vina_vs_wt.dropna().values)
        box_labels.append(f"{c}\n(n={len(sub)})")
    if box_data:
        ax.boxplot(box_data, labels=box_labels, showfliers=True)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_ylabel("Δ Vina score vs WT")
    ax.set_title("v3 Δ Vina by mutation category (well-docked only)")
    plt.xticks(rotation=30, ha="right", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(ANA_DIR, "delta_vina_by_category.png"), dpi=130)
    plt.close()

    cor = float(np.corrcoef(apo_arr[valid_mask], holo_arr[valid_mask])[0, 1]) \
          if valid_mask.sum() > 2 else None

    summary = {
        "n_rows": int(len(df)),
        "n_mutants": int(df.mutant.nunique()),
        "n_apo": int((df.condition == "apo").sum()),
        "n_holo": int((df.condition == "holo").sum()),
        "n_mis_docked": int(df.mis_docked.fillna(True).sum()),
        "wt_apo_aff": json.load(open(os.path.join(WT_DIR, "wt_apo.json")))["top_affinity"],
        "wt_holo_aff": json.load(open(os.path.join(WT_DIR, "wt_holo.json")))["top_affinity"],
        "top5_destab_apo_clean": df_apo_clean[["mutant", "category", "delta_vina_vs_wt",
                                                "rmsd_to_native", "n_modes"]].head(5).to_dict("records"),
        "top5_destab_holo_clean": df_holo_clean[["mutant", "category", "delta_vina_vs_wt",
                                                  "rmsd_to_native", "n_modes"]].head(5).to_dict("records"),
        "apo_holo_correlation": cor,
        "sign_convention": "delta_vina_vs_wt = top_aff_mut - top_aff_wt; positive = destabilising",
        "mean_topk_rule": "mean(affinities[:min(3, n_modes)])",
        "mis_docked_threshold_A": 3.0,
    }
    with open(os.path.join(ANA_DIR, "summary_v3.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)
    log(f"summary: top3 apo (well-docked) = "
        f"{[(r['mutant'], r['delta_vina_vs_wt']) for r in summary['top5_destab_apo_clean'][:3]]}")
    log("Stage 8 v3 DONE")


if __name__ == "__main__":
    main()
