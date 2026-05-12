#!/usr/bin/env python3
"""Stage 8 v2: Analysis & plots — apo vs holo, ddG ranking, vs v1 comparison."""
import os, sys, json, math
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = os.path.expanduser("~/conserved_site_project")
ANA_DIR = os.path.join(PROJECT, "08b_analysis_v2")
MUT_DIR = os.path.join(PROJECT, "07b_mut_docking_v2")
WT_DIR = os.path.join(PROJECT, "06b_docking_wt_v2")
WT_V1_DIR = os.path.join(PROJECT, "06_docking_wt")
MUT_V1_DIR = os.path.join(PROJECT, "07_mut_docking")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v2_08_analysis.log")


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V2] STAGE8: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def main():
    os.makedirs(ANA_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 8 v2 starting")

    df = pd.read_csv(os.path.join(MUT_DIR, "mutant_results_v2.csv"))
    log(f"loaded {len(df)} rows")

    # Pivot apo / holo
    pv = df.pivot_table(index=["mutant", "category"], columns="condition",
                        values=["top_affinity", "ddG_vs_wt", "rmsd_to_native"],
                        aggfunc="first").reset_index()
    pv.columns = ["_".join([str(x) for x in c if x]).strip("_") for c in pv.columns]
    pv.to_csv(os.path.join(ANA_DIR, "mutant_pivot_v2.csv"), index=False)
    log(f"wrote pivot: {pv.columns.tolist()}")

    # Ranked tables
    df_apo = df[df.condition == "apo"].sort_values("ddG_vs_wt", ascending=False).reset_index(drop=True)
    df_holo = df[df.condition == "holo"].sort_values("ddG_vs_wt", ascending=False).reset_index(drop=True)
    df_apo.to_csv(os.path.join(ANA_DIR, "ranked_apo.csv"), index=False)
    df_holo.to_csv(os.path.join(ANA_DIR, "ranked_holo.csv"), index=False)
    log(f"top apo destabilizers: {df_apo[['mutant','ddG_vs_wt']].head(5).values.tolist()}")
    log(f"top holo destabilizers: {df_holo[['mutant','ddG_vs_wt']].head(5).values.tolist()}")

    # Plot 1: bar chart of ddG, side-by-side apo/holo
    fig, ax = plt.subplots(figsize=(14, 5))
    mutants_ordered = df_apo.mutant.tolist()
    apo_d = []
    holo_d = []
    for m in mutants_ordered:
        a = df[(df.mutant == m) & (df.condition == "apo")]
        h = df[(df.mutant == m) & (df.condition == "holo")]
        apo_d.append(float(a.ddG_vs_wt.iloc[0]) if len(a) else float("nan"))
        holo_d.append(float(h.ddG_vs_wt.iloc[0]) if len(h) else float("nan"))
    x = np.arange(len(mutants_ordered))
    w = 0.4
    ax.bar(x - w/2, apo_d, w, label="apo", color="#3498db")
    ax.bar(x + w/2, holo_d, w, label="holo", color="#e67e22")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xticks(x); ax.set_xticklabels(mutants_ordered, rotation=70, ha="right", fontsize=8)
    ax.set_ylabel("ΔΔG vs WT (kcal/mol, +ve = destabilising)")
    ax.set_title("v2 mutant docking: ΔΔG vs WT (apo vs holo)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(ANA_DIR, "ddg_apo_holo.png"), dpi=130)
    plt.close()
    log("wrote ddg_apo_holo.png")

    # Plot 2: scatter apo vs holo ddG
    fig, ax = plt.subplots(figsize=(7, 7))
    apo_arr = np.array(apo_d, dtype=float)
    holo_arr = np.array(holo_d, dtype=float)
    valid_mask = (~np.isnan(apo_arr)) & (~np.isnan(holo_arr))
    ax.scatter(apo_arr[valid_mask], holo_arr[valid_mask], s=50, alpha=0.7)
    for i, m in enumerate(mutants_ordered):
        if valid_mask[i]:
            ax.annotate(m, (apo_arr[i], holo_arr[i]), fontsize=7, alpha=0.8)
    lo = min(np.nanmin(apo_arr), np.nanmin(holo_arr)) - 0.5
    hi = max(np.nanmax(apo_arr), np.nanmax(holo_arr)) + 0.5
    ax.plot([lo, hi], [lo, hi], "k--", alpha=0.3, label="y=x")
    ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
    ax.set_xlabel("ΔΔG apo (kcal/mol)")
    ax.set_ylabel("ΔΔG holo (kcal/mol)")
    ax.set_title("Apo vs Holo ΔΔG correlation")
    if valid_mask.sum() > 2:
        r = np.corrcoef(apo_arr[valid_mask], holo_arr[valid_mask])[0, 1]
        ax.text(0.05, 0.95, f"r = {r:.2f} (n={valid_mask.sum()})", transform=ax.transAxes,
               fontsize=10, va="top")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(ANA_DIR, "ddg_apo_vs_holo.png"), dpi=130)
    plt.close()
    log("wrote ddg_apo_vs_holo.png")

    # Plot 3: category boxplot
    fig, ax = plt.subplots(figsize=(11, 5))
    cats = sorted(df.category.unique())
    box_data = []
    box_labels = []
    for c in cats:
        sub = df[df.category == c]
        box_data.append(sub.ddG_vs_wt.dropna().values)
        box_labels.append(f"{c}\n(n={len(sub)})")
    ax.boxplot(box_data, labels=box_labels, showfliers=True)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_ylabel("ΔΔG vs WT")
    ax.set_title("v2 ΔΔG by mutation category (apo+holo combined)")
    plt.xticks(rotation=30, ha="right", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(ANA_DIR, "ddg_by_category.png"), dpi=130)
    plt.close()
    log("wrote ddg_by_category.png")

    # Plot 4: v1 vs v2 comparison (only on shared mutants)
    v1_csv = os.path.join(MUT_V1_DIR, "mutant_results.csv")
    v1_df = None
    if os.path.exists(v1_csv):
        try:
            v1_df = pd.read_csv(v1_csv)
            log(f"loaded v1 mutant CSV: {len(v1_df)} rows, cols={v1_df.columns.tolist()}")
        except Exception as e:
            log(f"could not load v1 CSV: {e}")

    if v1_df is not None and "mutant" in v1_df.columns:
        # Find common mutants
        v2_apo = df[df.condition == "apo"].set_index("mutant")
        v1_top_col = None
        for cand in ["top_affinity", "vina_top", "best_affinity"]:
            if cand in v1_df.columns:
                v1_top_col = cand
                break
        if v1_top_col is None and "ddG_vs_wt" in v1_df.columns:
            v1_top_col = "ddG_vs_wt"
        if v1_top_col:
            v1_set = v1_df.set_index("mutant")
            common = set(v2_apo.index) & set(v1_set.index)
            log(f"common mutants v1∩v2: {len(common)}")
            if common:
                fig, ax = plt.subplots(figsize=(7, 7))
                xs, ys, lbls = [], [], []
                for m in sorted(common):
                    if "ddG_vs_wt" in v1_df.columns and "ddG_vs_wt" in v2_apo.columns:
                        xs.append(float(v1_set.loc[m, "ddG_vs_wt"]))
                        ys.append(float(v2_apo.loc[m, "ddG_vs_wt"]))
                        lbls.append(m)
                if xs:
                    ax.scatter(xs, ys, s=50)
                    for i, m in enumerate(lbls):
                        ax.annotate(m, (xs[i], ys[i]), fontsize=7)
                    lo, hi = min(xs+ys) - 0.5, max(xs+ys) + 0.5
                    ax.plot([lo, hi], [lo, hi], "k--", alpha=0.3)
                    ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
                    ax.set_xlabel("v1 ΔΔG (chain A only, exh=16, box=22)")
                    ax.set_ylabel("v2 ΔΔG apo (dimer, exh=32, box=18)")
                    ax.set_title("v1 vs v2 ΔΔG (common mutants)")
                    if len(xs) > 2:
                        r = np.corrcoef(xs, ys)[0, 1]
                        ax.text(0.05, 0.95, f"r = {r:.2f}", transform=ax.transAxes, fontsize=10, va="top")
                    plt.tight_layout()
                    plt.savefig(os.path.join(ANA_DIR, "v1_vs_v2_ddg.png"), dpi=130)
                    plt.close()
                    log("wrote v1_vs_v2_ddg.png")

    # Stats summary
    summary = {
        "n_mutants_total": int(len(df) // 2) if df.condition.nunique() == 2 else int(len(df.mutant.unique())),
        "n_apo": int((df.condition == "apo").sum()),
        "n_holo": int((df.condition == "holo").sum()),
        "wt_apo_aff": json.load(open(os.path.join(WT_DIR, "wt_apo.json")))["top_affinity"],
        "wt_holo_aff": json.load(open(os.path.join(WT_DIR, "wt_holo.json")))["top_affinity"],
        "top5_destabilising_apo": df_apo[["mutant","ddG_vs_wt"]].head(5).to_dict("records"),
        "top5_destabilising_holo": df_holo[["mutant","ddG_vs_wt"]].head(5).to_dict("records"),
        "apo_holo_correlation": float(np.corrcoef(apo_arr[valid_mask], holo_arr[valid_mask])[0, 1])
                                 if valid_mask.sum() > 2 else None,
    }
    with open(os.path.join(ANA_DIR, "summary_v2.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)
    log(f"summary: {summary}")

    log("Stage 8 v2 DONE")


if __name__ == "__main__":
    main()
