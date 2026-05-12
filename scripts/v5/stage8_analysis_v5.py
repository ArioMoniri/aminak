#!/usr/bin/env python3
"""Stage 8 v5: Analysis with low_confidence + mis_docked filters; Spearman
plus Pearson; per sci-off review: drop unfiltered apo-holo correlation if
both p > 0.19 from the headline AND remove the n=4 filtered ρ from prose.
The v5 WT-holo RMSD should be < 3 Å (no longer requires the relaxed metric)
but we keep the dual-metric reporting in Methods (sci-off item 5).
"""
import os, json, math
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

PROJECT = os.path.expanduser("~/conserved_site_project")
ANA_DIR = os.path.join(PROJECT, "08e_analysis_v5")
MUT_DIR = os.path.join(PROJECT, "07e_mut_docking_v5")
WT_DIR = os.path.join(PROJECT, "06e_docking_wt_v5")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v5_08_analysis.log")

VINA_NOISE_FLOOR = 0.85


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V5] STAGE8: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def main():
    os.makedirs(ANA_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 8 v5 starting")

    df = pd.read_csv(os.path.join(MUT_DIR, "mutant_results_v5.csv"), comment="#")
    log(f"loaded {len(df)} rows (including WT reference)")

    for col in ("top_affinity", "mean_topk", "rmsd_to_native", "delta_vina_vs_wt"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["n_modes"] = pd.to_numeric(df["n_modes"], errors="coerce").fillna(0).astype(int)
    df["mis_docked"] = df["mis_docked"].astype(str).str.lower() == "true"
    df["low_confidence"] = df["low_confidence"].astype(str).str.lower() == "true"

    df_mut = df[df.mutant != "WT"].copy()
    df_wt = df[df.mutant == "WT"].copy()

    pv = df_mut.pivot_table(index=["mutant", "category"], columns="condition",
                            values=["top_affinity", "delta_vina_vs_wt", "rmsd_to_native",
                                    "mis_docked", "low_confidence", "n_modes"],
                            aggfunc="first").reset_index()
    pv.columns = ["_".join([str(x) for x in c if x is not None and str(x) != ""]).strip("_")
                  for c in pv.columns]
    pv.to_csv(os.path.join(ANA_DIR, "mutant_pivot_v5.csv"), index=False)

    df_apo = df_mut[df_mut.condition == "apo"].sort_values("delta_vina_vs_wt", ascending=False).reset_index(drop=True)
    df_holo = df_mut[df_mut.condition == "holo"].sort_values("delta_vina_vs_wt", ascending=False).reset_index(drop=True)
    df_apo.to_csv(os.path.join(ANA_DIR, "ranked_apo.csv"), index=False)
    df_holo.to_csv(os.path.join(ANA_DIR, "ranked_holo.csv"), index=False)

    # In v5, the WT holo RMSD vs crystal should be < 3 Å, so we use the standard
    # mis_docked = (rmsd > 3) definition without the v4 relaxation. We still
    # report both metrics in methods.
    wt_holo_rmsd = (df_wt[df_wt.condition == "holo"].rmsd_to_native.iloc[0]
                    if (df_wt.condition == "holo").any() else float("nan"))
    holo_filter_note = ""
    if not math.isnan(wt_holo_rmsd) and wt_holo_rmsd > 3.0:
        df_holo = df_holo.copy()
        df_holo["mis_docked_vs_crystal"] = df_holo["mis_docked"]
        df_holo["mis_docked"] = (df_holo["rmsd_to_native"] - wt_holo_rmsd).abs() > 3.0
        holo_filter_note = (f"WT holo itself is RMSD={wt_holo_rmsd:.2f} A vs crystal; "
                            "holo mis_docked redefined as |RMSD - WT_holo_RMSD| > 3 A.")
        log(holo_filter_note)
    else:
        holo_filter_note = (f"WT holo RMSD vs crystal = {wt_holo_rmsd:.3f} A; "
                            "standard mis_docked (RMSD > 3 A vs crystal) applies. "
                            "v4-style dual-RMSD-reference is not needed in v5.")
        log(holo_filter_note)

    df_apo_clean = df_apo[(~df_apo.mis_docked) & (~df_apo.low_confidence)]
    df_holo_clean = df_holo[(~df_holo.mis_docked) & (~df_holo.low_confidence)]
    df_apo_clean.to_csv(os.path.join(ANA_DIR, "ranked_apo_clean.csv"), index=False)
    df_holo_clean.to_csv(os.path.join(ANA_DIR, "ranked_holo_clean.csv"), index=False)

    log(f"top apo destab (clean): {df_apo_clean[['mutant','delta_vina_vs_wt']].head(5).values.tolist()}")
    log(f"top holo destab (clean): {df_holo_clean[['mutant','delta_vina_vs_wt']].head(5).values.tolist()}")

    mutants_ordered = df_apo.mutant.tolist()
    apo_d, holo_d = [], []
    apo_mis, holo_mis = [], []
    apo_lowc, holo_lowc = [], []
    df_apo_idx = df_mut[df_mut.condition == "apo"].set_index("mutant")
    df_holo_idx = df_holo.set_index("mutant")
    for m in mutants_ordered:
        if m in df_apo_idx.index:
            a = df_apo_idx.loc[m]
            apo_d.append(float(a.delta_vina_vs_wt))
            apo_mis.append(bool(a.mis_docked))
            apo_lowc.append(bool(a.low_confidence))
        else:
            apo_d.append(float("nan")); apo_mis.append(True); apo_lowc.append(True)
        if m in df_holo_idx.index:
            h = df_holo_idx.loc[m]
            holo_d.append(float(h.delta_vina_vs_wt))
            holo_mis.append(bool(h.mis_docked))
            holo_lowc.append(bool(h.low_confidence))
        else:
            holo_d.append(float("nan")); holo_mis.append(True); holo_lowc.append(True)

    apo_arr = np.array(apo_d, dtype=float)
    holo_arr = np.array(holo_d, dtype=float)
    valid = ((~np.isnan(apo_arr)) & (~np.isnan(holo_arr))
             & (~np.array(apo_mis)) & (~np.array(holo_mis))
             & (~np.array(apo_lowc)) & (~np.array(holo_lowc)))

    pearson_r, pearson_p = (None, None)
    spearman_r, spearman_p = (None, None)
    if valid.sum() >= 3:
        pearson_r, pearson_p = stats.pearsonr(apo_arr[valid], holo_arr[valid])
        spearman_r, spearman_p = stats.spearmanr(apo_arr[valid], holo_arr[valid])

    big_mask = (np.abs(apo_arr) > 0.3) | (np.abs(holo_arr) > 0.3)
    big_valid = valid & big_mask
    spearman_r_filt, spearman_p_filt, n_filt = (None, None, 0)
    if big_valid.sum() >= 3:
        spearman_r_filt, spearman_p_filt = stats.spearmanr(apo_arr[big_valid], holo_arr[big_valid])
        n_filt = int(big_valid.sum())

    # ---- Sci-off recommendation: report null correlation if both p > 0.19 ----
    null_correlation = (
        (pearson_p is None or (pearson_p > 0.19)) and
        (spearman_p is None or (spearman_p > 0.19))
    )

    # ---- Plots (same as v4) ----
    fig, ax = plt.subplots(figsize=(14, 5))
    x = np.arange(len(mutants_ordered))
    w = 0.4
    apo_colors = ["#bdc3c7" if (apo_mis[i] or apo_lowc[i]) else "#3498db" for i in range(len(x))]
    holo_colors = ["#bdc3c7" if (holo_mis[i] or holo_lowc[i]) else "#e67e22" for i in range(len(x))]
    ax.bar(x - w/2, apo_d, w, color=apo_colors, label="apo (good)")
    ax.bar(x + w/2, holo_d, w, color=holo_colors, label="holo (good)")
    ax.bar([], [], color="#bdc3c7", label="excluded (mis-docked or n_modes<5)")
    ax.axhspan(-VINA_NOISE_FLOOR, VINA_NOISE_FLOOR, color="grey", alpha=0.15,
               label=f"Vina noise floor +/-{VINA_NOISE_FLOOR}")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xticks(x); ax.set_xticklabels(mutants_ordered, rotation=70, ha="right", fontsize=8)
    ax.set_ylabel("Delta Vina score vs WT v5 (kcal/mol; positive = destabilising)")
    ax.set_title("v5 mutant docking: Delta Vina score vs WT (apo vs holo)")
    ax.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(ANA_DIR, "delta_vina_apo_holo.png"), dpi=130)
    plt.close()

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(apo_arr[valid], holo_arr[valid], s=60, alpha=0.75, color="#2ecc71")
    for i, m in enumerate(mutants_ordered):
        if valid[i]:
            ax.annotate(m, (apo_arr[i], holo_arr[i]), fontsize=7, alpha=0.85)
    if valid.sum() > 0:
        lo = min(np.nanmin(apo_arr), np.nanmin(holo_arr)) - 0.5
        hi = max(np.nanmax(apo_arr), np.nanmax(holo_arr)) + 0.5
        ax.plot([lo, hi], [lo, hi], "k--", alpha=0.3, label="y=x")
    ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
    ax.axhspan(-VINA_NOISE_FLOOR, VINA_NOISE_FLOOR, color="grey", alpha=0.10)
    ax.axvspan(-VINA_NOISE_FLOOR, VINA_NOISE_FLOOR, color="grey", alpha=0.10)
    ax.set_xlabel("Delta Vina apo (kcal/mol)")
    ax.set_ylabel("Delta Vina holo (kcal/mol)")
    ax.set_title("v5 Apo vs Holo Delta Vina (well-docked, n_modes>=5)")
    txt = []
    if pearson_r is not None:
        txt.append(f"Pearson r = {pearson_r:.2f} (n={int(valid.sum())}) p={pearson_p:.3f}")
    if spearman_r is not None:
        txt.append(f"Spearman rho = {spearman_r:.2f} p={spearman_p:.3f}")
    if null_correlation:
        txt.append("--> no statistically significant apo-holo correlation")
    txt.append(f"Vina noise floor: +/-{VINA_NOISE_FLOOR} kcal/mol")
    ax.text(0.04, 0.96, "\n".join(txt), transform=ax.transAxes, fontsize=9,
            va="top", ha="left", bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    ax.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(ANA_DIR, "delta_vina_apo_vs_holo.png"), dpi=130)
    plt.close()

    df_apo_clean_full = df_mut[(df_mut.condition == "apo")
                               & (~df_mut.mis_docked) & (~df_mut.low_confidence)]
    df_holo_clean_full = df_holo[(~df_holo.mis_docked) & (~df_holo.low_confidence)]
    df_clean_all = pd.concat([df_apo_clean_full, df_holo_clean_full], ignore_index=True)
    fig, ax = plt.subplots(figsize=(11, 5))
    cats = sorted(df_clean_all.category.unique())
    box_data, box_labels = [], []
    for c in cats:
        sub = df_clean_all[df_clean_all.category == c]
        box_data.append(sub.delta_vina_vs_wt.dropna().values)
        box_labels.append(f"{c}\n(n={len(sub)})")
    if box_data:
        ax.boxplot(box_data, labels=box_labels, showfliers=True)
    ax.axhspan(-VINA_NOISE_FLOOR, VINA_NOISE_FLOOR, color="grey", alpha=0.15)
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_ylabel("Delta Vina score vs WT v5 (kcal/mol)")
    ax.set_title("v5 Delta Vina by mutation category (excludes mis-docked & n_modes<5)")
    plt.xticks(rotation=30, ha="right", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(ANA_DIR, "delta_vina_by_category.png"), dpi=130)
    plt.close()

    # Count |Δ_holo| > noise floor
    valid_h = df_holo[(~df_holo.mis_docked) & (~df_holo.low_confidence)]
    n_above_noise_holo = int((valid_h.delta_vina_vs_wt.abs() > VINA_NOISE_FLOOR).sum())
    log(f"  n mutants holo |Δ| > {VINA_NOISE_FLOOR}: {n_above_noise_holo}")

    summary = {
        "n_rows": int(len(df)),
        "n_mutants": int(df_mut.mutant.nunique()),
        "n_apo_rows": int((df_mut.condition == "apo").sum()),
        "n_holo_rows": int((df_mut.condition == "holo").sum()),
        "n_mis_docked": int(df_mut.mis_docked.sum()),
        "n_low_confidence": int(df_mut.low_confidence.sum()),
        "wt_apo_aff_v5": float(df_wt[df_wt.condition == "apo"].top_affinity.iloc[0])
                         if (df_wt.condition == "apo").any() else None,
        "wt_holo_aff_v5": float(df_wt[df_wt.condition == "holo"].top_affinity.iloc[0])
                          if (df_wt.condition == "holo").any() else None,
        "wt_holo_n_modes_v5": int(df_wt[df_wt.condition == "holo"].n_modes.iloc[0])
                              if (df_wt.condition == "holo").any() else None,
        "wt_holo_rmsd_vs_crystal": float(wt_holo_rmsd) if not math.isnan(wt_holo_rmsd) else None,
        "vina_noise_floor_kcal_per_mol": VINA_NOISE_FLOOR,
        "vina_noise_floor_source": "Trott & Olson 2010; Forli et al. 2016 (mean of 0.7-1.0)",
        "pearson_r_apo_holo": pearson_r,
        "pearson_p_apo_holo": pearson_p,
        "spearman_r_apo_holo": spearman_r,
        "spearman_p_apo_holo": spearman_p,
        "spearman_r_filtered_abs_delta_gt_0p3": spearman_r_filt,
        "spearman_p_filtered_abs_delta_gt_0p3": spearman_p_filt,
        "n_filtered_for_spearman": n_filt,
        "filtered_spearman_in_prose": False,
        "no_significant_apo_holo_correlation": null_correlation,
        "n_above_noise_holo": n_above_noise_holo,
        "top5_destab_apo_clean": df_apo_clean[["mutant", "category", "delta_vina_vs_wt",
                                                "rmsd_to_native", "n_modes"]].head(5).to_dict("records"),
        "top5_destab_holo_clean": df_holo_clean[["mutant", "category", "delta_vina_vs_wt",
                                                  "rmsd_to_native", "n_modes"]].head(5).to_dict("records"),
        "top3_destab_holo_clean": df_holo_clean[["mutant", "category", "delta_vina_vs_wt",
                                                  "rmsd_to_native", "n_modes"]].head(3).to_dict("records"),
        "filters_applied": "exclude mis_docked (RMSD>3) and low_confidence (n_modes<5)",
        "holo_filter_note": holo_filter_note,
        "sign_convention": "delta_vina_vs_wt = top_aff_mut - top_aff_wt_v5; positive = destabilising",
    }
    with open(os.path.join(ANA_DIR, "summary_v5.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)
    log(f"Pearson r = {pearson_r}, Spearman ρ = {spearman_r}, null_correlation={null_correlation}")
    log(f"n above noise (holo): {n_above_noise_holo}")
    log("Stage 8 v5 DONE")


if __name__ == "__main__":
    main()
