#!/usr/bin/env python3
"""Phase 14 — aggregate per-strategy CSVs into master.csv + generate the 6 headline figures."""
from __future__ import annotations
import csv, json, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
PHASE = REPO / "14_inhibitor_design"
FIG = PHASE / "figures"
FIG.mkdir(parents=True, exist_ok=True)
AGG = PHASE / "05_aggregate"
AGG.mkdir(parents=True, exist_ok=True)

NOISE = 0.85   # Vina noise floor (Trott & Olson 2010)


def load_csv(p: Path):
    if not p.exists(): return []
    return list(csv.DictReader(p.open()))


def to_float(x):
    try: return float(x) if x not in (None, "") else None
    except (ValueError, TypeError): return None


def main():
    s1 = load_csv(PHASE / "01_active_site" / "results_summary.csv")
    s2 = load_csv(PHASE / "02_cofactor_site" / "results_summary.csv")
    s3 = load_csv(PHASE / "03_dimer_interface" / "results_summary.csv")
    s4 = load_csv(PHASE / "04_allosteric" / "results_summary.csv")
    print(f"loaded: S1={len(s1)} S2={len(s2)} S3={len(s3)} S4={len(s4)}")

    # tag each row with its strategy
    for r in s1: r["strategy"] = "1_active_site"
    for r in s2: r["strategy"] = "2_cofactor_site"
    for r in s3: r["strategy"] = "3_dimer_interface"
    for r in s4: r["strategy"] = "4_allosteric"

    # union of columns
    all_rows = s1 + s2 + s3 + s4
    if not all_rows:
        print("!!! no data in any strategy"); return
    all_cols = set()
    for r in all_rows: all_cols.update(r.keys())
    all_cols = sorted(all_cols)
    master = AGG / "master.csv"
    with master.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=all_cols)
        w.writeheader(); w.writerows(all_rows)
    print(f"  → {master}")

    # ===== FIGURE 1: per-strategy distribution of top1 (violin/strip) =====
    fig, ax = plt.subplots(1, 1, figsize=(9, 5))
    data, labels, refs = [], [], {}
    for tag, rows, ref_col in [
        ("S1 active-site\n(apo)", [r for r in s1 if r.get("state")=="apo" and to_float(r.get("top1_mean")) is not None], "delta_vs_dump"),
        ("S1 active-site\n(holo)", [r for r in s1 if r.get("state")=="holo" and to_float(r.get("top1_mean")) is not None], "delta_vs_dump"),
        ("S2 cofactor\n(apo)", [r for r in s2 if r.get("state")=="apo" and to_float(r.get("top1_mean")) is not None], "delta_vs_raltitrexed"),
        ("S2 cofactor\n(holo)", [r for r in s2 if r.get("state")=="holo" and to_float(r.get("top1_mean")) is not None], "delta_vs_raltitrexed"),
        ("S3 dimer-iface\n(peptides)", [r for r in s3 if to_float(r.get("top1_mean")) is not None], None),
        ("S4 allosteric\n(fragments)", [r for r in s4 if to_float(r.get("top1") or r.get("top1_mean")) is not None], None),
    ]:
        vals = [to_float(r.get("top1_mean") or r.get("top1")) for r in rows]
        vals = [v for v in vals if v is not None]
        if not vals: continue
        data.append(vals); labels.append(f"{tag}\n(n={len(vals)})")
    if data:
        parts = ax.violinplot(data, showmeans=False, showmedians=True)
        ax.scatter([i+1 for i, ds in enumerate(data) for _ in ds],
                   [v for ds in data for v in ds], alpha=0.5, s=18, color="navy")
        ax.set_xticks(range(1, len(labels)+1))
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel("Vina top1 (kcal/mol)  — more negative = tighter")
        ax.set_title("Phase 14: per-strategy distribution of best-pose Vina scores")
        ax.axhline(y=-8.78, ls="--", color="red", alpha=0.5, label="dUMP apo reference (-8.78)")
        ax.axhline(y=-9.08, ls="--", color="green", alpha=0.5, label="raltitrexed apo reference (-9.08)")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(True, alpha=0.3)
        plt.tight_layout(); plt.savefig(FIG / "fig1_distributions.png", dpi=140); plt.close()
        print(f"  → {FIG/'fig1_distributions.png'}")

    # ===== FIGURE 2: Δ vs strategy reference (ranked bar) =====
    # only meaningful for S1 (Δ vs dUMP) and S2 (Δ vs raltitrexed)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9))
    for ax, rows, dcol, refname in [
        (ax1, [r for r in s1 if r.get("state")=="apo"], "delta_vs_dump", "dUMP"),
        (ax2, [r for r in s2 if r.get("state")=="apo"], "delta_vs_raltitrexed", "raltitrexed"),
    ]:
        rows = [r for r in rows if to_float(r.get(dcol)) is not None]
        rows.sort(key=lambda r: to_float(r[dcol]))
        names = [r["compound"] for r in rows]
        deltas = [to_float(r[dcol]) for r in rows]
        colors = ["green" if d <= -NOISE else ("red" if d >= NOISE else "grey") for d in deltas]
        ax.barh(range(len(names)), deltas, color=colors, alpha=0.8)
        ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=9)
        ax.axvline(x=0, color="k", linewidth=0.5)
        ax.axvline(x=-NOISE, ls=":", color="green", alpha=0.5, label=f"Δ ≤ -{NOISE} = beats ref ({NOISE} = Vina noise floor)")
        ax.axvline(x=+NOISE, ls=":", color="red", alpha=0.5, label=f"Δ ≥ +{NOISE} = significantly worse")
        ax.set_xlabel(f"Δ vs {refname} (kcal/mol)  — negative = tighter")
        ax.set_title(f"Apo Δ vs {refname} (strategy reference)")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(True, axis="x", alpha=0.3)
    plt.tight_layout(); plt.savefig(FIG / "fig2_delta_ranking.png", dpi=140); plt.close()
    print(f"  → {FIG/'fig2_delta_ranking.png'}")

    # ===== FIGURE 3: apo-vs-holo gap (cryptic-pocket indicator) =====
    fig, ax = plt.subplots(1, 1, figsize=(9, 5))
    for strat, rows, color in [("S1 active-site", s1, "blue"), ("S2 cofactor-site", s2, "orange")]:
        compounds = sorted(set(r["compound"] for r in rows))
        apo_holo = []
        for c in compounds:
            apo = next((to_float(r["top1_mean"]) for r in rows if r["compound"] == c and r.get("state")=="apo" and to_float(r["top1_mean"]) is not None), None)
            holo = next((to_float(r["top1_mean"]) for r in rows if r["compound"] == c and r.get("state")=="holo" and to_float(r["top1_mean"]) is not None), None)
            if apo is None or holo is None: continue
            apo_holo.append((c, apo, holo, apo - holo))
        if not apo_holo: continue
        names = [x[0] for x in apo_holo]
        gaps = [x[3] for x in apo_holo]
        x = np.arange(len(names))
        ax.bar(x + (-0.2 if strat=="S1 active-site" else 0.2), gaps, width=0.4, color=color, alpha=0.7, label=strat)
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=45, ha="right", fontsize=8)
    ax.axhline(y=0, color="k", linewidth=0.5)
    ax.axhline(y=-2, ls=":", color="purple", alpha=0.5, label="cryptic-pocket threshold (apo < holo by ≥ 2 kcal/mol)")
    ax.set_ylabel("apo − holo top1 (kcal/mol)\nnegative = apo binds tighter (cofactor competes)")
    ax.set_title("Apo-minus-holo gap per compound (cryptic-pocket / induced-fit indicator)")
    ax.legend(loc="best", fontsize=8); ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout(); plt.savefig(FIG / "fig3_apo_holo_gap.png", dpi=140); plt.close()
    print(f"  → {FIG/'fig3_apo_holo_gap.png'}")

    # ===== FIGURE 4: Tier-1 vs Tier-2 separation (enrichment proxy) =====
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, rows, title in [(axes[0], [r for r in s1 if r.get("state")=="apo"], "S1 active-site (apo)"),
                            (axes[1], [r for r in s2 if r.get("state")=="apo"], "S2 cofactor-site (apo)")]:
        t1 = sorted([to_float(r["top1_mean"]) for r in rows if str(r.get("tier"))=="1" and to_float(r["top1_mean"]) is not None])
        t2 = sorted([to_float(r["top1_mean"]) for r in rows if str(r.get("tier"))=="2" and to_float(r["top1_mean"]) is not None])
        if not t1 or not t2:
            ax.text(0.5, 0.5, f"n_T1={len(t1)} n_T2={len(t2)} — insufficient data", ha="center", va="center")
            ax.set_title(title); continue
        ax.boxplot([t1, t2], labels=[f"Tier 1 (actives)\nn={len(t1)}", f"Tier 2 (decoys)\nn={len(t2)}"], widths=0.5)
        for x, vals in enumerate([t1, t2]):
            ax.scatter([x+1 + (np.random.rand()-0.5)*0.1 for _ in vals], vals, alpha=0.6, color="navy" if x==0 else "darkred", s=25)
        ax.set_ylabel("Vina top1 (kcal/mol)")
        ax.set_title(f"{title}\nΔ(median) = {np.median(t1)-np.median(t2):+.2f}")
        ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout(); plt.savefig(FIG / "fig4_tier_separation.png", dpi=140); plt.close()
    print(f"  → {FIG/'fig4_tier_separation.png'}")

    print("\nDone.")

if __name__ == "__main__":
    main()
