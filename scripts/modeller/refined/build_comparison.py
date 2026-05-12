"""Build the comparative plots:
1. comparison_before_after.png — bar chart of %fav/%all/%out for:
   - 1HVY crystal (reference)
   - 10 baseline models (mean ± SD)
   - 10 refined models (mean ± SD)
   - best loop-refined model
2. outlier_position_map.png — per-residue heat-map of outlier counts
   before vs after refinement, indexed by P04818 residue number
   (model residue + 25).
"""
from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT_DIR = Path(os.environ.get("PROJECT_DIR",
                                  str(Path.home() / "conserved_site_project")))
BASE_DIR = PROJECT_DIR / "10b_modeller_refined/01_baseline_lovell"
REFINED_DIR = PROJECT_DIR / "10b_modeller_refined/04_refined_lovell"


def load_summary(csv_path: Path) -> List[Dict]:
    with open(csv_path) as fh:
        return list(csv.DictReader(fh))


def _stats_set(rows, labels_must_contain=None, label_in=None):
    fav, al, out = [], [], []
    for r in rows:
        lab = r["label"]
        if labels_must_contain and labels_must_contain not in lab:
            continue
        if label_in is not None and lab not in label_in:
            continue
        fav.append(float(r["favoured_pct"]))
        al.append(float(r["allowed_pct"]))
        out.append(float(r["outlier_pct"]))
    return np.array(fav), np.array(al), np.array(out)


def build_bar_chart(out_path: Path):
    base = load_summary(BASE_DIR / "summary.csv")
    refined = load_summary(REFINED_DIR / "summary.csv")

    # 1HVY crystal
    crystal = [r for r in base if "1HVY" in r["label"]]
    crystal_fav = float(crystal[0]["favoured_pct"]) if crystal else 0.0
    crystal_all = float(crystal[0]["allowed_pct"]) if crystal else 0.0
    crystal_out = float(crystal[0]["outlier_pct"]) if crystal else 0.0

    # 10 baseline models (exclude best_model.pdb duplicate and crystal)
    base_models = [r for r in base
                   if r["label"].startswith("target.B999")]
    bf, ba, bo = _stats_set(base_models)

    # Refined models (filenames refined_B99990001 ... etc.)
    refined_models = [r for r in refined
                      if r["label"].startswith("refined_B")]
    rf, ra, ro = _stats_set(refined_models)

    # Best loop-refined
    loop = [r for r in refined if "best_loop_refined" in r["label"]]
    if loop:
        lf = float(loop[0]["favoured_pct"])
        la = float(loop[0]["allowed_pct"])
        lo = float(loop[0]["outlier_pct"])
    else:
        lf = la = lo = 0.0

    groups = ["1HVY\ncrystal", "Baseline\n(n=10)",
              "Refined\n(very_slow, n=10)", "Best\nloop-refined"]
    fav_means = [crystal_fav, bf.mean() if len(bf) else 0,
                 rf.mean() if len(rf) else 0, lf]
    all_means = [crystal_all, ba.mean() if len(ba) else 0,
                 ra.mean() if len(ra) else 0, la]
    out_means = [crystal_out, bo.mean() if len(bo) else 0,
                 ro.mean() if len(ro) else 0, lo]
    fav_sd = [0, bf.std() if len(bf) else 0,
              rf.std() if len(rf) else 0, 0]
    all_sd = [0, ba.std() if len(ba) else 0,
              ra.std() if len(ra) else 0, 0]
    out_sd = [0, bo.std() if len(bo) else 0,
              ro.std() if len(ro) else 0, 0]

    x = np.arange(len(groups))
    width = 0.27

    fig, ax = plt.subplots(figsize=(10, 6))
    b1 = ax.bar(x - width, fav_means, width, yerr=fav_sd, label="Favoured",
                color="#198754", capsize=4)
    b2 = ax.bar(x, all_means, width, yerr=all_sd, label="Allowed",
                color="#fd7e14", capsize=4)
    b3 = ax.bar(x + width, out_means, width, yerr=out_sd, label="Outlier",
                color="#dc3545", capsize=4)

    for bars, vals in [(b1, fav_means), (b2, all_means), (b3, out_means)]:
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                    f"{v:.1f}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(groups)
    ax.set_ylabel("% of residues")
    ax.set_title("Ramachandran statistics (Lovell scheme): "
                 "before vs after Modeller refinement")
    ax.set_ylim(0, 105)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"Wrote {out_path}", flush=True)

    return {
        "crystal": {"fav": crystal_fav, "all": crystal_all, "out": crystal_out},
        "baseline_mean": {"fav": bf.mean() if len(bf) else 0,
                          "all": ba.mean() if len(ba) else 0,
                          "out": bo.mean() if len(bo) else 0,
                          "fav_sd": bf.std() if len(bf) else 0,
                          "out_sd": bo.std() if len(bo) else 0,
                          "n": len(bf)},
        "refined_mean": {"fav": rf.mean() if len(rf) else 0,
                         "all": ra.mean() if len(ra) else 0,
                         "out": ro.mean() if len(ro) else 0,
                         "fav_sd": rf.std() if len(rf) else 0,
                         "out_sd": ro.std() if len(ro) else 0,
                         "n": len(rf)},
        "loop": {"fav": lf, "all": la, "out": lo},
    }


def build_outlier_heatmap(out_path: Path):
    """Per-residue outlier count map.

    For each P04818 residue number (model resnum + 25), count how many of
    the 10 baseline models flagged that residue as outlier, and how many
    of the 10 refined models did. Display as two stacked rows.
    """
    base_dir = BASE_DIR
    refined_dir = REFINED_DIR

    def collect(dir_, prefix):
        # Per-residue outlier counts across all matching CSVs
        counts: Dict[int, int] = {}
        files = sorted(dir_.glob(f"lovell_stats_{prefix}*.csv"))
        for f in files:
            with open(f) as fh:
                for row in csv.DictReader(fh):
                    if row["classification"] == "outlier":
                        rn = int(row["resnum"])
                        counts[rn] = counts.get(rn, 0) + 1
        return counts, len(files)

    base_counts, n_base = collect(base_dir, "target.B999")
    ref_counts, n_ref = collect(refined_dir, "refined_B")

    # x-axis: model residue range (1..287)
    if not base_counts and not ref_counts:
        print("No outliers in either set — skipping heatmap", flush=True)
        return
    all_res = sorted(set(list(base_counts.keys()) + list(ref_counts.keys())))
    lo, hi = min(all_res), max(all_res)
    span = list(range(lo, hi + 1))
    base_arr = np.array([base_counts.get(r, 0) for r in span])
    ref_arr = np.array([ref_counts.get(r, 0) for r in span])

    fig, ax = plt.subplots(figsize=(14, 3.2))
    data = np.vstack([base_arr, ref_arr])
    im = ax.imshow(data, aspect="auto", cmap="Reds", interpolation="nearest",
                   vmin=0, vmax=max(1, data.max()))
    ax.set_yticks([0, 1])
    ax.set_yticklabels([f"Baseline\n(n={n_base})",
                        f"Refined\n(n={n_ref})"])
    # x ticks every 20 residues, labeled with both model# and P04818
    n = len(span)
    tick_every = max(1, n // 20)
    xticks = list(range(0, n, tick_every))
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{span[i]}\n({span[i]+25})" for i in xticks],
                       fontsize=7)
    ax.set_xlabel("Model residue # (P04818 residue # in parentheses)")
    ax.set_title("Outlier positions — counts across 10 models, "
                 "baseline vs refined")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label("# models with outlier here")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f"Wrote {out_path}", flush=True)


def main() -> int:
    REFINED_DIR.mkdir(parents=True, exist_ok=True)
    stats = build_bar_chart(REFINED_DIR / "comparison_before_after.png")
    build_outlier_heatmap(REFINED_DIR / "outlier_position_map.png")

    # Final numeric summary
    import json
    (REFINED_DIR / "comparison_summary.json").write_text(
        json.dumps(stats, indent=2, default=float))
    print(json.dumps(stats, indent=2, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
