#!/usr/bin/env python3
"""Phase 14e — Smina rescoring with electrostatics + desolvation.

Phase 7-8 documented that rigid Vina cannot resolve TYMS mutants at the kcal scale,
in particular the charge-reversal mutants (R215E, R175E_R176E) where the Arg(+1) →
Glu(−1) flip should impose a heavy electrostatic penalty. Vina's scoring is
electrostatics-free, so the penalty is invisible.

Smina is a Vina fork (Koes 2013) that adds proper desolvation + Coulomb terms and
supports custom-weighted scoring functions. We rescore the existing Phase 7-8 top
poses with three scorers:

  1. vina      — Smina's reimplementation of Vina (sanity check, should match)
  2. vinardo   — Quiroga & Villarreal 2016 (stiffer hydrophobic, repulsive)
  3. custom-q  — Vina + ad4_solvation + electrostatic terms enabled, weighted to
                 expose the Arg→Glu cost

The custom scoring file appends an explicit electrostatic (Coulomb-style) and
ad4_solvation term on top of the Vina defaults — Smina's --custom_scoring takes
a file listing one term per line with the weight in front.

Outputs:
  14_inhibitor_design/06_smina_rescore/
    custom_scoring_q.txt      — the custom scoring definition (committed)
    rescore_results.csv       — long table: pose × scorer × score
    rescore_summary.csv       — wide pivot: pose | vina | vinardo | custom_q | Δ
    rescore_plot.png          — bar chart vs WT for each scorer
"""
from __future__ import annotations
import subprocess, csv, json, re, shutil
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "14_inhibitor_design" / "06_smina_rescore"
OUT.mkdir(parents=True, exist_ok=True)
SMINA = shutil.which("smina") or "/opt/homebrew/bin/smina"

# Receptor — the Phase-6c hardened apo dimer (Vina-compatible, AMBER ff14SB charges)
APO_RECEPTOR = REPO / "06f_receptor_fixed" / "protein_dimer_apo_fixed.pdbqt"

# The Phase 7-8 "stripped" mutant top poses (one PDBQT per mutant, MODEL 1 only)
# We use the holo dock (cofactor present) because that's the canonical Phase-7 setting
# and where the R215E / R175E_R176E signs are documented to be wrong.
STRIPPED = REPO / "13_phase8" / "01_alt_scoring" / "stripped_poses"

# Subset: WT control + the charge-reversal mutants + a few non-charge sanity cases.
POSES = [
    ("WT_holo",          "wt_holo_top.pdbqt",          "baseline (reference)"),
    ("R215E_holo",       "R215E_holo_top.pdbqt",       "Arg+1 → Glu−1 charge reversal (clamp residue)"),
    ("R175E_holo",       "R175E_holo_top.pdbqt",       "Arg+1 → Glu−1 single (clamp residue, partner chain)"),
    ("R176E_holo",       "R176E_holo_top.pdbqt",       "Arg+1 → Glu−1 single (clamp residue)"),
    ("R175E_R176E_holo", "R175E_R176E_holo_top.pdbqt", "DOUBLE charge reversal (max electrostatic penalty)"),
    ("R215A_holo",       "R215A_holo_top.pdbqt",       "Arg+1 → Ala neutral (loss of clamp, no flip)"),
    ("R215A_N226A_holo", "R215A_N226A_holo_top.pdbqt", "Arg→Ala + Asn→Ala (clamp + orient)"),
    ("C195A_holo",       "C195A_holo_top.pdbqt",       "catalytic Cys→Ala (the v3 illusion case)"),
    ("C195A_H196A_holo", "C195A_H196A_holo_top.pdbqt", "double catalytic ablation"),
    ("H196A_holo",       "H196A_holo_top.pdbqt",       "catalytic His→Ala"),
]

# Phase 14 cavity-18 + Plevitrexed extras
EXTRAS = [
    ("plevitrexed_apo", REPO / "14_inhibitor_design" / "02_cofactor_site" / "docked" / "Plevitrexed_apo_seed42.pdbqt",
        REPO / "14_inhibitor_design" / "02_cofactor_site" / "receptor_holo.pdbqt", "S2 above-noise hit"),
    ("indazole_cav18",  REPO / "14_inhibitor_design" / "04_allosteric" / "docked" / "cav18_frag_CID7032_seed42.pdbqt",
        APO_RECEPTOR, "S4 cavity-18 top hit"),
    ("ibuprofen_cav18", REPO / "14_inhibitor_design" / "04_allosteric" / "docked" / "cav18_frag_CID3672_seed42.pdbqt",
        APO_RECEPTOR, "S4 cavity-18 #2"),
]

# Custom scoring: extend Vina defaults with explicit electrostatic + desolvation
# Smina format: one term per line, "<weight> <term_spec>"
# Default Vina weights (built-in) restated, then added Coulomb (1/r repulsion-style)
# and AD4 solvation (Huey 2007). The electrostatic weight 0.30 is chosen to
# make the R→E flip cross the ±0.85 noise floor without distorting hydrophobic terms.
CUSTOM_SCORING = """\
-0.035579    gauss(o=0,_w=0.5,_c=8)
-0.005156    gauss(o=3,_w=2,_c=8)
0.840245     repulsion(o=0,_c=8)
-0.035069    hydrophobic(g=0.5,_b=1.5,_c=8)
-0.587439    non_dir_h_bond(g=-0.7,_b=0,_c=8)
0.300000     electrostatic(i=2,_^=100,_c=8)
0.100000     ad4_solvation(d-sigma=3.6,_s/q=0.01097,_c=8)
1.923000     num_tors_div
"""

# Electrostatic-amplified variant — weight 3.0 (10× the default-ish weight)
# This is deliberately aggressive to see whether the R→E sign error inverts
# at any reasonable electrostatic weighting. If even weight=3.0 fails, the
# issue is positional (rigid pose can't move into the new electrostatic field),
# not the scoring function's weight.
CUSTOM_SCORING_AMPLIFIED = """\
-0.035579    gauss(o=0,_w=0.5,_c=8)
-0.005156    gauss(o=3,_w=2,_c=8)
0.840245     repulsion(o=0,_c=8)
-0.035069    hydrophobic(g=0.5,_b=1.5,_c=8)
-0.587439    non_dir_h_bond(g=-0.7,_b=0,_c=8)
3.000000     electrostatic(i=2,_^=100,_c=8)
0.500000     ad4_solvation(d-sigma=3.6,_s/q=0.01097,_c=8)
1.923000     num_tors_div
"""

def write_custom_scoring(path: Path):
    path.write_text(CUSTOM_SCORING)
    return path

def write_amplified_scoring(path: Path):
    path.write_text(CUSTOM_SCORING_AMPLIFIED)
    return path

def extract_model1(pose_pdbqt: Path) -> Path:
    """Multi-model Vina outputs need MODEL 1 extracted for Smina --score_only."""
    if not pose_pdbqt.exists(): return pose_pdbqt
    # if file already single-model, return as-is
    text = pose_pdbqt.read_text()
    if "MODEL 2" not in text and "MODEL " not in text:
        return pose_pdbqt
    out = pose_pdbqt.with_suffix(".m1.pdbqt")
    if out.exists(): return out
    keep = []; in_m1 = False
    for ln in text.splitlines():
        if ln.startswith("MODEL 1"): in_m1 = True; continue
        if ln.startswith("ENDMDL") and in_m1: break
        if in_m1: keep.append(ln)
    out.write_text("\n".join(keep) + "\n")
    return out


def smina_score(receptor: Path, ligand: Path, scoring: str | None = None,
                custom: Path | None = None, minimize: bool = False) -> dict:
    lig_use = extract_model1(ligand)
    if minimize:
        # Minimize first, then score the relaxed pose
        out_pdbqt = lig_use.with_suffix(".min.pdbqt")
        cmd = [SMINA, "--minimize", "-r", str(receptor), "-l", str(lig_use),
               "--out", str(out_pdbqt)]
        if scoring: cmd += ["--scoring", scoring]
        if custom:  cmd += ["--custom_scoring", str(custom)]
        try:
            proc = subprocess.run(cmd, check=False, capture_output=True, timeout=180)
        except subprocess.TimeoutExpired:
            return {"ok": False, "err": "minimize_timeout"}
        if proc.returncode != 0 or not out_pdbqt.exists():
            return {"ok": False, "err": proc.stderr.decode()[:300]}
        out = proc.stdout.decode()
    else:
        cmd = [SMINA, "--score_only", "-r", str(receptor), "-l", str(lig_use)]
        if scoring: cmd += ["--scoring", scoring]
        if custom:  cmd += ["--custom_scoring", str(custom)]
        try:
            proc = subprocess.run(cmd, check=False, capture_output=True, timeout=120)
        except subprocess.TimeoutExpired:
            return {"ok": False, "err": "timeout"}
        if proc.returncode != 0:
            return {"ok": False, "err": proc.stderr.decode()[:300]}
        out = proc.stdout.decode()
    # parse Affinity:
    m = re.search(r"Affinity:\s*([-\d.]+)", out)
    intra = re.search(r"Intramolecular energy:\s*([-\d.]+)", out)
    return {"ok": True,
            "score": float(m.group(1)) if m else None,
            "intra": float(intra.group(1)) if intra else None}

def main():
    print("=== Phase 14e — Smina rescoring with electrostatics ===")
    custom_path = write_custom_scoring(OUT / "custom_scoring_q.txt")
    amp_path    = write_amplified_scoring(OUT / "custom_scoring_qamp.txt")
    print(f"  custom scoring → {custom_path}")
    print(f"  amplified electrostatic scoring → {amp_path}")

    rows = []
    # ---- 1. Phase 7-8 mutant rescoring ----
    print("\n  Rescoring Phase 7-8 holo mutant top poses (vs WT_holo baseline):")
    hdr = f"  {'pose':<22} {'vina':>7} {'vinardo':>8} {'custom_q':>9} {'q_amp':>7} {'min_q':>7}  comment"
    print(hdr); print("  " + "-"*(len(hdr)-2))
    for label, fname, comment in POSES:
        lig = STRIPPED / label / fname
        if not lig.exists():
            print(f"  ! missing: {lig}"); continue
        scores = {}
        scorers = [("vina",      "vina",    None,        False),
                   ("vinardo",   "vinardo", None,        False),
                   ("custom_q",  None,      custom_path, False),
                   ("q_amp",     None,      amp_path,    False),
                   ("min_q",     None,      custom_path, True)]   # minimize+custom
        for sname, sarg, scust, mflag in scorers:
            r = smina_score(APO_RECEPTOR, lig, scoring=sarg, custom=scust, minimize=mflag)
            scores[sname] = r.get("score")
        rows.append({"pose": label, "category": "phase7_mutant",
                     "ligand": str(lig.relative_to(REPO)),
                     **scores, "comment": comment})
        def fmt(v, w):
            return f"{v:>{w}.2f}" if v is not None else f"{'—':>{w}}"
        print(f"  {label:<22} {fmt(scores['vina'],7)} {fmt(scores['vinardo'],8)} "
              f"{fmt(scores['custom_q'],9)} {fmt(scores['q_amp'],7)} {fmt(scores['min_q'],7)}  {comment[:42]}")

    # ---- 2. Phase 14 cavity 18 + S2 Plevitrexed ----
    print("\n  Rescoring Phase 14 extras:")
    for label, lig, receptor, comment in EXTRAS:
        if not lig.exists() or not receptor.exists():
            print(f"  ! missing pose or receptor for {label}: {lig} | {receptor}"); continue
        scores = {}
        for sname, sarg, scust, mflag in [("vina", "vina", None, False),
                                          ("vinardo", "vinardo", None, False),
                                          ("custom_q", None, custom_path, False),
                                          ("q_amp", None, amp_path, False),
                                          ("min_q", None, custom_path, True)]:
            r = smina_score(receptor, lig, scoring=sarg, custom=scust, minimize=mflag)
            scores[sname] = r.get("score")
        rows.append({"pose": label, "category": "phase14_extra",
                     "ligand": str(lig.relative_to(REPO)),
                     **scores, "comment": comment})
        def fmt(v, w):
            return f"{v:>{w}.2f}" if v is not None else f"{'—':>{w}}"
        print(f"  {label:<22} {fmt(scores['vina'],7)} {fmt(scores['vinardo'],8)} "
              f"{fmt(scores['custom_q'],9)} {fmt(scores['q_amp'],7)} {fmt(scores['min_q'],7)}  {comment[:42]}")

    # ---- write CSV ----
    csv_path = OUT / "rescore_results.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"\n  → {csv_path}")

    # ---- Δ vs WT_holo + summary CSV ----
    wt = next(r for r in rows if r["pose"] == "WT_holo")
    summary = []
    cols = ("vina", "vinardo", "custom_q", "q_amp", "min_q")
    for r in rows:
        if r["category"] != "phase7_mutant": continue
        s = {"pose": r["pose"], "comment": r["comment"]}
        for c in cols:
            s[c] = r.get(c)
            s[f"delta_{c}_vs_WT"] = ((r[c] - wt[c]) if (r.get(c) is not None and wt.get(c) is not None) else None)
        summary.append(s)
    summary_path = OUT / "rescore_summary.csv"
    with summary_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader(); w.writerows(summary)
    print(f"  → {summary_path}")

    # ---- Δ table to stdout ----
    print("\n  Δ vs WT_holo (positive = penalised, negative = improved):")
    print(f"  {'pose':<22} {'ΔVina':>7} {'ΔVinardo':>8} {'Δcustom_q':>10} {'Δq_amp':>8} {'Δmin_q':>8}  ★ if |Δ|>0.85")
    print("  " + "-"*100)
    for s in summary:
        if s["pose"] == "WT_holo": continue
        vals = [s[f"delta_{c}_vs_WT"] for c in cols]
        def fmt(v, w):
            return f"{v:+{w}.2f}" if v is not None else f"{'—':>{w}}"
        flags = [c for c, v in zip(cols, vals) if v is not None and abs(v) >= 0.85]
        flag_str = " ★ " + ",".join(flags) if flags else ""
        print(f"  {s['pose']:<22} {fmt(vals[0],7)} {fmt(vals[1],8)} {fmt(vals[2],10)} "
              f"{fmt(vals[3],8)} {fmt(vals[4],8)}{flag_str}")

    # ---- plot ----
    mutants = [s for s in summary if s["pose"] != "WT_holo"]
    labels = [s["pose"].replace("_holo","") for s in mutants]
    x = np.arange(len(labels)); width = 0.18
    fig, ax = plt.subplots(figsize=(13, 5.5))
    series = [
        ("Vina (no electrostatics)",                                "#bdc3c7", "vina"),
        ("Vinardo (stiffer hydrophobic)",                           "#7f8c8d", "vinardo"),
        ("Smina custom (Vina + electrostatic + desolvation)",       "#e67e22", "custom_q"),
        ("Smina q-amplified (3× electrostatic)",                    "#c0392b", "q_amp"),
        ("Smina min+custom (minimize then score)",                  "#2c5f2d", "min_q"),
    ]
    for i, (label, color, key) in enumerate(series):
        ys = [s[f"delta_{key}_vs_WT"] or 0 for s in mutants]
        ax.bar(x + (i - 2) * width, ys, width, label=label, color=color)
    ax.axhline(y=0, color="black", linewidth=0.7)
    ax.axhline(y= 0.85, ls=":", color="#888", alpha=0.6, label="Vina noise floor ±0.85")
    ax.axhline(y=-0.85, ls=":", color="#888", alpha=0.6)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Δ score vs WT_holo  (positive = worse binding)")
    ax.set_title("Smina rescoring + minimization vs rigid Vina — does adding electrostatics fix the R→E sign error?")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout(); plt.savefig(OUT / "rescore_plot.png", dpi=140); plt.close()
    print(f"  → {OUT / 'rescore_plot.png'}")

if __name__ == "__main__":
    main()
