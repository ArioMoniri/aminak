#!/usr/bin/env python3
"""Phase 8c: integrate 8a + 8b results into a master comparison plot and
print the C195A / R215E / holo-Δ-magnitude diagnostic interpretation."""
from __future__ import annotations

import csv
import os
import statistics
import sys
from pathlib import Path

PROJECT = Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
ROOT = PROJECT / "13_phase8"
ALT = ROOT / "01_alt_scoring" / "alt_scoring_results.csv"
FLEX = ROOT / "02_flexres" / "flexres_compare.csv"


def load_alt() -> list[dict]:
    rows = []
    if not ALT.exists():
        return rows
    with ALT.open() as fh:
        for r in csv.DictReader(fh):
            for k in ("vina_score", "vinardo_score", "ad4_score"):
                try:
                    r[k] = float(r[k]) if r.get(k) else None
                except ValueError:
                    r[k] = None
            rows.append(r)
    return rows


def load_flex() -> dict[str, dict]:
    res = {}
    if not FLEX.exists():
        return res
    with FLEX.open() as fh:
        for r in csv.DictReader(fh):
            for k in ("rigid_vina_score", "flex_vina_score", "delta_flex"):
                try:
                    r[k] = float(r[k]) if r.get(k) else None
                except ValueError:
                    r[k] = None
            res[r["label"]] = r
    return res


def master_html(alt_rows: list[dict], flex_data: dict[str, dict]) -> None:
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        print("plotly not available; skip master HTML")
        return

    # Build a per-label table: holo rows from alt + flex map
    by_label = {}
    for r in alt_rows:
        if r["condition"] != "holo":
            continue
        by_label[r["label"]] = {
            "vina": r["vina_score"],
            "vinardo": r["vinardo_score"],
            "ad4": r["ad4_score"],
        }
    for lbl, fr in flex_data.items():
        by_label.setdefault(lbl, {}).update({"flex_vina": fr["flex_vina_score"]})

    # 4-panel: bar chart per scoring column
    labels = list(by_label.keys())

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=("Rigid Vina", "Vinardo", "AD4", "Flex Vina (priority 8)"),
        horizontal_spacing=0.10, vertical_spacing=0.18,
    )

    # Sort by Vina score
    sorted_labels = sorted(labels, key=lambda lb: (by_label[lb].get("vina") or 0))

    def bar(col: str, row: int, col_n: int, color: str):
        ys = [by_label[lb].get(col) for lb in sorted_labels]
        valid = [(lb, y) for lb, y in zip(sorted_labels, ys) if y is not None]
        if not valid:
            return
        fig.add_trace(
            go.Bar(
                x=[v[0] for v in valid],
                y=[v[1] for v in valid],
                marker=dict(color=color),
                showlegend=False,
                text=[f"{v[1]:.2f}" for v in valid],
                textposition="outside",
                textfont=dict(size=9),
            ),
            row=row, col=col_n,
        )
        fig.update_yaxes(title_text="kcal/mol", row=row, col=col_n)

    bar("vina", 1, 1, "#1f77b4")
    bar("vinardo", 1, 2, "#ff7f0e")
    bar("ad4", 2, 1, "#bdbdbd")  # likely all empty
    bar("flex_vina", 2, 2, "#2ca02c")

    fig.update_layout(
        title="Phase 8 master comparison: holo scores across scoring schemes",
        width=1280, height=900, template="plotly_white",
    )
    fig.update_xaxes(tickangle=-40, tickfont=dict(size=9))
    out = ROOT / "master_comparison.html"
    fig.write_html(str(out), include_plotlyjs="cdn")
    print(f"wrote {out}")


def diagnose(alt_rows: list[dict], flex_data: dict[str, dict]) -> str:
    """Compose the C195A / R215E / holo-Δ-magnitude interpretation."""
    # Index by (label, condition)
    by_key = {(r["label"], r["condition"]): r for r in alt_rows}

    def get(label, cond, key):
        r = by_key.get((label, cond))
        return r.get(key) if r else None

    wt_vina_holo = get("WT", "holo", "vina_score")
    wt_vinardo_holo = get("WT", "holo", "vinardo_score")
    wt_ad4_holo = get("WT", "holo", "ad4_score")

    def delta_vs_wt(label, cond, key, wt):
        v = get(label, cond, key)
        if v is None or wt is None:
            return None
        return v - wt

    # C195A: vina says it's a *better* binder (delta_vina_vs_wt = -2.246 in v5);
    # does Vinardo/AD4/flex agree it's better?
    c195a_vina_h = get("C195A", "holo", "vina_score")
    c195a_vinardo_h = get("C195A", "holo", "vinardo_score")
    c195a_flex = flex_data.get("C195A", {}).get("flex_vina_score")
    c195a_rigid = flex_data.get("C195A", {}).get("rigid_vina_score")

    # R215E: charge reversal at conserved Arg; should be destabilising
    r215e_vina_h = get("R215E", "holo", "vina_score")
    r215e_vinardo_h = get("R215E", "holo", "vinardo_score")
    r215e_flex = flex_data.get("R215E", {}).get("flex_vina_score")

    # Holo Δ magnitudes across schemes (avg absolute Δ vs WT)
    holo_rows = [r for r in alt_rows if r["condition"] == "holo" and r["label"] != "WT"]
    def avg_abs_delta(key, wt):
        if wt is None:
            return None
        vals = [abs(r[key] - wt) for r in holo_rows if r.get(key) is not None]
        if not vals: return None
        return statistics.mean(vals)

    avg_abs_vina_h = avg_abs_delta("vina_score", wt_vina_holo)
    avg_abs_vinardo_h = avg_abs_delta("vinardo_score", wt_vinardo_holo)

    lines = []
    lines.append("=== Phase 8 diagnostic interpretation ===")
    lines.append("")
    lines.append(f"WT holo scores: vina={wt_vina_holo}, vinardo={wt_vinardo_holo}, ad4={wt_ad4_holo}")
    lines.append("")
    lines.append("1) C195A illusion check (was 'too good' under rigid Vina):")
    if c195a_vina_h is not None and wt_vina_holo is not None:
        lines.append(f"   rigid Vina (holo):    C195A {c195a_vina_h:+.2f}  vs WT {wt_vina_holo:+.2f}  Δ={c195a_vina_h - wt_vina_holo:+.2f}")
    if c195a_vinardo_h is not None and wt_vinardo_holo is not None:
        lines.append(f"   Vinardo (holo):       C195A {c195a_vinardo_h:+.2f}  vs WT {wt_vinardo_holo:+.2f}  Δ={c195a_vinardo_h - wt_vinardo_holo:+.2f}")
    if c195a_flex is not None and c195a_rigid is not None:
        lines.append(f"   flex Vina (priority): C195A flex={c195a_flex:+.2f}  vs C195A rigid={c195a_rigid:+.2f}  Δ={c195a_flex - c195a_rigid:+.2f}")
    # Did the illusion get fixed? Check sign of Δ vs WT
    delta_v = (c195a_vina_h - wt_vina_holo) if (c195a_vina_h is not None and wt_vina_holo is not None) else None
    delta_vd = (c195a_vinardo_h - wt_vinardo_holo) if (c195a_vinardo_h is not None and wt_vinardo_holo is not None) else None
    fixed = []
    if delta_v is not None and delta_v < 0:
        fixed.append(f"rigid Vina still shows C195A as 'better' (Δ={delta_v:+.2f})")
    elif delta_v is not None:
        fixed.append(f"rigid Vina has C195A worse than WT (Δ={delta_v:+.2f})")
    if delta_vd is not None and delta_vd < 0:
        fixed.append(f"Vinardo also says C195A is 'better' (Δ={delta_vd:+.2f}) — illusion NOT fixed")
    elif delta_vd is not None:
        fixed.append(f"Vinardo says C195A is worse than WT (Δ={delta_vd:+.2f}) — illusion FIXED by Vinardo")
    lines.append("   verdict: " + "; ".join(fixed))
    lines.append("")
    lines.append("2) R215E charge-reversal penalty check:")
    if r215e_vina_h is not None and wt_vina_holo is not None:
        lines.append(f"   rigid Vina: R215E {r215e_vina_h:+.2f}  vs WT {wt_vina_holo:+.2f}  Δ={r215e_vina_h - wt_vina_holo:+.2f}")
    if r215e_vinardo_h is not None and wt_vinardo_holo is not None:
        lines.append(f"   Vinardo:    R215E {r215e_vinardo_h:+.2f}  vs WT {wt_vinardo_holo:+.2f}  Δ={r215e_vinardo_h - wt_vinardo_holo:+.2f}")
    if r215e_flex is not None:
        lines.append(f"   flex Vina:  R215E flex {r215e_flex:+.2f}")
    lines.append("")
    lines.append("3) Holo Δ magnitudes (mean |Δ vs WT| across non-WT holo rows):")
    if avg_abs_vina_h is not None:
        lines.append(f"   rigid Vina mean |Δ|: {avg_abs_vina_h:.3f} kcal/mol")
    if avg_abs_vinardo_h is not None:
        lines.append(f"   Vinardo   mean |Δ|: {avg_abs_vinardo_h:.3f} kcal/mol")
    if avg_abs_vina_h is not None and avg_abs_vinardo_h is not None:
        ratio = avg_abs_vinardo_h / avg_abs_vina_h
        more_or_less = "LARGER (more discrimination)" if ratio > 1.0 else "SMALLER (less discrimination)"
        lines.append(f"   Vinardo Δ ratio vs Vina = {ratio:.2f}x ({more_or_less})")

    return "\n".join(lines)


def write_readme(diagnostic: str) -> None:
    out = ROOT / "README.md"
    text = f"""# Phase 8: Beyond rigid Vina

This phase asks: *can we do better than the rigid-receptor Vina baseline?*
Three sub-tasks, all using AutoDock Vina 1.2.7 (the only Apple-Silicon-native
docking binary available; GNINA does not ship for arm64-darwin).

## 8a -- Alternative scoring functions (Vinardo + AD4)

Directory: `01_alt_scoring/`

Re-score the existing Vina top poses with Vina's other built-in scoring
functions:

- **Vinardo** (Quiroga & Villarreal 2016) is a refined empirical scoring
  function with stiffer hydrophobic and repulsive terms; it is known to
  rank actives vs decoys better than default Vina on several benchmarks.
- **AD4** is the AutoDock 4 force-field-style scoring.  In Vina 1.2 AD4
  scoring **requires precomputed affinity maps** (via `autogrid4`).  The
  `autogrid4` binary is not available on this Apple Silicon host (no
  Homebrew formula, no arm64 build from Scripps).  AD4 column is therefore
  empty in `alt_scoring_results.csv`; this is a documented limitation, not
  a script bug.

Driver: `scripts/v8/task_8a_alt_scoring.py`
Plots: `scripts/v8/task_8a_plots.py`

The mutant **apo** receptor PDBQTs were destructively overwritten by Vina
docking output in the v3 pipeline (see `scripts/v3/stage7_mutants_v3.py` --
`out_pdbqt` overwrites `rec_apo`).  8a regenerates them on-the-fly from
`{{mut}}_mut_h.pdb` using `obabel -xr -p 7.4 --partialcharge gasteiger`
(same recipe as v3's `prepare_receptor_with_charges` first method).

Outputs:
- `alt_scoring_results.csv`  -- 42 rows (WT_apo + WT_holo + 20 mut x apo,holo)
- `alt_scoring_compare.png`  -- static Vinardo-vs-Vina scatter
- `alt_scoring_compare.html` -- interactive Plotly scatter

## 8b -- Flexible-residue Vina re-dock (8 priority mutants)

Directory: `02_flexres/`

For each of the 8 priority mutants (same list as `scripts/v7/task_a_replicas.py`),
re-dock with the 14-residue active-site panel made flexible via Vina's
built-in `--flex` flag.

The standard tooling (`prepare_flexreceptor4.py` from MGLTools; meeko's
`mk_prepare_receptor.py -f`) was unusable on this dataset:
- MGLTools is not installed (no arm64 build).
- Meeko fails with `RuntimeError: Updated 1 H positions but deleted 7`
  on every PDB we tried (apo, holo, with and without explicit H, after
  prody clean-up, after pdb2pqr clean-up).

`scripts/v8/flex_split.py` is a self-contained replacement that emits the
exact Vina flex PDBQT format (`BEGIN_RES`/`ROOT`/`BRANCH`/`ENDBRANCH`/
`END_RES`).  It uses hardcoded chi-rotation topology templates for each
amino acid (ALA through VAL) and walks the templates to emit one BRANCH
per chi torsion.

Receptor build pipeline per mutant:
1. `obabel {{mut}}_mut_h.pdb -xr -p 7.4 --partialcharge gasteiger`
   -> apo PDBQT (chains A/B preserved; residues renumbered 1..N).
2. Shift residue numbers by +25 to restore source numbering
   (source PDBs start at residue 26 in chain A).
3. Concatenate cofactor PDBQT atoms (`06f_receptor_fixed/cofactor_A.pdbqt`
   and `cofactor_B.pdbqt`) to make holo PDBQT.
4. `flex_split.split_clean_pdbqt(holo_text, FLEX_PANEL)`
   -> rigid PDBQT (everything except the 14 panel side chains)
      + flex PDBQT (BEGIN_RES blocks for each panel side chain).
5. `vina --receptor rigid --flex flex --ligand dump.pdbqt --exhaustiveness 32
   --num_modes 20 --seed 42 --center/size {{...same box as v5...}}`

Outputs:
- `{{mutant}}_rigid.pdbqt`     -- rigid portion of the holo receptor
- `{{mutant}}_flexres.pdbqt`   -- the 14 panel side chains in flex format
- `{{mutant}}_flex.pdbqt`      -- Vina docking output (up to 20 modes)
- `{{mutant}}_flex.log`        -- Vina stdout/stderr
- `flexres_compare.csv`        -- summary: label, rigid_vina_score, flex_vina_score, delta_flex
- `flex_vs_rigid.png`/`.html`  -- scatter of flex vs rigid scores

## 8c -- Documentation + master comparison

`master_comparison.html` -- four-panel bar chart of holo scores across
rigid Vina, Vinardo, AD4 (empty), and flex Vina.

The diagnostic interpretation (C195A illusion, R215E charge reversal,
holo-Δ magnitudes) is printed to stdout at the end of
`scripts/v8/task_8c_integrate.py`.

## Why these three are "better than rigid Vina"

| Aspect                            | Rigid Vina | Vinardo | AD4 | Flex Vina |
|-----------------------------------|:---------:|:-------:|:---:|:---------:|
| Different energy function         |     -     |   yes   | yes |    -      |
| Better hydrophobic terms          |     -     |   yes   |  -  |    -      |
| Receptor side-chain rearrangement |     -     |    -    |  -  |   yes     |
| Same pose space                   |    yes    |   yes   | yes |    -      |
| Apple-Silicon-native              |    yes    |   yes   | no  |   yes     |

Each addresses a different failure mode of rigid Vina.  Vinardo asks
"would a better scoring function change the ranking on the *same* poses?"
Flex Vina asks "would induced fit of side chains allow a different pose
that the rigid receptor blocked?"  AD4 would have given a force-field
view, but is unavailable here.

## Limitations

- No induced-fit *backbone* movement (only side-chain rotamers in 8b).
- No proper continuum electrostatics (Vinardo and Vina both use a
  distance-dependent dielectric; R215E charge effects are still
  approximated, not Poisson-Boltzmann).
- The 14-residue flex panel adds ~30-40 rotatable degrees of freedom,
  pushing the conformational search hard; `--exhaustiveness 32` may be
  too low for full convergence.  A more thorough study would scan
  exhaustiveness or use multi-replica seeds (Phase 7's recipe).
- AD4 column omitted (no autogrid4 binary on this Apple Silicon host).
- The custom `flex_split.py` topology templates cover the 20 standard
  amino acids but assume canonical atom names (PDB v3.3 conventions);
  non-canonical residues or alternate names (e.g. legacy 1HB vs HB1) are
  best-effort.

## Diagnostic interpretation (auto-generated)

{diagnostic}
"""
    out.write_text(text)
    print(f"wrote {out}")


def main() -> int:
    alt_rows = load_alt()
    flex_data = load_flex()
    if not alt_rows:
        print("ERROR: alt_scoring_results.csv missing; run task_8a_alt_scoring.py first")
        return 1
    master_html(alt_rows, flex_data)
    diag = diagnose(alt_rows, flex_data)
    write_readme(diag)
    print()
    print(diag)
    return 0


if __name__ == "__main__":
    sys.exit(main())
