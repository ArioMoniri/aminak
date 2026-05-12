# v2 Code Review

Agent role: Code reviewer (round 2 — verify v1 fixes plus look for new bugs).
Date: 2026-05-12.

## v1-gap closure summary

| # | v1 issue | Status | Evidence |
|---|----------|--------|----------|
| 1 | JSD weighted window | PASS | stage1:183-198 |
| 2 | Exclude >50% gap from percentile | PASS | stage1:315-318 |
| 3 | Real TYMS orthologs + bifunctional trim | PASS | stage1:23-35, 104-113 |
| 4 | Keep chains A & B; CME43→CYS | PASS | stage3:75, 49-68 |
| 5 | Vina parse from PDBQT | PASS | stage6:75-85 |
| 6 | Identical box/seed/exh WT vs mutants | PASS | stage6:166 = stage7:207 |
| 7 | G217W drop on clash <1.8 Å | PASS | stage7:276-279 |
| 8 | Apo + holo conditions | PASS | stage6:206, stage7:306 |
| 9 | Centralised config | FAIL | P04818/1HVY hardcoded across 4+ files |
| 10 | Box 18³ + exh ≥32 | PASS | stage6:166-167 |

**Verdict: 9/10 v1 gaps closed.** Only #9 (centralised config) is unaddressed — pipeline is still TYMS-specific.

## Per-script details

### `stage1_msa_v2.py`
- JSD weighted window correctly implemented per Capra & Singh.
- Eligible mask `gap_fraction <= 0.5` filters before percentile ranking.
- Ortholog list hardcoded (P04818, P0A884, P00469, P06785 etc.); the comment notes P07807 was wrong (DHFR) and substitutes P06785. Bifunctional trimming via `TS_MOTIF` works on Plasmodium/Arabidopsis.

### `stage3_structure_v2.py`
- `accept_chain` returns `chain.id in ("A","B")` ✓
- CME43→CYS preservation: pre-pass rewrites `CME` HETATM lines into `ATOM ... CYS` records, drops 2-hydroxyethyl atoms, keeps only canonical CYS heavy atoms.

### `stage6_dock_wt_v2.py`
- REMARK VINA RESULT parsing from PDBQT ✓ (no MODEL 1 vs 10 ordering bug).
- Box 18³ + exhaustiveness 32 + seed 42 ✓
- Apo + holo loop ✓

### `stage7_mutants_v2.py`
- G217W drop on heavy-atom clash <1.8 Å ✓
- Both apo and holo per mutant ✓
- New issue: `mean_top3 = nan` rows are EXPECTED — Vina returned <3 modes for tight holo pockets (verified R50A_holo PDBQT has only 2 REMARK lines).

### `stage9_report_v2.py`
- **Operator-precedence bug at line 209**: `f"...n>0)" if cor is not None else f""` evaluates the conditional ONLY on the trailing `f"...n>0)"` fragment. When `cor is None`, the catalytic-summary text disappears entirely.

### `build_mutation_plot.py`
- **TYPE_STYLE keys mismatch real categories**: keys are `single_ala`, `single_opposite`, `arg_clamp`, `double_catalytic`, `double`, `control_surface` — but stage 7 produces `ala_scan`, `opposite`, `arg_clamp`, `double_dyad`, `double_phosclamp`, `double_polar_neutral`, `double_substrate_orient`, `double_aromatic_swap`, `control_surface`, `explore_g217w`. The substring matcher (`if k in typ`) returns the default grey for almost every point. Real bug — legend collapses.
- Docstring says `results_full.csv` but reads `mutant_results_v2.csv`.

### `build_viewers.py`
- Backtick / backslash escaping is correct (escape `\` first then `` ` ``).
- `${` is not escaped, but PDB text never contains it — safe in practice.
- Active-site residue list (line 87) is hardcoded.

## New v2 issues

- `build_mutation_plot.py` `style_for` substring-match keys don't match stage 7's real category names → legend collapses to default grey for almost every point.
- `stage9_report_v2.py:209` operator precedence bug in `v2_vs_v1_summary` f-string when `cor is None`.
- `mean_top3=nan` for holo runs is **not** a parser bug — Vina returned <3 poses (verified PDBQT for R50A_holo has 2 REMARK lines).
- `build_mutation_plot.py:12` docstring/path drift.
- `split_top_pose_pdbqt` literal `MODEL 1` match remains brittle for MODEL 10+.

## Cross-cutting findings
- Receptor prep: meeko `--no-flexible` flag is rejected (`error: unrecognized arguments`); silent fallback to `obabel -xr` is used for **every** receptor across stages 6 and 7. Documented in pipeline.log.
