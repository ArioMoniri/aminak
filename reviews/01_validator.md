# Validator Review (read-only verification)

Agent role: Validator (file integrity, number consistency, reproducibility).
Date: 2026-05-12.

## ✅ PASSED

- **Report files exist**: `09_report/report.html` (3.7 MB) and `09_report/report.pdf` (2.65 MB) on disk. All images in the HTML are inlined as `data:image/png;base64,...` (no external `src` references), so no broken links possible.
- **Stage-4 PyMOL screenshots present and decodable** (1600×1200 RGBA): `01_overview.png` (562 KB), `02_closeup.png` (668 KB), `03_conservation.png` (517 KB), `04_cavity.png` (61 KB — see FLAGS).
- **All 21 mutant PNGs in `07_mut_docking/`** are real images, all ≥351 KB (smallest D218K.png).
- **CSV row count correct**: `07_mut_docking/results_full.csv` has 23 lines = 1 header + 22 data rows (1 WT + 8 single-Ala + 7 single-opposite + 5 doubles + 1 control). Matches log line `Wrote ... results_full.csv with 22 rows` and the report's `Total mutants tested: 21` (excludes WT).
- **WT top affinity reproducible**: re-parsing `06_docking_wt/wt_poses.pdbqt` → 20 `REMARK VINA RESULT` lines, mode 1 = **-7.729 kcal/mol**. Matches CSV (-7.729), `wt_result.json`, `vina_wt.log` mode 1 row, and report's "-7.73 kcal/mol".
- **D218K Δ verified**: CSV `delta_vs_wt = +1.12` matches report ("Δ=+1.12") and log (`Δ=+1.12 rmsd=7.34`). RMSD 7.339 Å.
- **Y258A RMSD verified**: CSV `rmsd_top_to_native = 4.5130` matches report ("RMSD 4.51 Å") and log (`rmsd=4.51`).
- **Ortholog count = 8** in `01_msa/aligned.fa` (Homo, Mus, E. coli, L. casei, S. cerevisiae, C. elegans, T4 phage, P. falciparum) — matches report and ≥6 spec.
- **Stage-end markers present** for stages 1, 2, 3, 4, 5, 6, 7, 8, 9 in `pipeline.log`.
- **Timestamps monotonic forward** across the entire `pipeline.log` (0 violations).
- **All 22 screenshot paths in CSV** point to existing PNGs.

## ❌ FAILED

- **Cys195 conservation peak is wrong / misleading**: the report claims Cys195 is part of the catalytic core but `01_msa/conservation_scores.csv` line for position 195 shows `C, JS=0.0861, percentile=36.7` — i.e. it is in the **bottom-third** of conserved positions, not a peak. The augmentation rule should be stated as "Cys195 is included from the literature *despite* its sub-25% JSD".
- **MSA reference-row gap fraction is extreme**: P04818 in `aligned.fa` is **657 / 970 = 67.7 % gaps** (alignment length 970, reference ungapped length 313). Caused by one input being the bifunctional `P21520` (P. falciparum), which the report mentions but downstream JSD/numbering steps appear to ignore.

## ⚠️ FLAGS

- `04_pymol/04_cavity.png` is only **60.96 KB**, an order of magnitude smaller than the other three Stage-4 PNGs (517–668 KB). Still a valid 1600×1200 RGBA image; PyMOL `cavity_mode` rendering is sparse.
- **STAGE2 logged "starting" three times** and "DONE" only twice — Stage 2 was clearly re-run after a failure. Not flagged in the report.
- **Every single mutant docking logged `meeko failed: mk_prepare_receptor.py: error: unrecognized arguments: --no-flexible`** (22 occurrences). The pipeline transparently fell back to `obabel -xr`; the failures are silent in the report's analysis.
- **Identical-row coincidence**: `R175E_R176E` and `CTRL_T170A` produce **identical** `top_affinity (-7.694)`, `delta (+0.035)`, `mean_top3 (-7.528)`, and `rmsd (1.05098)` to 5+ decimals. Statistically improbable for two different mutants — worth investigating.
- **`conservation_scores.csv` has rows 1–313 only** (the ungapped reference numbering), so the JSD column numbering matches the human protein, not the 970-column MSA. Internally consistent but worth a one-line docstring.

## 🔧 RECOMMENDATIONS (severity-ordered)

1. **High** — Add a per-mutant prep-status row in the report's Stage 5/6 section so the silent meeko failures are visible.
2. **High** — Investigate why `R175E_R176E` and `CTRL_T170A` produce 5-decimal-identical Vina results.
3. **Medium** — Replace any wording that implies "Cys195 conservation peak"; the empirical percentile is 36.7.
4. **Medium** — Re-render `04_cavity.png` with a tighter clip plane.
5. **Low** — Emit explicit `STAGE2: ABORT/SUCCESS` markers for grep-level integrity checks.
6. **Low** — Document in the CSV header that `js_score` numbering is the ungapped P04818 reference.
