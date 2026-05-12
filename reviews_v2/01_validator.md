# v2 Validator Review

Agent role: Validator (file integrity, number reproduction, viewer reality).
Date: 2026-05-12 (round 2, against v2 outputs).

## ✅ PASSED
- All v2 folders populated; 80 viewer PDBs in `07b_mut_docking_v2/viewer_files/`; 85 viewer HTMLs in `viewers/` matching the 85 anchors in `viewers/index.html`.
- 3 PNGs decode OK with PIL; `workflow_diagram_v2.png` is 439 KB (well over the 100 KB floor).
- `mutant_results_v2.csv` = 41 lines (header + 40 = 20 mutants × 2 conditions). `ranked_apo.csv` and `ranked_holo.csv` = 21 lines each. Zero NaN in `top_affinity` across all three CSVs.
- WT holo top affinity reproducible: `wt_holo.pdbqt` first `REMARK VINA RESULT = -3.215`, matches `wt_holo.json`, `summary_v2.json.wt_holo`, and report HTML.
- MSA: 10 sequences in `aligned.fa`; reference (P04818) gap fraction **21.6 %** (v1 was 67.7 %). The conserved `LPxMALPP[CR]` motif is present in human/mouse/rat full block, and in the core for casei/yeast/fly/arabidopsis/T4/plasmodium.
- Conservation: Cys195 percentile = **100.0** (v1 was 36.7 %). H196 96.1, R175 89.5, R176 90.8, R215 94.4, N226 97.1.
- G217W: NO `G217W*` files in `07b_mut_docking_v2/viewer_files/`; absent from all CSVs; `summary_v2.json.skipped = ["G217W", "g217w_clash:8"]`.

## ❌ FAILED
None.

## ⚠️ FLAGS
1. **G217W stub directory still present** at `07b_mut_docking_v2/G217W/` (PDB + PML), even though all downstream artefacts honour the drop. Either delete or rename to `G217W_skipped/`.
2. **`n_panel` vs panel size discrepancy**: `summary_v2.json` reports `n_panel=21` while only 20 mutants have docking rows.
3. **`mean_top3 = nan` cluster is holo-only** (10/40 rows). Plausible (low pose count in tight holo pocket) but worth confirming Vina is asked for ≥3 poses uniformly.
4. **Conservation column resolution**: CSV mixes `ref_position` (UniProt numbering) and `aln_col` (alignment column) — ensure consumers always key on `ref_position`.

## 🔧 RECOMMENDATIONS (severity-ordered)
1. (low) Remove or relabel `07b_mut_docking_v2/G217W/`.
2. (low) Document the `mean_top3=nan` rule inline in the CSV header or `summary_v2.json`.
3. (low) Show WT holo affinity to 3 decimals (-3.215) in the report.
4. (info) Add MD5 checksum manifest for viewer files so future audits can verify integrity in O(1).
