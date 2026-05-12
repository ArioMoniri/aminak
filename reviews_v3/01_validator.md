# V3 Validator Review

Agent role: Validator (round 3, against v3 outputs).
Date: 2026-05-12.

## ✅ PASSED

1. **File integrity** — All v3 directories exist with expected contents:
   - `09c_report_v3/` (3 files: docx 314 KB, html 414 KB, pdf 377 KB) — all > 100 KB ✓
   - `08c_analysis_v3/` (12 files including 5 PNGs, ranked CSVs, summary)
   - `07c_mut_docking_v3/viewer_files/` (80 PDB files)
   - `06c_docking_wt_v3/` (28 files)
   - `scripts/v3/` (5 stage scripts)
   - PIL spot-check: `mutation_effect_plot.png` 2203×977 RGBA, `delta_vina_apo_holo.png` 1820×650, `apo_vs_holo_concordance.png` 1172×939 — all valid

2. **Receptor charges fixed** — Both PDBQTs have non-zero Gasteiger charges:
   - APO: max |q| = **0.507**, 4650/5657 atoms non-zero
   - HOLO: max |q| = **0.507**, 4712/5782 atoms non-zero

3. **WT holo dock fixed** — `wt_holo.json`:
   - top_affinity = **−8.314** (vs v2's −3.215, much more negative ✓)
   - n_modes = **3** ✓ (≥ 3)
   - rmsd_top_to_native = **2.08 Å** ✓ (< 2.5)

4. **Sign convention correct**:
   - T170A apo Δ = +0.110, holo Δ = +0.029 ✓ (control near 0)
   - H196A holo Δ = **+0.503** ✓ (positive = destabilising)

5. **Mis-dock annotation present**: `mis_docked` column exists; True=2, False=38.

6. **mean_topk fix verified**: 9 rows have n_modes=2; mean_topk over 2 affinities, NOT NaN. Formula `mean(top min(3, n_modes))` confirmed.

7. **G217W absent** from v3 mutant list ✓

8. **Viewers regenerated**: `viewers/index.html` exists; **86 HTML files**; `wt_holo_complex.html` is 754 KB and contains **9,274 embedded ATOM lines** (v3 receptor is embedded).

9. **CSV header documents sign convention**: First line has the convention header.

10. **Report files** all present with sizes > 100 KB.

## ❌ FAILED
None — all 10 verification criteria satisfied.

## ⚠️ FLAGS
- WT row not in mutant CSV (delta computed against external WT).
- mis_docked threshold spot-check not exhaustive.
- n_modes=2 frequency in holo (9 of 20 holo rows return only 2 poses) — likely pocket geometry, but ensemble statistics will have higher variance.
- wt_apo top affinity (−9.127) is more negative than wt_holo (−8.314) — physically plausible (cofactor occupies part of pocket) but worth a one-line note.

## 🔧 RECOMMENDATIONS (severity-ordered)
1. Low: Add a WT reference row (delta=0) to `mutant_results_v3.csv`.
2. Low: Verify the 2 `mis_docked=True` rows have RMSD > 3.0.
3. Low: Add a 1-sentence note explaining apo > holo affinity.
4. Info: Document n_modes < 5 issue.

**Verdict: v3 outputs PASS all 10 criteria.**
