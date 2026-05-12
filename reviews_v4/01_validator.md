# V4 Validator Review

Round 4. Date: 2026-05-12.

## ✅ PASSED — 10 / 10
1. File integrity — report.docx 198 KB, .pdf 252 KB, .html 245 KB; all v4 dirs populated.
2. Cofactor MD5 differs from v2/v3 (`e6f06a08…` vs `c3a555f6…`).
3. WT holo n_modes = 30 (≥ 10).
4. Affinity-based seed selection (best seed 256, lowest top affinity −5.242).
5. UMP atom-name preservation in `wt_holo_top_pose.pdb`.
6. CSV integrity: 42 data rows incl. WT row in both conditions; `low_confidence` column present.
7. T170A control Δ ≈ 0 in both conditions (apo +0.180, holo +0.171).
8. Limitations strings present: `Vina noise floor` ×5; `polyglutamylat` ×2.
9. AD4 polar-H zero-charge note present.
10. Spearman correlation in `summary_v4.json` (ρ = −0.4507, p = 0.191).

## ⚠️ FLAGS
- WT holo RMSD vs crystal = 12.95 Å (the reported "ionised cofactor expels dUMP" interpretation; flagged transparently).
- Top destabilising holo Δ = +0.204 — well below Vina's 0.85 noise floor.
- 21 of 42 rows mis_docked — half the holo dockings lost the pocket.
- Filtered Spearman ρ = +0.60 has n = 4; should not be cited.

## 🔧 RECOMMENDATIONS
1. Lead with the noise-floor caveat in the abstract.
2. Consider co-folding (AlphaFold3 or Boltz-1) for the holo state.
3. Report apo as primary, holo as exploratory.
4. Add n=4 disclaimer next to filtered Spearman.
5. No code changes needed for v4.
