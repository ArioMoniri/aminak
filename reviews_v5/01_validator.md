# V5 Validator Review

Round 5 — pipeline APPROVED.
Date: 2026-05-12.

## ✅ PASSED — 10 / 10
1. Files: report.docx 202 KB, .pdf 257 KB, .html 252 KB; all v5 dirs populated; 80 viewer PDBs.
2. Cofactor in-place fix verified: heavy-atom RMSD 0.000000 Å vs 1HVY chain-A D16 (32/32 atoms exact); MD5 differs from v4.
3. Carboxylates deprotonated: SMILES contains 2× `[O-]`; provenance JSON confirms.
4. Zero protein clashes: 0 contacts <1.8 Å (also 0 <2.0 Å) for both chains.
5. WT holo top affinity −8.249 kcal/mol; RMSD 0.334 Å vs crystal (<3 Å). All 3 selected seeds converge to within 0.036 kcal/mol.
6. Viewers: index.html 8.9 KB; wt_holo_complex.html 753 KB.
7. Headline finding present in report HTML ("noise floor" ×6, "principal finding" hit).
8. C195A pink-flagged in HTML; suppressed from `ranked_holo_clean.csv`.
9. Top destabilisers match doer claim (R215A_N226A +0.772, H196A +0.491, R215E +0.470 holo).
10. No mutant exceeds Vina ±0.85 noise floor in clean tables (max |Δ| = 0.772).

## ❌ FAILED
None.

## ⚠️ FLAGS
- R50A holo (Δ +0.828) and R50E holo (Δ +0.482) are mis_docked; correctly excluded from clean ranking but remain in raw `ranked_holo.csv`.
- R215A_N226A apo Δ = +1.15 DOES exceed 0.85 noise floor — the "no mutant exceeds noise" claim holds for HOLO ONLY; apo signal may legitimately be flagged.
- WT holo n_modes is only 2 even at exh=128 — funnel collapse explanation needs hedging.

## 🔧 RECOMMENDATIONS
- Add a sanity assertion that apo `ranked_clean` is bounds-checked against noise floor.
- Lock the v5 provenance JSON into the docx/pdf appendix.
- Boost exhaustiveness/num_modes for v6 holo (current top RMSD 0.33 Å is excellent; non-blocking).

**Overall: v5 cleanly resolves the v4 cofactor-placement blocker. Pipeline approved.**
