# Phase-6 Validator Review

Round 6. Date: 2026-05-12.

## ✅ PASSED — 13 / 13

1. File integrity: all subdirs + scripts present.
2. Step 1: chain A only, 0 HETATMs.
3. Step 2 templates: 3IHI 92.71%, 6K7Q 75.96%, 5H39 72.28% — all 30-95%; 1HVY excluded; templates/_A.pdb present.
4. Step 3 PIR: exactly 10 colon-separated fields per header. ✓
5. Step 4 models: 10 + best_model.pdb; ATOM count >2000 each.
6. Step 4 scores CSV: 10 rows + header, no NaN.
7. Step 5 RMSD: 0.367–0.390 Å.
8. Step 5 figures: pairwise ×10 + all_models_overlay; PIL-decode OK.
9. Step 6 Ramachandran: favoured 83.5–85.3% (just over 80% floor).
10. Step 6 SAVES doc: explicit upload instructions.
11. Viewers: modeller_model{1..10}.html present; index.html updated.
12. DOCX: report_PHASE6.docx 6.7 MB.
13. <100% educational rule: max identity 92.71%.

## ⚠️ FLAGS

- Folder naming: `07_viewers/` (not `07_summary` per spec wording — functionally OK).
- Duplicate viewers: `10_modeller/07_viewers/` and `viewers/modeller_model*.html` (mirrored).
- Ramachandran %favoured borderline low (83-85% vs typical 90% for Modeller homology models).
- Modeller intermediate outputs (`target.D*`, `target.V*`) consume disk, not documented.
- Stray non-chain template PDBs (`1RTS.pdb`, `3IHI.pdb` etc.) alongside `_A.pdb` versions.

## 🔧 RECOMMENDATIONS
- Run UCLA SAVES on `best_model.pdb` for canonical PROCHECK Ramachandran.
- Delete duplicate `10_modeller/07_viewers/` or have index.html reference it.
- Document `target.D*` / `target.V*` intermediates.
- Clean stray template PDBs.
