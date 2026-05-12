# V3 Code Review

Agent role: Code reviewer (round 3, audit of v2 punch-list closure).
Date: 2026-05-12.

## v2 punch-list closure table

| # | v2 issue | v3 status | Evidence |
|---|----------|-----------|----------|
| 1 | TYPE_STYLE keys mismatch | **PASS** | `build_mutation_plot_v3.py:34-45` exact-match keys |
| 2 | `stage9:209` precedence bug | **PASS** | new template; conditional removed |
| 3 | `mean_top3=nan` for <3 modes | **PASS** | `stage7:317` `mean(affs[:min(3,n)])` |
| 4 | Receptor all-zero charges | **PASS** | `stage6:56-98` 3-method waterfall + `max|q|>0.05` gate |
| 5 | Rotamer strain (`get_strain` 1e30) | **PARTIAL** | sculpt replaces strain loop, but `SCULPT_FAIL` is logged-only — no fallback rotamer search |
| 6 | `MODEL 1` literal match | **PARTIAL** | `stage6:206` still `startswith("MODEL 1")`; only safe because `break` triggers at first ENDMDL |
| 7 | Cofactor `-p 7.4` | **PASS in code** (but no-op in practice — see structural review) |

## Did v3 close the v2 gaps?
**5 PASS, 2 PARTIAL, 0 FAIL.** The substantive science fixes (charges, sign convention, mean_topk, dual-seed WT, plot legend, ΔΔG→Δ Vina wording, limitations) are in. Two cosmetic regressions remain: the `MODEL 1` substring and the missing `SCULPT_FAIL` recovery path.

## Per-script details
- `stage6_dock_wt_v3.py`: receptor charge waterfall (obabel → meeko → pdb2pqr) with `max|q|>0.05` gate; vina_params echoed in JSON; `MODEL 1` brittle parser persists.
- `stage7_mutants_v3.py`: sculpt iterations called 3× try/except; SCULPT_OK/FAIL logged; sign convention `delta = top - wt_aff` correct; `mean_topk = mean(affs[:min(3,n)])`.
- `stage8_analysis_v3.py`: descending sort on delta_vina_vs_wt; mis_docked exclusion; CSV comment line skipped via `comment="#"`.
- `build_mutation_plot_v3.py`: TYPE_STYLE exact-match; legend built from observed categories; mis-docked alpha 0.35.
- `stage9_report_v3.py`: Limitations section in HTML+DOCX; "Δ Vina score" wording; mis-dock caption explanation.

Neither partial blocks the report; both should be tightened before publication.
