# V4 Code Review

Round 4. Date: 2026-05-12.

## v3 punch-list closure

| # | Item | Status | Evidence |
|---|------|--------|----------|
| 1 | Real cofactor reprotonation + MD5 assert | PASS | `stage3_cofactor_v4.py` L60-99, L277-279 |
| 2 | WT ≥5-seed sweep, affinity selection | PASS | `stage6_dock_wt_v4.py` L126, L168-170 |
| 3 | Atom-name preservation (index map) | PASS (greedy mapping caveat) | `_common_v4.py` L158-213 |
| 4 | `low_confidence` column + ranking exclusion | PASS | `stage7` L158, `stage8` L86-87 |
| 5 | Spearman correlation | PASS | `stage8` L131, L138 |
| 6 | Noise floor / polyglutamylation / AD4 polar H | PASS | `stage9` L164-193 |
| 7 | WT row in mutant CSV | PASS | `stage7` L241-264 |

## New v4 issues (non-blocking)

1. **stage3 numpy dtype** — `np.array([conf.GetAtomPosition(...)])` may produce `dtype=object`; could silently corrupt Kabsch. Should coerce to `float`.
2. **stage6 vina_params** — when fallback exh=128 fires, `vina_params_base` still records exh=96. Reproducibility metadata understates exhaustiveness.
3. **stage6/7 protocol asymmetry** — WT apo uses 5-seed/exh=96, mutant apo reuses v3 single-seed/exh=32. Δ_apo mixes two protocols.
4. **`restore_atom_names` greedy mapping** — no symmetry handling; can mis-name equivalent atoms (phosphate Os, ring symmetry). Visualisation-grade only.

All scripts pass closure; no code regressions. The protocol-asymmetry note belongs in the methods section, not as a code fix.
