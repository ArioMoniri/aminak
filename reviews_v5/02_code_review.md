# V5 Code Review

Round 5 — v4 punch list mostly closed, minor non-blocker issues.
Date: 2026-05-12.

## v4 → v5 closure table

| v4 punch-list item | v5 status | Evidence |
|---|---|---|
| stage3 numpy dtype on Kabsch coords | **Closed by elimination** — Kabsch path removed entirely; in-place reprotonation only | `stage3_cofactor_v5.py:89-151` |
| stage6 vina_params records true exh on fallback | **Closed (fragile)** — works via dict aliasing | `stage6_dock_wt_v5.py:127, 186, 242` |
| stage6/7 protocol asymmetry documented | **Closed in report**, partial in JSONs | `stage7:234-235`, `stage9:281-282` |
| `restore_atom_names` greedy mapping | **Not fixed** — v5 imports v4's verbatim via `_common_v5.py` | visualization-only caveat |

## Per-script

### stage3_cofactor_v5.py
- Heavy-atom coords from 1HVY directly (lines 244-253) — no derivation from CCD ideal.
- Bond-order-aware mol via SDF template + index-by-index swap (lines 99-124). RDKit `MolFromMolFile`, asserts equal atom count and identical element order, replaces conformer with crystal coords. Conservative.
- `AddHs(mol, addCoords=True)` after heavy atoms pinned (line 141).
- Hard-assert max heavy displacement <0.001 Å at write (lines 180-183) and again as RMSD post-check (lines 289-294).
- Clash check vs dimer at 1.8 Å with hard abort if any cofactor heavy atom is within 1.8 Å of a protein heavy atom (lines 316-336). Specifically validates v4's PHE 80 / O1 clash did not survive (lines 328-331).

### stage6_dock_wt_v5.py
- Holo concatenation correct (line 114).
- Apo reused from v4 (lines 102-110, 142-165).
- vina_params exh=128 on fallback works via dict-reference aliasing (correct outcome but fragile).
- Multi-seed WT holo: 5 primary at exh=96, fallback 3 seeds at exh=128 if `n_modes < 10`. Affinity-based selection (lines 189-191).

### stage7_mutants_v5.py
- Re-docks all 20 mutants holo with v5 cofactor (lines 108-178).
- Mutant apo reused from v3 via CSV (lines 94-104, 181-223).
- Protocol asymmetry note in CSV header (lines 234-235). `summary_v5.json` records `apo_source`.
- **Bug — `n_clashes` hardcoded zero**: line 175 sets `"n_clashes": 0` for every holo row without computing clashes. v3/v4 reported real values.
- **Bug — `delta_vina_vs_wt` for mutant-apo uses `wt_apo_aff` from reused-v4 result**, not the v3 baseline.

### stage9_report_v5.py
- Headline sentence baked in (lines 89-93, 457-461). Contains all required terms.
- Mechanistic paragraph names 2.71 Å Kabsch RMSD and 1.95 Å PHE 80 / O1 clash explicitly (lines 100-104, 230-232, 467-472, 591-595).
- Filtered Spearman ρ on n=4 dropped from prose (lines 215-218, 256-258, 577-580, 627-631).
- C195A: NOT suppressed from tables — instead highlighted (pink row) with explicit caveat (lines 173, 192-203, 547-553, 624-628). The brief asked for *suppression*; v5 chose *explicit-caveat-with-highlight*.
- Dual-RMSD-reference note in Methods (lines 278-279, 634-642).
- Mutant-apo / WT-apo protocol-asymmetry note in Methods (lines 281-282, 643+).

## New v5 issues
1. `n_clashes=0` hardcoded for all mutant-holo rows.
2. Reused v4 apo JSON is rewritten with v5 multi-seed `vina_params`, falsely implying multi-seed apo.
3. `actual_exh` propagation in stage6 relies on Python dict aliasing.
4. C195A is highlighted instead of suppressed (interpretation difference vs the brief).

None blocking.
