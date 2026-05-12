# Phase-6 Code Review

Round 6. Date: 2026-05-12.

## Verdict: CONDITIONAL PASS

Functional and largely correct. Minor improvements recommended.

## Per-script

### `_common.py`
Clean, minimal. Hard-coded Homebrew binary paths — Linux/Intel Mac users will break. Suggest `shutil.which(...)` fallback.

### `step1_clean_pdb.py`
Identity check uses Biopython `PairwiseAligner` with substitution_matrix=None — crude. Magic constants `TYMS_P04818`, `1hvy`, hard-coded path.

### `step2_blast.py`
BLAST `hitlist_size=250, expect=10.0, word_size=3`. `EXCLUDE_PDB={"1HVY"}` correctly removes source. Identity 30-95%, coverage ≥0.80, resolution ≤2.5 Å filters enforced. **Concern**: resolution check is permissive — `if res is not None and res > MAX_RES: skip` means unparsable REMARK 2 is silently accepted. Recommend `if res is None or res > MAX_RES: skip`.

### `step3_clustalw_pir.py`
ClustalW with default scoring. PIR header `sequence:target:26:A:313:A:::-1.00:-1.00` — verified 10 colon-separated fields. **The 11-field bug WAS fixed.** PASS.

### `step4_run_modeller.py`
`AutoModel(env, alnfile, knowns, sequence, assess_methods=(assess.DOPE, assess.GA341))` — correct Modeller 10.x API. Score parsing reads from `a.outputs` directly. **No `random_seed` set** — Modeller seeds refinement with system entropy by default. Recommend `env.rand_seed = -12345`.

### `step5_pymol_compare.py`
`cmd.align("mdl and name CA", "crystal and name CA")` — sequence-based alignment on Cα. RMSDs are Cα-only, after outlier rejection. Not documented in CSV header. Recommend renaming to `ca_rmsd_refined_A`.

### `step6_validate.py`
Ramachandran regions are **hand-drawn polygons** (lines 29-41), not derived from a Lovell reference distribution. Glycine/proline are scatter-coloured but classified against general polygons — technically wrong (each has its own Ramachandran). DOPE per-residue uses `Selection(model).assess_dope` — proper API.

### `step7_viewers_report.py`
Inlines full PDB text into HTML — fine for ~20 KB structures. Backtick-escaping defensible.

## Cross-cutting findings
- **Reproducibility**: no `env.rand_seed` anywhere — Modeller refinement is non-deterministic.
- **Magic strings**: `1HVY`, chain `A`, residue 26/313, hard-coded paths.
- **Silent fallbacks**: step 2 unparsable resolution accepted; step 6 missing per-residue DOPE only logs a warning.

## Required fixes before publication-grade sign-off
1. Set `env.rand_seed`.
2. Treat unparsable resolution as fail in step 2.
3. Document RMSD provenance (Cα refined-align).
4. Move hand-drawn Ramachandran caveat into a visible label on the plot.
