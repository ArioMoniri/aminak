# Code Review

Agent role: Code reviewer (correctness, robustness, reproducibility of the Python scripts).
Date: 2026-05-12.

## stage1_msa.py
- **JSD vs Capra & Singh 2007 — partially correct**: Uses Robinson-Robinson background frequencies and the `(1-λ)·p + λ·BG` pseudocount with λ=0.5 (✓). Gap penalty `(1 - gap_frac)` is applied (✓). Reference-row gap exclusion is correct.
- **Window weighting NOT per spec**: Capra & Singh use a weighted window `0.5·s[i] + 0.25·(s[i-1]+s[i+1])`; `windowed()` does an unweighted mean of three positions.
- **BLOSUM62-derived background claim**: Comment says "BLOSUM62-derived background" but uses Robinson & Robinson 1991 frequencies — mislabeled.
- `subprocess.run(["cp", aligned_fa, final_fa])` — use `shutil.copyfile`.

## stage2_active_site.py
- Hardcoded `REF="P04818"`, `PDB="1hvy"`, `target_ligands={"UMP","D16"}`, `catalytic_must_include=[195,196]`. Not parameterised.
- Assumes `1HVY` author residue numbering == UniProt numbering. Acknowledged but not programmatically checked.
- No retry/backoff on PDBe HTTP calls.

## stage3_structure.py
- Hardcoded `"UMP"`, `"D16"`, chain `"A"`, accession `"1HVY"`.
- Selected-residue verification logs `"WARNING"` but does not exit — downstream stages will silently dock against a structure missing requested residues.
- `obabel` rc nonzero is logged but not fatal.

## stage4_pymol.py
- `int(line[22:26])` ignores PDB insertion codes (col 27).
- `score*100` packed into B-factor field — fine for visualization, should clamp to fit `%6.2f`.

## stage5_6_dock_wt.py
- **`parse_vina_stdout` is fragile**: regex matches a 4-column row `mode | affinity | rmsd_l.b. | rmsd_u.b.`. Vina 1.2.x sometimes emits 3 columns. Should parse `REMARK VINA RESULT` lines from the output PDBQT instead.
- **Top-pose split is brittle**: `line.startswith("MODEL 1")` matches `MODEL 10..19` too. Need `line.rstrip() == "MODEL 1"` or regex `^MODEL\s+1\s*$`.
- `rmsd_to_native` greedy nearest-neighbour fallback is not symmetric and not optimal-assignment; use Hungarian (`scipy.optimize.linear_sum_assignment`).
- `rmsd_to_native` does NOT superpose; it computes raw coordinate RMSD. Correct for "redocked vs crystal in same frame" but undocumented.

## stage7_mutants.py
- **Same box / seed / exhaustiveness as WT**: ✓ correct.
- **RMSD comparison**: mutant top pose vs WT crystal dUMP — ✓ correct.
- **Receptor re-prep**: per mutant — ✓.
- `add_h_obabel` return value ignored — silent failure → empty `mut_h` → confusing downstream error.
- `pymol_mutate` does not check that the selection actually matched a residue.
- `opposite_for` dict has dead entries (Y→A but Y branch comment says Y→D).
- Hardcoded double-mutant residue ids assume residues are in `selected`; if not, `res_letter.get(p, "A")` falls back to `"A"` and silently mislabels.
- Top-pose split has same `MODEL 1` prefix bug as stage 5_6.

## stage8_analysis.py
- `mutation_id` regex `([A-Z])(\d+)([A-Z])` will match `CTRL_T170A` weirdly. Filtering by `type` first prevents the bug.
- Conservation lookup defaults missing positions to 0 — silently misclassifies missing-data points as low-conservation.

## stage9_report.py
- `md_to_html` is a hand-rolled markdown subset; emits `<li>` without enclosing `<ul>` (invalid HTML).
- `n_uniprot=int((asdf["source"]!="pdbe").sum())`, `n_pdbe=int((asdf["source"]!="uniprot").sum())` — these double-count "both". Counts will sum to more than `n_db`.
- Bare `except Exception: pass` around control-distance recomputation hides real errors.

## Cross-cutting findings
- **Error swallowing**: stage9 (control distance), stage5_6 (`pdbqt_models` parsing per-line `except Exception: pass`).
- **Path safety**: all subprocess calls use list form (no `shell=True`). User input is not interpolated.
- **Reproducibility seeds**: Vina seed=42 ✓ everywhere. MAFFT `--auto` is non-deterministic across versions — consider pinning algorithm.
- **Hardcoded protein-specific constants** appear in: stage1, stage2, stage3, stage4, stage7, stage9. To swap proteins, edits required in ≥6 files.

## Punch list (impact-ordered)
1. **Fix Vina affinity parsing** — switch to `REMARK VINA RESULT` lines from the output PDBQT.
2. **Fix `MODEL 1` prefix bug** — regex-anchor model split.
3. **Stage7 hardcoded double-mutant residue numbers** must be validated against `selected` and `res_letter` before use.
4. **Capra & Singh windowing** — implement weighted window or update docstring.
5. **Stage9 UniProt/PDBe counts** double-count "both" — fix the source-count expressions.
6. **Optimal-assignment RMSD** instead of greedy nearest-neighbour.
7. **Stage3 missing-residue check should be fatal**, not warning.
8. **Centralize protein-specific config** into one JSON read by every stage.
9. **Replace `subprocess.run(["cp", …])`** with `shutil.copyfile`.
10. **Stage9** wrap `<li>` blocks in `<ul>`; replace bare `except: pass`.
11. **Stage7** check `add_h_obabel` and `pymol_mutate` return codes.
