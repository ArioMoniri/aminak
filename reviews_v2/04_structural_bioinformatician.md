# v2 Structural Bioinformatics Audit

Agent role: Structural bioinformatician — round 2 (deep methods audit).
Date: 2026-05-12.

## Verdict: **CONDITIONAL PASS**

The two blocking v1 issues — broken MSA and CME-contaminated structure — are genuinely fixed and the docking pipeline now runs to completion across the planned mutant set with consistent box/seed/exhaustiveness. However:
1. **The receptor PDBQT contains all-zero partial charges** (silent meeko-CLI fallback);
2. **The rotamer-strain selection silently degrades to "always rotamer 0"**;
3. The JSON-reported pose RMSD does not reproduce by atom-name alignment.

These three issues invalidate the affinity numbers as quantitative ΔG estimates, so rigor is "pass for geometry, conditional on rerun for energetics."

## Per-audit-point findings

1. **MSA quality** — PASS. 10 orthologs at 399 columns, reference gap fraction **0.216** (vs v1 0.677). The `LPxMALPPCHAL` block aligns cleanly across all 10 species.
2. **PfDHFR-TS trimming** — PASS. P13922 trimmed to 315 aa via `trim_bifunctional_ts()`. Arabidopsis Q05762 also trimmed.
3. **Dimer prep** — PASS. `protein_dimer_h.pdb` has chain A (4626 ATOM) + chain B (4626 ATOM). Position 43 chain A is `CYS A 43` with full standard atoms; **0 occurrences of `CME`** anywhere.
4. **Box parameters** — PASS. Vina log confirms grid 18×18×18 Å, exhaustiveness 32, seed 42.
5. **Mutant heavy-atom clashes (<1.8 Å)** — PASS for the three spec mutants (R175E, Y258F_F225Y, N226D show 0 clashes).
6. **G217W dropped** — PASS. No `_apo*` / `_holo*` / `_top*` files.
7. **Pose RMSD** — FAIL for atom-name reproducibility. WT apo top pose has only 4 unique atom names (obabel stripped of `'`), crystal has 20 distinct names (`N1, C2, ... C5', O5', P, OP1, OP2, OP3`). Index-aligned heavy-atom RMSD ≈ 5.37 Å vs JSON 5.60 Å — close but not identical.
8. **Ligand multi-format** — PASS. `dump.{pdb,mol2,sdf,pdbqt}` all exist and non-empty. PDBQT has 5 active torsions, 22 atoms, Gasteiger charges.
9. **Receptor charges** — **FAIL**. **All 5,661 protein ATOM records in `protein_dimer_apo.pdbqt` carry `+0.000`.** Root cause: `mk_prepare_receptor.py --no-flexible` is rejected, meeko exits rc=2, silent fallback to `obabel -xr` writes no charges. Vina therefore evaluates electrostatics as **zero on the receptor side**. Neither AM1-BCC nor Gasteiger was applied to the protein.

## New issues found

- **Rotamer-strain selection is a no-op.** Pipeline.log shows `ROTAMER_PICK ... best_idx=0 strain=1e+30` for **every single mutation in every mutant**. PyMOL `get_strain()` returns `1e30` (or raises) in headless mode, so the loop never updates `best_idx`. Every mutant uses rotamer 0 (the default), with no relaxation.
- **`mean_top3` NaN** is a counting bug. Flagged `*_holo` rows generated only 1–3 docking modes (e.g., R50A holo MODELs=2, Q214A holo MODELs=1, N226D holo MODELs=3), too few for a top-3 mean → NaN propagates.
- **Cofactor (raltitrexed/D16) protonation** uses `obabel -h` only; **no `-p 7.4` flag**, so the carboxylate tail is neutralised rather than carrying its physiological −1 charge. With zero receptor charges this is doubly problematic.

## Specific corrections needed

1. **Receptor charges** — fix `scripts/v2/stage6_dock_wt_v2.py:38`: drop `--no-flexible` (replace with `--default_altloc A` or omit) so meeko runs; OR call `mk_prepare_receptor.py -i {prot_in} -o {prot_out_basename}` (newer meeko writes both `.pdbqt` and `.json`); OR call `obabel <prot> -O <prot.pdbqt> -p 7.4 --partialcharge gasteiger -xr`. Re-run all stage 6/7 dockings.
2. **Rotamer pick** — call `cmd.local_minimize` after `apply()`, OR explicit chi-angle scanning (`cmd.set_dihedral`) + `cmd.get_model("byres ...").energy`.
3. **Cofactor protonation** — change to `obabel src -O dst -h -p 7.4`.
4. **`mean_top3`** — change to `mean(top_k)` with `k = min(3, len(affs))`; report `n_modes` alongside.
5. **Pose-RMSD reproducibility** — preserve native atom names through obabel (use `meeko export_docking_poses.py`).

## Lower-priority improvements
- Echo box/exhaustiveness/seed into `wt_*.json` so audit doesn't need the log.
- E. coli ortholog gap fraction (0.338) is high; consider gap-trimming alignment before scoring.
- D218N_N226D may have been silently dropped (only 19 mutants in CSV vs 21 attempted; G217W is one).
- Add sanity-check assertion in stage 6 that fails if `max(|charge|) < 0.01`.
