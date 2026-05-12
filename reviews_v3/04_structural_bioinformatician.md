# V3 Structural Bioinformatics Audit

Agent role: Structural bioinformatician (round 3).
Date: 2026-05-12.

## VERDICT: CONDITIONAL PASS

Three v2 blockers fixed; **one fix is fake (cofactor pH 7.4) and one fix is incomplete (atom-name preservation)**.

## Per audit point

**1. Receptor charges — PASS.**
- `protein_dimer_apo.pdbqt`: n=5657, mean|q|=0.144, max|q|=0.507
- `protein_dimer_holo.pdbqt`: n=5782, mean|q|=0.142, max|q|=0.507
- All 1070 zero-charge atoms in holo are AutoDock `HD` polar hydrogens (correctly folded into parent in AD4 united-atom convention) plus 2 stray `C` atoms.

**2. Rotamer minimisation (sculpt) — PASS.**
- `pipeline.log` shows `sculpt: OK` for all 20 mutants.
- `R175E_mut.pdb` spot-check: 0 heavy-atom clashes <2.0 Å around GLU175. Carboxylate geometry correct (CD–OE1 1.18 Å, CD–OE2 1.25 Å, OE1–CD–OE2 = 120.1°). Sculpt actually ran and produced sane geometry.

**3. Atom-name preservation / pose-RMSD reproducibility — FAIL.**
- `wt_holo_top_pose.pdb` atom names are still generic Open Babel labels (`C`, `N`, `O`, `P` with no primes), not the UMP PDB standard names.
- Atom counts also differ (22 v3 vs 20 v2 reference).
- An RMSD-by-atom-name reproducer is impossible.

**4. Cofactor protonation at pH 7.4 — FAIL (no-op fix).**
- `md5(cofactor_chainA_h_v3.pdb) == md5(cofactor_chainA_h.pdb v2)` — byte-identical.
- `obabel -p 7.4` was called but did not change the file: input PDB lacked bond-order info, so obabel preserved existing wrong H placement. **Carboxylates still protonated as neutral COOH.**
- At pH 7.4 a glutamate γ-carboxylate (pKa ≈ 4.3) should be deprotonated COO⁻.
- Cofactor electrostatics in `protein_dimer_holo.pdbqt` still wrong → all holo Vina scores carry the same bias as v2.

**5. mean_topk NaN bug — PASS.**
- 0 NaN/empty in 40 mutant rows; rows with `n_modes=2,3,4,19,20` all produced finite means.

**6. Box / seed / exhaustiveness echoed — PASS.**
- Both `wt_*.json` record `exhaustiveness: 96`, `num_modes: 32`, `box_size: 22`, `seed_primary: 42`, `seed_sanity: 7`, centroid coords and source.

## New issues found

**A. Best-seed selection is circular (medium severity).**
`stage6_dock_wt_v3.py:307-313` picks `best_seed = min(... key=rmsd_to_native)`. Choosing the seed whose top pose is closest to the crystal pose, then reporting "RMSD-to-crystal = 2.08 Å" as a sanity metric, is self-fulfilling.

**B. WT holo n_modes = 3 is genuinely low (low severity).**
Even at exhaustiveness 96 and box 22³, only 3 unique poses survive Vina's clustering. mean_topk over 3 modes is noisy.

## Required corrections before PASS

1. **Re-protonate cofactor properly**: feed obabel an SDF/MOL with bond orders (e.g. fetch SMILES for D16 from PDB-CCD). All holo dockings must be re-run.
2. **Preserve UMP atom names** through the docking pipeline: dock the residue-name PDB directly without obabel renaming, or maintain an atom-index-to-name mapping JSON.

## Lower-priority improvements
- Replace circular RMSD-based seed selection with affinity-based or consensus selection.
- Flag `n_modes < 5` rows in the analysis CSV.
- Document that 1070 zero-charge atoms are AD4 polar hydrogens (avoid future false alarms).
- Add a CI assertion that `cofactor_chainA_h_v3.pdb` is **not** identical to v2 after `reprotonate_cofactor` runs.
