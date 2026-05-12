# Stage 8 Analysis – Human Thymidylate Synthase (P04818) / dUMP

_Generated 2026-05-12 11:56_

## Wild-type baseline

- Re-docking of dUMP into 1HVY chain A reproduces the crystal pose with RMSD **1.08 Å** (heavy-atom).
- Top affinity: **-7.73 kcal/mol**, mean of top-3 poses -7.56 kcal/mol.
- The pose recapitulation gives confidence in the box geometry (22³ Å, centred on the catalytic-residue centroid) and parameter set.

## Mutational panel summary

- Total mutants tested: **21** (8 Ala-scan, 7 chemically opposite singles, 5 doubles, 1 distant-surface control).

- Δ affinity range: **-0.38** to **+1.12** kcal/mol vs WT.

- Mutants with pose displacement (RMSD > 3 Å): D218K, N226A, Y258A, Y258F_F225Y.


## Most disruptive mutations (largest positive Δ)

- **D218K** (single_opposite): Δ=+1.12 kcal/mol, top-pose RMSD 7.34 Å.
- **Y258A** (single_ala): Δ=+0.84 kcal/mol, top-pose RMSD 4.51 Å.
- **N226A** (single_ala): Δ=+0.73 kcal/mol, top-pose RMSD 5.88 Å.

## Most affinity-enhancing mutations (largest negative Δ)

- **C195S** (single_opposite): Δ=-0.38 kcal/mol, top-pose RMSD 0.55 Å.
- **C195A** (single_ala): Δ=-0.27 kcal/mol, top-pose RMSD 0.89 Å.
- **C195S_H196N** (double_polar_neutral): Δ=-0.25 kcal/mol, top-pose RMSD 0.65 Å.

## Catalytic-residue mutants

- **C195A**: Δ=-0.27 kcal/mol, RMSD 0.89 Å. Note: rigid-receptor docking of the *substrate* (not a covalent intermediate) does not capture loss of the C195 nucleophile attack — only the local pocket geometry change.
- **H196A**: Δ=+0.12 kcal/mol, RMSD 1.08 Å. Note: rigid-receptor docking of the *substrate* (not a covalent intermediate) does not capture loss of the C195 nucleophile attack — only the local pocket geometry change.
- **C195A_H196A**: Δ=-0.10 kcal/mol, RMSD 0.87 Å. Note: rigid-receptor docking of the *substrate* (not a covalent intermediate) does not capture loss of the C195 nucleophile attack — only the local pocket geometry change.

## Distant-surface control
- **CTRL_T170A**: Δ=+0.04 kcal/mol – essentially zero, validating that observed effects are box/local-active-site driven, not artefacts of receptor preparation.


## Conservation vs effect

- Pearson correlation between per-residue JSD conservation and |Δ affinity| across single mutants: **r = 0.03** (n=15). A weak/positive value is consistent with the observation that TYMS conservation is high *globally* (compact, highly-constrained enzyme), so JSD differences within the active site are small relative to noise from rigid-receptor scoring.


## Caveats

- Vina is rigid-receptor; PyMOL mutagenesis wizard provides only a rotamer pick (no minimization). Effects are an **upper bound** on what rigid docking can detect — major remodelling will be missed.

- Substrate-level docking does not model the covalent Michaelis intermediate or the methylene-tetrahydrofolate cofactor; biological loss of activity for C195 mutants will exceed the modest Δ affinity reported here.

- Doubles whose constituent singles independently disturb the pocket can compensate (e.g., C195S_H196N's geometry-preserving polar substitution).
