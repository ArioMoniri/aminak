# Task D: SASA per residue and correlation with Vina

## Method
- `freesasa` 2.2.1 (Python binding) on chain A only.
- **Classifier**: `freesasa.Classifier()` default = NACCESS-style protein radii (Lee-Richards / Ali et al. 1983 atom-radii table). Default probe radius = 1.4 A. Default Lee-Richards algorithm with 20 slices/atom.
- **Hydrogens**: stripped by the default protein classifier (verified empirically: 0/4650 H atoms are loaded from `protein_dimer_h.pdb`). This matches the AD4 / Naccess / DSSP convention of *heavy-atom-only* SASA.
- WT input: `protein_dimer_h.pdb` (the protonated dimer used everywhere downstream — H atoms present in the file but ignored by the classifier).
- Mutant inputs: each `<mut>_mut_h.pdb` from `07e_mut_docking_v5/<mut>/`.
- Per-residue values are summed over all heavy atoms in that residue.
- "6 A neighbours" in the README sentence is a *sequence-distance* proxy (residues at positions +/-3 from the mutated site), not a 3D contact map. This is documented in `task_d_sasa.py` and is a deliberate simplification.

## Output schema
`sasa_<mut>.csv` columns: `residue_position, sasa_A2, wt_sasa_A2, dsasa_A2`.
`sasa_vs_dvina.csv` aggregates per mutant:
- `dsasa_at_mut_A2`: SASA change *at* the mutated residue(s).
- `dsasa_neigh_pm3_A2`: SASA change at +/- 3 sequence neighbours.
- `dsasa_total_focus_A2`: sum of the above.
- `delta_vina_vs_wt`: from `07e_mut_docking_v5/mutant_results_v5.csv`.

## Result
Pearson r between dSASA(focus) and dVina(holo) across 20 mutants:
**r = -0.192**.

No strong monotonic SASA-Vina relationship in this set.

## Interpretation
For the alanine-scan mutants in this panel, removing a bulky sidechain
opens the pocket, but the binding affinity change depends on whether the
removed sidechain was donating productive H-bonds / electrostatics with
dUMP.  C195A is the cleanest example: nucleophile-knockout that *also*
removes a bulky thiol -> the dUMP pyrimidine slides deeper, picking up
~2.2 kcal/mol of "fake" affinity that does not reflect catalytic
competence.  R215A/R215E lose the dUMP-phosphate salt bridge and the
pocket opens up but the affinity goes the wrong way (less negative)
because Vina rewards the lost H-bond more than it credits the new
breathing room.

So the SASA -> dVina correlation here is **-0.19**: not flat,
but not deterministic either.  Pocket geometry alone does not predict
ligand affinity in this enzyme; specific polar contacts (Arg phosphate
clamps, Asn226 ribose H-bond, Tyr258 stacking) dominate.
