# V4 Structural Bioinformatics Audit

Round 4. Date: 2026-05-12.

## Verdict: CONDITIONAL PASS

The v3 FAIL items (cofactor protonation, atom-name preservation) are genuinely fixed. Multi-seed sweep, n_modes, and affinity-based selection are all correct. **However, the WT-holo result (top affinity −5.24, RMSD 12.9 Å) reflects a methodological flaw, not biology.** It must not be reported as the canonical holo affinity without large caveats — or should be redone.

## Per-audit findings

1. **Cofactor protonation is real.** MD5 differs from v2; both carboxylate Os have no H within 1.3 Å; provenance JSON confirms canonical SMILES with two `[O-]` carboxylates, net charge −2.
2. **Holo receptor rebuilt** — `protein_dimer_holo.pdbqt` mtime 2026-05-12 13:29 (newer than v3 12:56). Contains v4 cofactor coords.
3. **Atom-name preservation works** — `wt_apo_top_pose.pdb` contains all 20 crystal heavy atom names; named-pose RMSD = 0.9116 Å matches JSON's `rmsd_top_named` exactly.
4. **Multi-seed sweep verified** — 5 seeds for both wt_apo and wt_holo; n_modes range 27–31.
5. **Affinity-based selection** — `key=lambda s: (top, -n_modes)` confirmed; no `key=rmsd_to_native`.
6. **n_modes for WT holo = 30**, ≫ 3 (v3).

## Scientific consequence (the central question)

The −5.24 kcal/mol holo result is **not credibly an electrostatic-repulsion story**. The 1HVY crystal shows dUMP and raltitrexed coexisting in the active site; a "−2 cofactor repels −2 dUMP" framing is inconsistent with experiment.

**The dominant artefact is COFACTOR PLACEMENT, not protonation.** The placed cofactor has a **2.71 Å heavy-atom RMSD vs the 1HVY-bound conformer** because Kabsch superimposed the CCD-ideal D16 onto the bound conformer. Combined with the now correct double-deprotonation, the misplaced glutamate γ-carboxylate is sitting near the dUMP phosphate region and blocking it. There is also a real cofactor-protein clash (cofA O1 ↔ PHE 80 CD2 at 1.95 Å) confirming the placement is non-physical. **This is a placement problem, not a protonation problem; the v4 fix exposed it rather than created it.**

This should be a **blocker for any "holo affinity" claim**. Acceptable framings:
- Report apo only (−9.20 kcal/mol, 0.91 Å) as the credible WT result.
- Mark holo as "not interpretable: cofactor placement RMSD 2.71 Å introduces steric/electrostatic artifacts in the dUMP pocket".

## Specific corrections still needed

1. **Re-place cofactor by structural superposition of the TS protein backbone, not Kabsch on the ligand atoms.** Use the ORIGINAL `cofactor_chainA.pdb` from 1HVY (already correctly bound) and only re-protonate it IN PLACE — do not re-derive coordinates from the CCD ideal SDF. Add Hs to existing crystal heavy-atom coords; remove COOH Hs (so carboxyls become COO⁻); never move any heavy atom.
2. **Treat WT holo affinity as undefined** in the v4 report. Re-run holo only after fix #1.
3. **Re-derive mutant holo Δ values** after fix #1.

## Lower-priority improvements
- Stale `WT apo named-pose RMSD = 3.802` log line should be removed.
- `_common_v4.restore_atom_names` falls back silently to `05b_ligand_v2/dump.pdbqt` — log when fallback fires.
- Cofactor-protein clash check should be added to stage 3 as a gate.
