# Task C: AlphaFold vs Modeller vs crystal

## Files
- `AF-P04818-F1-model_v6.pdb` — AlphaFold prediction for human TYMS UniProt P04818
  (latest version 6, downloaded 2026-05-13).  Note: the URL in the
  Phase 7 spec referenced `_v4`; the AlphaFold-EBI cache no longer hosts that
  version, so we fetched `_v6` via the
  `https://alphafold.ebi.ac.uk/api/prediction/P04818` endpoint.
- `comparison.csv` — Ramachandran statistics + Cα RMSD vs 1HVY chain A,
  computed with both `cmd.super` (structure-based) and `cmd.align`
  (sequence-based) in PyMOL.
- `triple_overlay.png` — ray-traced PyMOL overlay (1600×1200).
- `../viewers/alphafold_overlay.html` — interactive 3Dmol.js viewer.

## How AlphaFold compares to a homology model when the crystal exists
The AlphaFold model and the Modeller homology model are independently
predicting the *same* sequence (P04818, residues 1–313).  The Modeller
model is templated on PDB entries that *include* 1HVY chain A, which is
why the Modeller-vs-1HVY Cα RMSD over the structurally aligned core is
typically <1 Å — Modeller essentially reproduces its template.  The
AlphaFold model, by contrast, was *not* templated on 1HVY at training
time in any way that is identifiable, yet it still recovers the
canonical TYMS fold to ~1 Å Cα RMSD over the well-modelled core.

## What AF's confidence (pLDDT) tells us about the active site
The bundled per-residue pLDDT scores in the B-factor column of the AF
PDB are global confidence values (0–100).  The Phase 5 active-site
panel — residues 50, 109, 175, 176, 195, 196, 214, 215, 225, 226, 256,
258 — is in the well-folded core of TYMS and AlphaFold reports very
high pLDDT (>90) at every one of these positions.  The disordered
N-terminal extension (residues 1–~25) carries low pLDDT (<70) and is
correctly predicted as flexible/disordered, so the AF model's
confidence map is consistent with the crystallographic reality.
