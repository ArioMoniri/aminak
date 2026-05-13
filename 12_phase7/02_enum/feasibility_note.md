# Feasibility note: enumeration of active-site mutations

We restrict the enumeration to the **14 residues** that line the dUMP /
folate active site of human TYMS (positions: 50, 109, 170, 175, 176, 195, 196, 214, 215, 216, 225, 226, 256, 258).

## Singles
- Count: 14 positions x 19 alternative AAs = **266 singles**.
- Per-mutation cost: build mutant model (PyMOL mutagenesis), prep PDBQT
  (PDB2PQR + Meeko), Vina dock 5x replicates @ exh 32 num_modes 20.
- Realistic wall-time on this 4-CPU laptop: ~5-10 s per single dock.
  Single-seed sweep: 266 x 5 s ~= 1330 s (22 min).
  5-replicate sweep: 6650 s (~110 min).
  Plus PyMOL mutagenesis + protonation pre-step: ~10 s/mutant -> ~44 min.
  Realistic: 1-2 hours for the full singles sweep with replicates.

## Doubles
- Count restricted to *pairs from this panel only*:
  C(14,2) x 19^2 = 91 x 361 = **32851 doubles**.
- At 5 s per single dock that is ~45.6 h for one seed each.
  With 5 replicates that becomes ~228.1 h (~9.5 days).
- **Infeasible** without GPU acceleration (e.g. Uni-Dock, Vina-GPU 2.1) or
  cluster batch.  Even GPU Vina at ~0.3 s/dock would be ~13.7 h.

## Why we still ship the IDs
Even without docking, the enumerated lists are useful for:
- prioritising biologically interesting subsets (charge-reversal, gain-of-aromatic,
  loss-of-thiol-nucleophile);
- cross-referencing against ClinVar / COSMIC / gnomAD population variants;
- guiding a smaller targeted sub-sweep (e.g. only the 56 singles whose
  functional_class is "charge_reversal" or "loss_of_charge").

## What we actually docked
Phase 5/7 dock 20 hand-picked mutants 1x and 8 priority mutants 5x.  Multi-replica
SD bracketed the per-mutation Vina noise floor (~0.2-0.5 kcal/mol typical;
see `12_phase7/01_replicas/multi_replica_results.csv`).
