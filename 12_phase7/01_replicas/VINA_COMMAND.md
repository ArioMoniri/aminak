# Exact Vina invocation used for Task A

All 50 docks (10 conditions x 5 seeds) used identical box / scoring
parameters.  The literal CLI for one example — `R215A_N226A_holo`, seed 42 — is:

```bash
/opt/homebrew/bin/vina \
  --receptor /Users/ario/conserved_site_project/07e_mut_docking_v5/R215A_N226A/R215A_N226A_holo.pdbqt \
  --ligand   /Users/ario/conserved_site_project/05b_ligand_v2/dump.pdbqt \
  --center_x -0.137 --center_y 4.232 --center_z 15.159 \
  --size_x   18.0   --size_y   18.0  --size_z   18.0 \
  --exhaustiveness 32 \
  --num_modes      20 \
  --seed           42 \
  --out /Users/ario/conserved_site_project/12_phase7/01_replicas/raw/R215A_N226A_holo_seed42.pdbqt
```

## Box parameters

| Parameter         | Value                                 |
|-------------------|---------------------------------------|
| Centroid x        | -0.137                                |
| Centroid y        | 4.232                                 |
| Centroid z        | 15.159                                |
| Box size x/y/z    | 18.0 A x 18.0 A x 18.0 A              |
| Exhaustiveness    | 32                                    |
| num_modes         | 20                                    |
| Seeds             | 42, 7, 13, 99, 256                    |
| CPU               | 4 (auto)                              |
| Receptor format   | PDB2PQR-corrected -> Meeko PDBQT      |
| Ligand            | dUMP, 5 active torsions (heavy=20)    |

The centroid is the chain-A active-site Cα centroid carried over from
the v3/v4/v5 canonical pipeline (the v5 mutant runs that produced
`07e_mut_docking_v5/mutant_results_v5.csv`).  The 18-A cubic box is
slightly smaller than the v5 22-A canonical box; this is intentional —
Phase 7 is specifically asking for stochastic noise floor of a tight
box around the dUMP pocket, not a wide rescore.

## Reproducing one run from scratch

1. Source the environment: `source 00_setup/env.sh`.
2. Confirm the receptor and ligand PDBQTs exist (paths above).
3. Run the CLI exactly as shown, keeping `--out` so the docked poses are
   captured.  The `--seed` flag fully determines the run, so identical
   commands are bit-reproducible on the same Vina build.
4. Parse the top affinity from stdout (line starting with `   1` after
   the dashed table header).
