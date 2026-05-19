# Phase 14g — HADDOCK3 for the Strategy-3 dimer-interface PPI question

HADDOCK3 is the right tool for the Phase-14 Strategy-3 question (does the
LR octapeptide block the TYMS dimer interface?) because:

  - it natively docks peptides ≥ 6 residues (Vina cannot, per Hassan 2017)
  - it accepts ambiguous interaction restraints (AIRs) — exactly the way
    a real medicinal chemist would specify "the peptide should touch one
    of these interface residues"
  - it returns a flexible, water-refined ensemble — not a single rigid pose

## Tooling required

```
# HADDOCK3 (BSD-licensed, Bonvin lab):
git clone https://github.com/haddocking/haddock3.git
cd haddock3 && pip install .

# CNS (Crystallography & NMR System) — required runtime dependency
# License: free academic, registration at https://cns-online.org/v1.3/
# Apple Silicon: compile from source (~20 min with --enable-aqua flag)
```

CNS is the non-trivial install — the legacy Fortran from Brünger's group.
The Bonvin lab maintains a build recipe at
https://github.com/haddocking/cns_solve_macos.

## The Phase-14 setup (committed configs)

### 1. Receptors and ligand
- Receptor:   `06f_receptor_fixed/protein_dimer_apo_fixed.pdbqt` (TYMS dimer)
- Peptide:    LR-octapeptide LSCQLYQR (Cardinale 2011, FEBS J 278:1487)
              built as PDB via RDKit `Chem.MolFromSequence` in Phase 14
              (`14_inhibitor_design/03_dimer_interface/ligands/LR8_LSCQLYQR.sdf`)

### 2. Active and passive residues for AIRs

Active residues on chain A (4 Å contact map output from Phase 14 A3):

  20, 21, 22, 23, 24, 25, 32, 33, 34, 35, 36, 37, 39,
  117, 118, 135, 148, 150, 151, 153, 157, 158, 159, 160,
  167, 168, 172, 173, 175, 177, 178, 199, 200, 202

Passive residues = chain-A residues within 6.5 Å of any active (HADDOCK
will treat these as restraint partners but with lower weight).

### 3. AIR (Ambiguous Interaction Restraint) generation

```bash
# Bonvin lab utility — ships with haddock3
haddock3-restraints active_passive_to_ambig \
    --active 20+21+22+23+24+25+32+33+34+35+36+37+39+117+118+135+148+150+151+153+157+158+159+160+167+168+172+173+175+177+178+199+200+202 \
    --target receptor \
    --output airs.tbl
```

### 4. config file — full HADDOCK3 pipeline

See `haddock3_config.cfg` in this folder. Three modules:
  - `topoaa` — generate topologies
  - `rigidbody` — initial docking (FFT-based)
  - `flexref` — semi-flexible refinement (peptide + interface side chains)
  - `mdref`   — short MD in explicit water
  - `caprieval` — CAPRI-style scoring (i_RMSD, l_RMSD, F_nat, DockQ)

### 5. Run

```bash
haddock3 haddock3_config.cfg     # ~3-6 hours on a single mac for one
                                  # peptide-receptor pair
```

### 6. Expected output structure

```
run1/
├── 0_topoaa/                    # built topologies
├── 1_rigidbody/                 # 1000 rigid-body decoys
├── 2_flexref/                   # 200 semi-flex refined
├── 3_mdref/                     # 50 explicit-water refined
├── 4_caprieval/                 # CAPRI scoring vs reference
└── analysis/
    ├── cluster.txt              # final clusters by HADDOCK score
    └── capri_ss.tsv             # per-pose CAPRI metrics
```

### 7. The scrambled control

Build a length-matched scrambled peptide (`numpy.random.default_rng(42).permutation`),
run the same HADDOCK3 protocol on it, and compare the cluster HADDOCK
scores. Specificity = top cluster score (canonical) − top cluster score (scrambled).

## Why we stop here

CNS + HADDOCK3 install + 2× ~6 h docking runs (canonical + scrambled) is
out of scope for one turn. The configs and restraint files are committed
for downstream execution.
