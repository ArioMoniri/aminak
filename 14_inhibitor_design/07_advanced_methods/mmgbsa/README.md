# Phase 14f — MM-GBSA rescoring (next step beyond Smina)

The Phase 14e Smina rescoring (Vina + electrostatic + desolvation, even at 10×
electrostatic weight) **could not resolve the R→E charge-reversal mutants**
because the rigid dUMP pose's phosphate is >6 Å from the Arg/Glu side-chain
centre and Smina's electrostatic term has an 8 Å cutoff. The fix needs:

  1. structural relaxation of the side chain at the mutated position
  2. Poisson-Boltzmann (or generalised-Born) electrostatic treatment
  3. thermodynamic-cycle ΔΔG vs WT

This is what AmberTools' MMPBSA.py provides.

## Tooling required

```
brew install --cask amber-tools     # NOT a Homebrew formula; obtain from
                                    # https://ambermd.org/AmberTools.php
                                    # (free, registration required)
# OR via conda:
conda install -c conda-forge ambertools=23
```

AmberTools is ~3 GB; license is free for academic use. On arm64-darwin
the conda-forge bottle compiles cleanly.

## Prep pipeline (per pose)

```
# 0. extract MODEL 1 of each docked PDBQT to a clean PDB
obabel mutant.pdbqt -O mutant_lig.pdb

# 1. parametrise the ligand with antechamber (GAFF2)
antechamber -i mutant_lig.pdb -fi pdb -o mutant_lig.mol2 -fo mol2 -c bcc -nc 0
parmchk2 -i mutant_lig.mol2 -f mol2 -o mutant_lig.frcmod

# 2. tleap — build receptor + ligand + solvate (TIP3P) + neutralise (Na+/Cl-)
tleap -f tleap.in    # tleap.in below

# 3. minimise → heat 0->300K → 1 ns NPT equilibration  → 10 ns NVT production
pmemd.cuda -O -i min.in   -p complex.prmtop -c complex.inpcrd -o min.out  -r min.rst
pmemd.cuda -O -i heat.in  -p complex.prmtop -c min.rst       -o heat.out -r heat.rst
pmemd.cuda -O -i eq.in    -p complex.prmtop -c heat.rst      -o eq.out   -r eq.rst
pmemd.cuda -O -i prod.in  -p complex.prmtop -c eq.rst        -o prod.out -r prod.rst -x prod.nc

# 4. extract ligand-only and receptor-only topologies (no waters/ions)
ante-MMPBSA.py -p complex.prmtop -c complex_nowat.prmtop -r receptor.prmtop -l ligand.prmtop \
                -s :WAT,Na+,Cl-

# 5. ΔΔG vs WT via MMPBSA.py (single-trajectory protocol)
MMPBSA.py -O -i mmpbsa.in -o mmpbsa.dat -sp complex.prmtop \
           -cp complex_nowat.prmtop -rp receptor.prmtop -lp ligand.prmtop \
           -y prod.nc
```

## Per-mutant cost on this hardware

On arm64-darwin **without CUDA** (pmemd.cuda needs an NVIDIA GPU): use
sander (CPU). Per-mutant wall-time estimate:
- minimise:     ~5 min
- heat + eq:    ~30 min
- 10 ns prod:   ~24 h on CPU (sander) — feasible per mutant but expensive
- MMPBSA.py:    ~30 min for GB; PB adds ~2 h

So for the 9 mutants × WT = 10 systems × 24 h = ~10 days on a single arm64
Mac without a GPU. With a CUDA box (RTX 4090, ~50× speedup) the whole
sweep is overnight.

## Why we stop here (signposted)

This is the right next step but is **out of scope for a single-laptop
arm64-darwin pipeline within this turn**. The ready-to-run scripts are
committed for downstream execution on appropriate hardware.

See: AmberTools MMPBSA.py manual §3, Wang et al. 2019 (J. Chem. Inf. Model.).
