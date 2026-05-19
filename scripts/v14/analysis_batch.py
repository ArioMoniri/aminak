#!/usr/bin/env python3
"""Phase 14 — analysis batch:
  A3: Extend Smina rescoring to ALL Phase-7 mutants (full panel)
  A4: Modeller vs AlphaFold comparison plot
  A5: PPI / BSA dimer-interface analysis
  A6: Ramachandran outlier-reduction story (already-computed comparison)
  A7: MM-GBSA prep script (committed, not run)
  A8: HADDOCK3 config + restraint generation (committed, not run)
"""
from __future__ import annotations
import csv, json, subprocess, shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
OUT_FIG = REPO / "14_inhibitor_design" / "presentation" / "figures"
OUT_FIG.mkdir(parents=True, exist_ok=True)
ANALYSIS = REPO / "14_inhibitor_design" / "07_advanced_methods"
ANALYSIS.mkdir(parents=True, exist_ok=True)
SMINA = shutil.which("smina") or "/opt/homebrew/bin/smina"


# ────────────────────────────────────────────────────────────────
# A4: Modeller vs AlphaFold — load existing summaries, render plot
# ────────────────────────────────────────────────────────────────
def a4_modeller_vs_alphafold():
    src = REPO / "12_phase7" / "03_alphafold" / "comparison.csv"
    if not src.exists():
        print("  ! Phase 7c comparison.csv missing")
        return
    rows = list(csv.DictReader(src.open()))
    print(f"  Phase 7c comparison rows: {len(rows)} — columns: {list(rows[0].keys()) if rows else []}")
    # Save copy
    (ANALYSIS / "modeller_vs_alphafold.csv").write_text(src.read_text())

    # Plot: per-model RMSD vs 1HVY chain A + Lovell %favoured side-by-side
    # We use comparison.csv if it has the right columns
    models = []
    rmsd_vs_crystal = []
    rama_favoured = []
    for r in rows:
        name = r.get("model", "") or r.get("structure", "")
        rmsd_ka = r.get("rmsd_kabsch_align_A") or r.get("rmsd_align_A") or r.get("rmsd_align")
        fav = r.get("rama_pct_favored") or r.get("ramachandran_favored") or r.get("lovell_favored")
        if name and rmsd_ka:
            try:
                models.append(name)
                rmsd_vs_crystal.append(float(rmsd_ka))
                rama_favoured.append(float(fav) if fav else None)
            except ValueError: continue
    if not models:
        # try to construct from the per-model lovell stats files
        af_csv = REPO / "12_phase7" / "03_alphafold" / "summary_alphafold_v6.csv"
        mod_3 = REPO / "12_phase7" / "03_alphafold" / "summary_modeller_B99990003.csv"
        mod_10 = REPO / "12_phase7" / "03_alphafold" / "summary_modeller_B99990010.csv"
        for label, p in [("AlphaFold v6", af_csv), ("Modeller B99990003", mod_3), ("Modeller B99990010", mod_10)]:
            if not p.exists(): continue
            r = next(csv.DictReader(p.open()), None)
            if r is None: continue
            models.append(label)
            for k in r:
                if "rmsd" in k.lower() and "align" in k.lower():
                    try: rmsd_vs_crystal.append(float(r[k])); break
                    except: pass
            else:
                rmsd_vs_crystal.append(0)
            for k in r:
                if "favor" in k.lower():
                    try: rama_favoured.append(float(r[k])); break
                    except: pass
            else:
                rama_favoured.append(None)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    colours = ["#e6553f" if "AlphaFold" in m else "#3a86c8" for m in models]
    ax1.barh(range(len(models)), rmsd_vs_crystal, color=colours)
    ax1.set_yticks(range(len(models))); ax1.set_yticklabels(models, fontsize=10)
    ax1.set_xlabel("Cα RMSD vs 1HVY chain A (Å)  — lower = closer to crystal")
    ax1.set_title("Backbone fidelity vs the experimental crystal")
    ax1.grid(True, axis="x", alpha=0.3)
    for i, v in enumerate(rmsd_vs_crystal):
        ax1.annotate(f"{v:.2f} Å", xy=(v, i), xytext=(4, 0), textcoords="offset points",
                     va="center", fontsize=9, color="black")

    if any(rama_favoured):
        favs = [f if f else 0 for f in rama_favoured]
        ax2.barh(range(len(models)), favs, color=colours)
        ax2.set_yticks(range(len(models))); ax2.set_yticklabels(models, fontsize=10)
        ax2.set_xlabel("Ramachandran % favoured  (Lovell 4-map)")
        ax2.set_title("Local geometry quality")
        ax2.axvline(x=92.2, ls="--", color="#888", alpha=0.6, label="1HVY crystal (92.2%)")
        ax2.set_xlim(0, 100)
        ax2.grid(True, axis="x", alpha=0.3); ax2.legend(loc="lower right", fontsize=9)
        for i, v in enumerate(favs):
            ax2.annotate(f"{v:.1f} %", xy=(v, i), xytext=(4, 0), textcoords="offset points",
                         va="center", fontsize=9, color="black")
    else:
        ax2.text(0.5, 0.5, "Lovell %favoured not in Phase-7c summary", ha="center", va="center")
        ax2.axis("off")

    plt.tight_layout()
    plt.savefig(OUT_FIG / "modeller_vs_alphafold.png", dpi=140, facecolor="white")
    plt.close()
    print(f"  → {OUT_FIG / 'modeller_vs_alphafold.png'}")


# ────────────────────────────────────────────────────────────────
# A5: PPI / BSA dimer-interface analysis
# ────────────────────────────────────────────────────────────────
def a5_ppi_bsa():
    import freesasa
    apo_pdb = REPO / "06f_receptor_fixed" / "dimer_noH.pdb"
    if not apo_pdb.exists():
        print("  ! dimer_noH.pdb missing"); return

    # Need separate chains A and B as standalone PDBs for delta-SASA
    text = apo_pdb.read_text().splitlines()
    n_atoms_total = sum(1 for ln in text if ln.startswith("ATOM"))
    half = n_atoms_total // 2
    chA_lines = []; chB_lines = []; cnt = 0
    for ln in text:
        if ln.startswith("ATOM"):
            (chA_lines if cnt < half else chB_lines).append(ln)
            cnt += 1
    chA_pdb = ANALYSIS / "chainA_only.pdb"; chB_pdb = ANALYSIS / "chainB_only.pdb"
    chA_pdb.write_text("\n".join(chA_lines) + "\nEND\n")
    chB_pdb.write_text("\n".join(chB_lines) + "\nEND\n")

    # Compute SASA for each
    def sasa_of(p):
        s = freesasa.calc(freesasa.Structure(str(p)))
        return float(s.totalArea())

    sA = sasa_of(chA_pdb); sB = sasa_of(chB_pdb); sAB = sasa_of(apo_pdb)
    bsa = (sA + sB - sAB) / 2.0   # standard BSA definition
    print(f"  SASA chain A: {sA:.1f} Å²")
    print(f"  SASA chain B: {sB:.1f} Å²")
    print(f"  SASA dimer:   {sAB:.1f} Å²")
    print(f"  → BSA (buried surface area, per side): {bsa:.1f} Å²")

    # Per-residue BSA: residue SASA in (A or B alone) − residue SASA in dimer
    # use freesasa.classify_by_residue
    def per_res(p, chain_label):
        s = freesasa.calc(freesasa.Structure(str(p)))
        out = {}
        ra = s.residueAreas()
        # ra is dict-of-dict {chain: {resid: ResidueArea}}
        for chain in ra:
            for resid_str, area in ra[chain].items():
                out[(chain_label, int(resid_str))] = float(area.total)
        return out

    apo_per = per_res(apo_pdb, "AB")
    a_alone = per_res(chA_pdb, "A")
    b_alone = per_res(chB_pdb, "B")

    # Walk residues — match by resid (chain id is artificial because the source PDB has only chain A label)
    rows = []
    # we treat the first-half resids as chain A and second-half as chain B
    half_resid = sorted(set(r for _, r in a_alone))
    for (_, r) in sorted(a_alone):
        in_dimer = apo_per.get(("AB", r))
        alone = a_alone[("A", r)]
        if in_dimer is None: continue
        delta = alone - in_dimer
        rows.append({"chain":"A", "resid": r, "SASA_alone": alone, "SASA_in_dimer": in_dimer,
                     "delta_SASA": delta})
    for (_, r) in sorted(b_alone):
        in_dimer = apo_per.get(("AB", r))
        alone = b_alone[("B", r)]
        if in_dimer is None: continue
        delta = alone - in_dimer
        rows.append({"chain":"B", "resid": r, "SASA_alone": alone, "SASA_in_dimer": in_dimer,
                     "delta_SASA": delta})

    csv_path = ANALYSIS / "ppi_per_residue_bsa.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"  → {csv_path}  ({len(rows)} rows)")

    # Top hot-spots — residues with the largest BSA contribution (per side)
    rows_sorted = sorted(rows, key=lambda r: -r["delta_SASA"])
    top_a = [r for r in rows_sorted if r["chain"]=="A"][:15]
    top_b = [r for r in rows_sorted if r["chain"]=="B"][:15]
    print(f"\n  Top hot-spots (chain A, by Δ SASA on dimerisation):")
    for r in top_a[:10]:
        print(f"    A/{r['resid']:>4d}  Δ {r['delta_SASA']:>6.1f} Å²")

    # Plot — bar chart of top 20 per side
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, rows_, title in [(axes[0], top_a, "Chain A — top 15 dimer-interface hot-spots"),
                              (axes[1], top_b, "Chain B — top 15 dimer-interface hot-spots")]:
        rs = [r["resid"] for r in rows_]
        ds = [r["delta_SASA"] for r in rows_]
        ax.barh(range(len(rs)), ds, color="#3a86c8")
        ax.set_yticks(range(len(rs))); ax.set_yticklabels([str(r) for r in rs], fontsize=9)
        ax.set_xlabel("Δ SASA on dimerisation (Å²)  — larger = more buried")
        ax.set_title(title, fontsize=11)
        ax.invert_yaxis()
        ax.grid(True, axis="x", alpha=0.3)
        for i, v in enumerate(ds):
            ax.annotate(f"{v:.1f}", xy=(v, i), xytext=(4, 0), textcoords="offset points",
                        va="center", fontsize=8)
    plt.suptitle(f"TYMS dimer interface — total BSA = {bsa:.0f} Å² per side", fontsize=13)
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    plt.savefig(OUT_FIG / "ppi_dimer_interface.png", dpi=140, facecolor="white")
    plt.close()
    print(f"  → {OUT_FIG / 'ppi_dimer_interface.png'}")

    # write summary JSON
    summary = {
        "sasa_chainA_alone_A2": round(sA, 1),
        "sasa_chainB_alone_A2": round(sB, 1),
        "sasa_dimer_A2": round(sAB, 1),
        "bsa_per_side_A2": round(bsa, 1),
        "bsa_total_A2": round(2 * bsa, 1),
        "top_a_resids": [{"resid": r["resid"], "delta_sasa": round(r["delta_SASA"], 1)} for r in top_a[:15]],
        "top_b_resids": [{"resid": r["resid"], "delta_sasa": round(r["delta_SASA"], 1)} for r in top_b[:15]],
    }
    (ANALYSIS / "ppi_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"  → {ANALYSIS / 'ppi_summary.json'}")


# ────────────────────────────────────────────────────────────────
# A6: Ramachandran improvement — surface the existing comparison
# ────────────────────────────────────────────────────────────────
def a6_ramachandran():
    src = REPO / "10b_modeller_refined" / "04_refined_lovell" / "comparison_before_after.png"
    if src.exists():
        shutil.copyfile(src, OUT_FIG / "ramachandran_before_after.png")
        print(f"  → {OUT_FIG / 'ramachandran_before_after.png'} (copied from Phase 6b)")
    # also save the stats CSV reference
    stats = REPO / "10b_modeller_refined" / "04_refined_lovell" / "lovell_stats_best_refined.csv"
    if stats.exists():
        shutil.copyfile(stats, ANALYSIS / "ramachandran_best_refined_stats.csv")
        print(f"  → {ANALYSIS / 'ramachandran_best_refined_stats.csv'}")


# ────────────────────────────────────────────────────────────────
# A7: MM-GBSA prep — write a ready-to-run script (do not run AmberTools)
# ────────────────────────────────────────────────────────────────
MMGBSA_README = """# Phase 14f — MM-GBSA rescoring (next step beyond Smina)

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
ante-MMPBSA.py -p complex.prmtop -c complex_nowat.prmtop -r receptor.prmtop -l ligand.prmtop \\
                -s :WAT,Na+,Cl-

# 5. ΔΔG vs WT via MMPBSA.py (single-trajectory protocol)
MMPBSA.py -O -i mmpbsa.in -o mmpbsa.dat -sp complex.prmtop \\
           -cp complex_nowat.prmtop -rp receptor.prmtop -lp ligand.prmtop \\
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
"""

TLEAP_IN = """# tleap.in — solvate complex for MM-GBSA / MM-PBSA
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2
loadamberparams mutant_lig.frcmod
lig = loadmol2 mutant_lig.mol2
rec = loadpdb mutant_apo.pdb         # apo receptor (no ligand, charged via pdb2pqr)
complex = combine {rec lig}
solvatebox complex TIP3PBOX 12.0
addions complex Na+ 0
addions complex Cl- 0
saveamberparm rec receptor.prmtop receptor.inpcrd
saveamberparm lig ligand.prmtop ligand.inpcrd
saveamberparm complex complex.prmtop complex.inpcrd
savepdb complex complex_solvated.pdb
quit
"""

MMPBSA_IN = """# mmpbsa.in — GB + PB single-trajectory protocol
&general
  startframe=1, endframe=1000, interval=10, verbose=2,
&end
&gb
  igb=5, saltcon=0.150,         # generalised Born, 150 mM physiological salt
&end
&pb
  istrng=0.150, fillratio=4.0,  # Poisson-Boltzmann at 150 mM
&end
"""

def a7_mmgbsa_setup():
    p = ANALYSIS / "mmgbsa"
    p.mkdir(parents=True, exist_ok=True)
    (p / "README.md").write_text(MMGBSA_README)
    (p / "tleap.in").write_text(TLEAP_IN)
    (p / "mmpbsa.in").write_text(MMPBSA_IN)
    print(f"  → {p}/README.md + tleap.in + mmpbsa.in")


# ────────────────────────────────────────────────────────────────
# A8: HADDOCK3 step-by-step config + restraint generation
# ────────────────────────────────────────────────────────────────
HADDOCK3_README = """# Phase 14g — HADDOCK3 for the Strategy-3 dimer-interface PPI question

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
haddock3-restraints active_passive_to_ambig \\
    --active 20+21+22+23+24+25+32+33+34+35+36+37+39+117+118+135+148+150+151+153+157+158+159+160+167+168+172+173+175+177+178+199+200+202 \\
    --target receptor \\
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
"""

HADDOCK3_CFG = """# haddock3_config.cfg — TYMS dimer + LR octapeptide
run_dir = "run1_TYMS_LR8"
postprocess = true
clean = false

molecules = [
    "receptor.pdb",
    "lig_LR8.pdb",
]

[topoaa]

[rigidbody]
ambig_fname = "airs.tbl"
sampling = 1000

[flexref]
ambig_fname = "airs.tbl"

[mdref]
ambig_fname = "airs.tbl"

[caprieval]
reference_fname = "receptor_chain_a_b_only.pdb"
"""

def a8_haddock3_setup():
    p = ANALYSIS / "haddock3"
    p.mkdir(parents=True, exist_ok=True)
    (p / "README.md").write_text(HADDOCK3_README)
    (p / "haddock3_config.cfg").write_text(HADDOCK3_CFG)
    # active-residue list as a separate file for clarity
    (p / "active_residues_chainA.txt").write_text(
        "20 21 22 23 24 25 32 33 34 35 36 37 39\n"
        "117 118 135 148 150 151 153 157 158 159 160\n"
        "167 168 172 173 175 177 178 199 200 202\n"
    )
    print(f"  → {p}/README.md + haddock3_config.cfg + active_residues_chainA.txt")


# ────────────────────────────────────────────────────────────────
# A3: Extend Smina rescoring to ALL Phase-7 mutants
# ────────────────────────────────────────────────────────────────
def a3_extended_smina():
    """Smina-rescore every Phase-7 stripped pose (not just the subset)."""
    import re
    stripped = REPO / "13_phase8" / "01_alt_scoring" / "stripped_poses"
    apo = REPO / "06f_receptor_fixed" / "protein_dimer_apo_fixed.pdbqt"
    custom_q = REPO / "14_inhibitor_design" / "06_smina_rescore" / "custom_scoring_q.txt"
    q_amp = REPO / "14_inhibitor_design" / "06_smina_rescore" / "custom_scoring_qamp.txt"

    rows = []
    for mut_dir in sorted(stripped.iterdir()):
        if not mut_dir.is_dir(): continue
        pose = next(mut_dir.glob("*.pdbqt"), None)
        if not pose: continue
        # need MODEL 1 only
        # (the stripped poses might be MODEL 1 already, but be safe)
        scores = {}
        for sname, sarg, scust in [("vina","vina",None),
                                    ("vinardo","vinardo",None),
                                    ("custom_q",None,custom_q),
                                    ("q_amp",None,q_amp)]:
            cmd = [SMINA, "--score_only", "-r", str(apo), "-l", str(pose)]
            if sarg: cmd += ["--scoring", sarg]
            if scust: cmd += ["--custom_scoring", str(scust)]
            try:
                r = subprocess.run(cmd, capture_output=True, timeout=60)
                m = re.search(r"Affinity:\s*([-\d.]+)", r.stdout.decode())
                scores[sname] = float(m.group(1)) if m else None
            except subprocess.TimeoutExpired:
                scores[sname] = None
        rows.append({"mutant": mut_dir.name, "pose": str(pose.relative_to(REPO)),
                     **scores})

    csv_path = ANALYSIS / "smina_full_panel.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"  → {csv_path}  ({len(rows)} mutants × 4 scorers)")

    # Plot — full-panel comparison
    rows_sorted = sorted(rows, key=lambda r: r["vina"] if r["vina"] else 0)
    labels = [r["mutant"].replace("_holo","").replace("_apo","") for r in rows_sorted]
    vina = [r["vina"] or 0 for r in rows_sorted]
    vinardo = [r["vinardo"] or 0 for r in rows_sorted]
    custom = [r["custom_q"] or 0 for r in rows_sorted]
    qamp = [r["q_amp"] or 0 for r in rows_sorted]
    x = np.arange(len(labels)); width = 0.20
    fig, ax = plt.subplots(figsize=(15, 6))
    ax.bar(x - 1.5*width, vina,    width, label="Vina",       color="#bdc3c7")
    ax.bar(x - 0.5*width, vinardo, width, label="Vinardo",    color="#7f8c8d")
    ax.bar(x + 0.5*width, custom,  width, label="custom_q (Vina + 0.3·elec + 0.1·desolv)", color="#e67e22")
    ax.bar(x + 1.5*width, qamp,    width, label="q_amp (10× electrostatic)", color="#c0392b")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Smina top-pose score (kcal/mol)  — more negative = tighter")
    ax.set_title("Phase 14e extended — full Phase-7 panel rescored with Smina ± electrostatics")
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(OUT_FIG / "smina_full_panel.png", dpi=140, facecolor="white")
    plt.close()
    print(f"  → {OUT_FIG / 'smina_full_panel.png'}")


if __name__ == "__main__":
    print("=== Phase 14 — analysis batch ===")
    print("\n--- A3: Extended Smina rescoring (full Phase-7 panel) ---")
    a3_extended_smina()
    print("\n--- A4: Modeller vs AlphaFold ---")
    a4_modeller_vs_alphafold()
    print("\n--- A5: PPI / BSA dimer interface ---")
    a5_ppi_bsa()
    print("\n--- A6: Ramachandran improvement (Phase 6b) ---")
    a6_ramachandran()
    print("\n--- A7: MM-GBSA setup (configs committed, not run) ---")
    a7_mmgbsa_setup()
    print("\n--- A8: HADDOCK3 setup (configs committed, not run) ---")
    a8_haddock3_setup()
    print("\nDone.")
