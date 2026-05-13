# Phase 8: Beyond rigid Vina

This phase asks: *can we do better than the rigid-receptor Vina baseline?*
Three sub-tasks, all using AutoDock Vina 1.2.7 (the only Apple-Silicon-native
docking binary available; GNINA does not ship for arm64-darwin).

## 8a -- Alternative scoring functions (Vinardo + AD4)

Directory: `01_alt_scoring/`

Re-score the existing Vina top poses with Vina's other built-in scoring
functions:

- **Vinardo** (Quiroga & Villarreal 2016) is a refined empirical scoring
  function with stiffer hydrophobic and repulsive terms; it is known to
  rank actives vs decoys better than default Vina on several benchmarks.
- **AD4** is the AutoDock 4 force-field-style scoring.  In Vina 1.2 AD4
  scoring **requires precomputed affinity maps** (via `autogrid4`).  The
  `autogrid4` binary is not available on this Apple Silicon host (no
  Homebrew formula, no arm64 build from Scripps).  AD4 column is therefore
  empty in `alt_scoring_results.csv`; this is a documented limitation, not
  a script bug.

Driver: `scripts/v8/task_8a_alt_scoring.py`
Plots: `scripts/v8/task_8a_plots.py`

The mutant **apo** receptor PDBQTs were destructively overwritten by Vina
docking output in the v3 pipeline (see `scripts/v3/stage7_mutants_v3.py` --
`out_pdbqt` overwrites `rec_apo`).  8a regenerates them on-the-fly from
`{mut}_mut_h.pdb` using `obabel -xr -p 7.4 --partialcharge gasteiger`
(same recipe as v3's `prepare_receptor_with_charges` first method).

Outputs:
- `alt_scoring_results.csv`  -- 42 rows (WT_apo + WT_holo + 20 mut x apo,holo)
- `alt_scoring_compare.png`  -- static Vinardo-vs-Vina scatter
- `alt_scoring_compare.html` -- interactive Plotly scatter

## 8b -- Flexible-residue Vina re-dock (8 priority mutants)

Directory: `02_flexres/`

For each of the 8 priority mutants (same list as `scripts/v7/task_a_replicas.py`),
re-dock with the 14-residue active-site panel made flexible via Vina's
built-in `--flex` flag.

The standard tooling (`prepare_flexreceptor4.py` from MGLTools; meeko's
`mk_prepare_receptor.py -f`) was unusable on this dataset:
- MGLTools is not installed (no arm64 build).
- Meeko fails with `RuntimeError: Updated 1 H positions but deleted 7`
  on every PDB we tried (apo, holo, with and without explicit H, after
  prody clean-up, after pdb2pqr clean-up).

`scripts/v8/flex_split.py` is a self-contained replacement that emits the
exact Vina flex PDBQT format (`BEGIN_RES`/`ROOT`/`BRANCH`/`ENDBRANCH`/
`END_RES`).  It uses hardcoded chi-rotation topology templates for each
amino acid (ALA through VAL) and walks the templates to emit one BRANCH
per chi torsion.

Receptor build pipeline per mutant:
1. `obabel {mut}_mut_h.pdb -xr -p 7.4 --partialcharge gasteiger`
   -> apo PDBQT (chains A/B preserved; residues renumbered 1..N).
2. Shift residue numbers by +25 to restore source numbering
   (source PDBs start at residue 26 in chain A).
3. Concatenate cofactor PDBQT atoms (`06f_receptor_fixed/cofactor_A.pdbqt`
   and `cofactor_B.pdbqt`) to make holo PDBQT.
4. `flex_split.split_clean_pdbqt(holo_text, FLEX_PANEL)`
   -> rigid PDBQT (everything except the 14 panel side chains)
      + flex PDBQT (BEGIN_RES blocks for each panel side chain).
5. `vina --receptor rigid --flex flex --ligand dump.pdbqt --exhaustiveness 32
   --num_modes 20 --seed 42 --center/size {...same box as v5...}`

Outputs:
- `{mutant}_rigid.pdbqt`     -- rigid portion of the holo receptor
- `{mutant}_flexres.pdbqt`   -- the 14 panel side chains in flex format
- `{mutant}_flex.pdbqt`      -- Vina docking output (up to 20 modes)
- `{mutant}_flex.log`        -- Vina stdout/stderr
- `flexres_compare.csv`        -- summary: label, rigid_vina_score, flex_vina_score, delta_flex
- `flex_vs_rigid.png`/`.html`  -- scatter of flex vs rigid scores

## 8c -- Documentation + master comparison

`master_comparison.html` -- four-panel bar chart of holo scores across
rigid Vina, Vinardo, AD4 (empty), and flex Vina.

The diagnostic interpretation (C195A illusion, R215E charge reversal,
holo-Δ magnitudes) is printed to stdout at the end of
`scripts/v8/task_8c_integrate.py`.

## Why these three are "better than rigid Vina"

| Aspect                            | Rigid Vina | Vinardo | AD4 | Flex Vina |
|-----------------------------------|:---------:|:-------:|:---:|:---------:|
| Different energy function         |     -     |   yes   | yes |    -      |
| Better hydrophobic terms          |     -     |   yes   |  -  |    -      |
| Receptor side-chain rearrangement |     -     |    -    |  -  |   yes     |
| Same pose space                   |    yes    |   yes   | yes |    -      |
| Apple-Silicon-native              |    yes    |   yes   | no  |   yes     |

Each addresses a different failure mode of rigid Vina.  Vinardo asks
"would a better scoring function change the ranking on the *same* poses?"
Flex Vina asks "would induced fit of side chains allow a different pose
that the rigid receptor blocked?"  AD4 would have given a force-field
view, but is unavailable here.

## Limitations

- No induced-fit *backbone* movement (only side-chain rotamers in 8b).
- No proper continuum electrostatics (Vinardo and Vina both use a
  distance-dependent dielectric; R215E charge effects are still
  approximated, not Poisson-Boltzmann).
- The 14-residue flex panel adds ~30-40 rotatable degrees of freedom,
  pushing the conformational search hard; `--exhaustiveness 32` may be
  too low for full convergence.  A more thorough study would scan
  exhaustiveness or use multi-replica seeds (Phase 7's recipe).
- AD4 column omitted (no autogrid4 binary on this Apple Silicon host).
- The custom `flex_split.py` topology templates cover the 20 standard
  amino acids but assume canonical atom names (PDB v3.3 conventions);
  non-canonical residues or alternate names (e.g. legacy 1HB vs HB1) are
  best-effort.

## Diagnostic interpretation (auto-generated)

=== Phase 8 diagnostic interpretation ===

WT holo scores: vina=-7.407, vinardo=-5.112, ad4=None

1) C195A illusion check (was 'too good' under rigid Vina):
   rigid Vina (holo):    C195A -10.50  vs WT -7.41  Δ=-3.10
   Vinardo (holo):       C195A -8.52  vs WT -5.11  Δ=-3.41
   flex Vina (priority): C195A flex=-6.39  vs C195A rigid=-10.49  Δ=+4.10
   verdict: rigid Vina still shows C195A as 'better' (Δ=-3.10); Vinardo also says C195A is 'better' (Δ=-3.41) — illusion NOT fixed

2) R215E charge-reversal penalty check:
   rigid Vina: R215E -7.77  vs WT -7.41  Δ=-0.37
   Vinardo:    R215E -5.53  vs WT -5.11  Δ=-0.42
   flex Vina:  R215E flex -2.49

3) Holo Δ magnitudes (mean |Δ vs WT| across non-WT holo rows):
   rigid Vina mean |Δ|: 0.978 kcal/mol
   Vinardo   mean |Δ|: 1.396 kcal/mol
   Vinardo Δ ratio vs Vina = 1.43x (LARGER (more discrimination))
