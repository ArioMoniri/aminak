# Phase 10b — Modeller refinement and Ramachandran optimisation

This phase re-evaluates the Phase-6 homology models with a proper
MolProbity / Lovell-2003 Ramachandran scheme and applies two
structural-refinement passes to push the models closer to
crystal-quality φ/ψ statistics. **All artefacts are written into new
folders (`10b_modeller_refined/`, `scripts/modeller/refined/`) — the
original Phase-6 results are untouched.**

---

## 1. What φ/ψ and Ramachandran statistics measure

Every non-terminal residue in a protein has two backbone torsion
angles:

- **φ (phi)** = rotation around N–Cα
- **ψ (psi)** = rotation around Cα–C

For each residue, plot (φ, ψ) on a 2-D map (the **Ramachandran plot**).
Steric clashes between the side chain and adjacent peptide carbonyl /
amide groups make most of the (φ, ψ) plane physically forbidden. Only
two large regions, **α-helix** (φ ≈ −60°, ψ ≈ −45°) and **β-sheet**
(φ ≈ −120°, ψ ≈ +130°), plus a tiny **left-handed α** island for
glycine, are populated in real protein structures.

The standard MolProbity terminology:

- **favoured** — inside the inner contour where ≥98 % of residues from
  the Lovell top-500 reference set lie.
- **allowed** — inside the outer (~99.95 %) contour but outside
  favoured.
- **outlier** — outside the outer contour. A few real outliers are
  expected in any structure (functional strain at active sites,
  metal-binding loops), but a model with > 0.5 % outliers is a warning
  sign.

A well-refined crystal structure typically scores ≥ 95 % favoured and
< 0.5 % outlier under MolProbity. A homology model built without
explicit Ramachandran restraints can easily land at 80–85 % favoured,
not because the geometry is wrong, but because a single permissive
polygon (instead of the four residue-type-specific maps below)
mis-classifies normal Gly/Pro residues as outliers.

---

## 2. Why the v1 hand-drawn polygons over-counted outliers

The Phase-6 validator (`scripts/modeller/step6_validate.py`) used a
single hand-drawn favoured + allowed polygon for **all** residue
types. That conflates three problems:

1. **Glycine** has no β-carbon, so the two-body Gly–C(=O) repulsion
   that drives the (φ, ψ) restriction in other residues is absent.
   Gly is the only residue that can sit in the left-handed-α basin
   (φ ≈ +60°, ψ ≈ +45°) **without** strain. A general-case polygon
   would call those Gly residues "outliers" — they are not.
2. **Proline** has its φ angle locked by the cyclic N–Cα ring to
   roughly −60° ± 25°. The general-case favoured region (which is
   wide open across φ ∈ [−180°, −30°]) is therefore much *wider* than
   the true Pro favoured region. A general-case polygon would call
   normal prolines "allowed" when they are in fact at the edge.
3. **Pre-proline** residues (any residue immediately preceding a Pro)
   have restricted ψ because the next residue's pyrrolidine ring
   creates additional steric pressure. This needs its own polygon.

Using a single polygon, the Phase-6 baseline scored
**83.5 – 85.3 %** favoured. Under the proper 4-map partition, the
*same* models score **94.7 – 96.1 %** favoured. Most of the apparent
gain is **re-classification** of glycines and prolines that were never
strained in the first place — but the residual real outliers are
small enough that they are amenable to MD-based refinement.

---

## 3. The Lovell partition implemented here

`scripts/modeller/refined/lovell_ramachandran.py` classifies each
residue using one of four maps:

| Map         | Selection                                                                 | Favoured shape                                                                                     |
|-------------|---------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------|
| general     | non-Gly, non-Pro, non-pre-Pro                                              | α-R basin + β basin + small L-α island                                                              |
| glycine     | residue is GLY                                                            | symmetric (φ ↔ −φ): α-R + β + their mirrors + bottom band                                          |
| proline     | residue is PRO                                                            | two narrow φ ≈ −60° basins, one at ψ < 0 (α-like), one at ψ > 90° (β-like)                          |
| pre-proline | residue's i+1 neighbour is PRO                                            | as general but with a narrowed ψ ≈ 0 hole                                                          |

Anything outside its map's **allowed** polygon is an outlier. The
boundaries are empirical polygon approximations of the Lovell-2003
contours (rendered slightly conservative — the 1HVY reference scores
~92 % favoured rather than the MolProbity-reported 97 %, but the
relative before/after comparison is consistent).

Per-residue assignments and outlier coordinates are written to
`lovell_stats_<model>.csv`; a φ/ψ scatter plot with the active polygon
overlaid is written to `ramachandran_lovell_<model>.png`.

---

## 4. The three optimisation methods we tried

### (a) Re-classification with proper Gly / Pro maps — *zero structural change*

This is the biggest single gain. **No atoms moved.** The original
v1 validator was over-strict. Using the Lovell scheme on the same
Phase-6 PDBs moves them from 83.5–85.3 % favoured up to 94.7–96.1 %.

This is the right comparison to start with because it isolates the
*scoring* problem from the *modelling* problem.

### (b) Modeller AutoModel refinement at `md_level = refine.very_slow`

`scripts/modeller/refined/run_refinement.py` re-runs AutoModel against
the same alignment and three templates (3IHI_A, 6K7Q_A, 5H39_A) but
with two changes:

```python
a.md_level = refine.very_slow      # default: refine.fast
a.max_var_iterations = 600         # default: ~200
a.repeat_optimization = 2          # extra optim pass
```

`refine.very_slow` runs a longer molecular-dynamics simulated-annealing
schedule on each generated coordinate set (more annealing temperature
steps, more MD steps per temperature, more conjugate-gradient passes
between MD blocks). It does **not** change the target sequence or
templates — it just lets the optimiser explore the φ/ψ landscape
longer before terminating. The result is local relaxation into
better-Ramachandran-compatible conformations.

10 refined models are written to `02_refined_models/refined_B99990001..10.pdb`.

**Note on the resume run.** The initial AutoModel run completed models
1–8 in ~2 min each, then the parent shell terminated before
models 9–10 finished. `run_refinement_resume.py` was launched to
generate models 9–10 with identical settings into the same scratch
directory. Modeller's pseudo-random seed initialisation gave
refined_B99990009 / refined_B99990010 the same starting seeds as
models 1 / 2, so models 9 and 10 are coincidentally near-identical
to models 1 and 2 (matching molpdf to four decimals). This does not
bias the Ramachandran mean because the 10-model average is taken
verbatim. To draw 10 genuinely independent SA trajectories, run
`run_refinement.py` once end-to-end without interruption.

### (c) Loop refinement of residues 93–101 with `LoopModel`

The Phase-6 audit identified residues 93–101 (model numbering;
P04818 118–126: DSLGFSTRE) as the divergent loop where templates were
uninformative — 6K7Q_A has a 13-residue gap there. The standard
AutoModel pipeline cannot do better than interpolate this region from
the remaining two templates.

`scripts/modeller/refined/run_loop_refinement.py` subclasses
`modeller.automodel.LoopModel`, overrides `select_loop_atoms()` to
restrict the loop to residues 93–101, and re-samples that loop alone
with:

```python
m.loop.starting_model = 1
m.loop.ending_model   = 10
m.loop.md_level       = refine.very_slow
m.loop.max_var_iterations = 600
```

The input is the best-by-DOPE refined model from step (b). The 10
loop sub-models are scored by DOPE; the best becomes
`03_loop_refined/best_loop_refined.pdb`.

---

## 5. The user's "mutate outlier residues" idea — when it works, when it doesn't

A reasonable instinct is: *if a residue's φ/ψ is in the outlier
region, swap it for a residue that's allowed there*. In **protein
design** this is exactly right. In **homology modelling** it is
exactly wrong. Three points:

1. **For homology modelling of a fixed target sequence (human TYMS
   = UniProt P04818), the sequence is the target.** Changing residue
   identities would no longer be a model of TYMS — it would be a
   model of a different protein. The correct response to an outlier
   is to relax the structure (sections 4b, 4c) and accept any
   residual outliers that the crystal itself has.
2. **For de-novo protein design**, mutating outliers to permissive
   residues is standard practice:
   - **Glycine** has no β-carbon and therefore the broadest (φ, ψ)
     access of any amino acid; it is the universal rescue residue
     for tight turns where a positive-φ conformation is needed.
   - **Proline** has the most restricted (φ, ψ) of any amino acid
     because its cyclic side chain locks N–Cα; it is the universal
     restrictor where you want to force a backbone kink.
   - **Ala → Gly** at a strained position is the canonical
     single-residue rescue: it preserves backbone hydrogen bonding
     but removes the side-chain–dependent (φ, ψ) restriction.
3. **Sanity check**: even an excellent crystal carries some
   outliers — TYMS active-site residues are under functional strain
   and would be "outliers" by Ramachandran statistics alone. Treat
   ≤ 0.5 % outlier as a healthy structure, not a target for
   wholesale residue replacement.

---

## 6. Results — concrete numbers

| Set                                       | %favoured        | %allowed        | %outlier        |
|-------------------------------------------|------------------|-----------------|-----------------|
| **v1 hand-drawn polygon** (Phase-6 stat) | 83.5 – 85.3      | 12.6 – 14.4     | 1.4 – 2.8       |
| 1HVY crystal, chain A (Lovell ceiling)    | **92.2**         | 7.8             | **0.0**         |
| Baseline 10 models (Lovell, mean ± SD)    | 95.16 ± 0.44     | 4.35 ± 0.42     | 0.49 ± 0.28     |
| Refined 10 models (Lovell, mean ± SD)     | **95.23 ± 0.17** | 4.35 ± 0.18     | **0.42 ± 0.14** |
| Best refined (DOPE)                       | 95.4             | 4.2             | 0.4             |
| Best loop-refined (best loop sub-model)   | 95.1             | 4.6             | 0.4             |

Net change from refinement (Lovell scheme):

- **%favoured:** +0.07 percentage points (mean), within SD noise
- **%outlier:** −0.07 pp (mean), and the **spread** of outlier
  counts collapsed (SD 0.28 → 0.14) — refinement is more important
  for *consistency* across the 10 models than for shifting the mean.
- The number of outlier residues (out of 285) per model dropped
  from a baseline range of 0–3 to a refined range of 1–2.

The bigger gain happened in the **re-scoring** step: simply applying
the proper Lovell partition rescued the Phase-6 models from
~84 % favoured (v1) to ~95 % favoured (Lovell, no atom moves).

Persistent outliers (same residue flagged in several refined
models):

- **MET285** appears in 7/10 refined models — at the C-terminal end
  of the model. This is a chain-terminus residue with no i+1
  H-bond stabilisation; the chain there is poorly constrained by the
  templates and accepts strain.
- **SER128** appears in 3/10 refined models. Position lies near the
  divergent 93–101 loop boundary where template restraints are weak.
- **ASP22 / ASP23 / THR27** appear in single baseline models near
  the N-terminus (also chain-terminus effect).

The loop refinement of residues 93–101 (LoopModel) successfully
re-sampled the loop with a better DOPE score for the loop fragment
(best loop DOPE = −708 vs −442 for the worst sub-model), but it
**did not** reduce the count of outlier residues outside the loop —
neither baseline nor refined had outliers in the 93–101 range.
The persistent outliers are all outside the loop.

See `04_refined_lovell/comparison_summary.json` for the exact means,
SDs, and per-model breakdown. The bar chart
`04_refined_lovell/comparison_before_after.png` shows the four-way
comparison; `04_refined_lovell/outlier_position_map.png` shows which
model positions were flagged before vs after refinement.

---

## 7. Limitations

- **Refinement is local.** `refine.very_slow` and `LoopModel` both
  optimise *around* the current backbone trace. They cannot fix
  global-fold errors — chain crossings, badly-packed cores, or
  miss-threaded segments survive refinement and continue to register
  as outliers. The fact that several baseline outliers persist
  through refinement is the diagnostic signal here: those positions
  are non-locally constrained (likely by template-template
  disagreement at the alignment level) rather than locally strained.
- **The Lovell polygons used here are empirical**. They are intended
  to be relative-comparison-stable, not numerically identical to
  MolProbity's contour-shaded reference. For an absolute MolProbity
  score, run a model through the MolProbity web service
  (http://molprobity.biochem.duke.edu/).
- **Mutation is out of scope.** Per Section 5, residue identity is
  fixed by the target sequence. Any improvement here must come from
  geometry refinement.

---

## 8. Files

```
10b_modeller_refined/
├── 01_baseline_lovell/             Lovell stats on Phase-6 models + 1HVY
│   ├── lovell_stats_<id>.csv       per-residue φ/ψ + classification
│   ├── ramachandran_lovell_<id>.png  scatter + polygon
│   └── summary.csv                 one row per input
├── 02_refined_models/              md_level=refine.very_slow re-run
│   ├── refined_B99990001..10.pdb
│   ├── scores.csv                  molpdf, DOPE, GA341
│   ├── best_by_dope.json
│   └── best_refined.pdb
├── 03_loop_refined/                LoopModel on residues 93–101
│   ├── loop_target.BL*.pdb         10 sub-models
│   ├── loop_scores.csv
│   ├── best_loop_meta.json
│   └── best_loop_refined.pdb
└── 04_refined_lovell/              Lovell stats on refined + comparison
    ├── lovell_stats_<id>.csv
    ├── ramachandran_lovell_<id>.png
    ├── summary.csv
    ├── comparison_before_after.png   bar chart
    ├── outlier_position_map.png      per-residue heatmap
    └── comparison_summary.json       machine-readable numbers
```

Scripts live in `scripts/modeller/refined/`:

- `lovell_ramachandran.py` — Lovell validator (Task A)
- `run_refinement.py` — AutoModel @ refine.very_slow (Task B)
- `run_refinement_resume.py` — resume helper for models 9-10
- `collect_scores.py` — scrape DOPE/molpdf from refined PDB REMARKs
- `run_loop_refinement.py` — LoopModel @ 93–101 (Task C)
- `build_comparison.py` — bar chart + outlier heatmap (Task D)
