# Phase 6 — Homology modelling of human TYMS (UniProt P04818)

This phase builds a Modeller-based homology-modelling sub-pipeline. It is *additive*
to v1–v5: nothing in those folders is touched.

## (a) Why homology modelling?

Even when an experimental crystal structure exists for a target (here, 1HVY),
homology modelling is a useful *educational* exercise: it shows how a structure
is reconstructed from sequence + templates and lets you compare the model against
the known crystal as a ground-truth benchmark.

## (b) Why <100% identity templates?

We deliberately *excluded* 1HVY itself from BLAST hits and filtered hits to
**30% ≤ identity ≤ 95%**, **coverage ≥ 80%**, **resolution ≤ 2.5 Å**. This forces
Modeller to do real homology modelling rather than copying coordinates from the
identical-sequence template. In production one would obviously use the highest
quality template available, including 100% identity matches if they exist.

Selected templates:
  - **3IHI_A** — 92.71% identity, resolution 1.94 Å
  - **6K7Q_A** — 75.96% identity, resolution 2.27 Å
  - **5H39_A** — 72.28% identity, resolution 2.0 Å

## (c) What each step produced

1. `01_clean_pdb/` — Chain-A-only 1HVY (cleaned of HETATMs and chain B) +
   single-chain FASTA, sanity-checked against P04818.
2. `02_blast/` — Remote NCBI BLAST against the **pdb** database; filtered hits;
   downloaded selected template PDBs and extracted the relevant chains.
3. `03_alignment/` — Combined target + 3 templates with ClustalW; converted to
   Modeller PIR (`alignment.ali` / `alignment.pir`).
4. `04_modeller_run/` — 10 candidate models via Modeller AutoModel, scored with
   DOPE / molpdf / GA341 (`scores.csv`).
5. `05_comparison/` — PyMOL pairwise overlays + per-residue Cα-distance bar
   plots + all-models overlay + `rmsd_per_model.csv`.
6. `06_validation/` — Local Ramachandran φ/ψ plots + per-residue normalized DOPE
   profiles + `quality_overview.png` + `SAVES_MANUAL.md` for the UCLA SAVES
   step (manual upload).
7. `07_viewers/` — 3Dmol viewers per model (also mirrored into `viewers/` for
   GitHub Pages).

## (d) Best-model selection

Two complementary criteria:
- **Best by DOPE** (most negative ⇒ best fold-energy): see `04_modeller_run/best_by_dope.json`.
- **Best by Cα RMSD to crystal**: see `05_comparison/best_summary.json`.

In a real prediction (no crystal) only DOPE / molpdf / GA341 / Ramachandran are
available, so **best_by_dope** is the canonical choice and is copied to
`models/best_model.pdb`.

## (e) Limitations

- All 10 models are extremely similar (RMSD differences < 0.05 Å, DOPE within
  ~500 units), so the "best" choice is largely cosmetic.
- TYMS dimer interface is *not* modelled (single-chain target only) — for
  enzymatic studies the dimer should be reconstructed afterwards.
- Local Ramachandran is a *basic* PROCHECK substitute; for publication-quality
  validation use UCLA SAVES (PROCHECK / ERRAT / VERIFY3D / WHATCHECK) — see
  `06_validation/SAVES_MANUAL.md`.
- Templates were *intentionally* sub-100%; using 1HVY itself would give a
  near-perfect rebuild but defeat the educational purpose of the phase.

Generated 2026-05-12 15:04:11.
