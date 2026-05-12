# CHANGELOG

Multi-agent doer↔verifier audit history. Every iteration here was driven by a critical issue surfaced by 4 specialised review agents (validator, code reviewer, scientific officer, structural bioinformatician) running in parallel against the previous version. New agents reading this file should consult the round-N report files in `reviews_vN/` for the verbatim verdicts.

## Phase 6 — Modeller homology modelling (round 6)

**Status:** under audit. 4 reviewers running.

**Built:** 10 homology models of TYMS using Modeller 10.8 against 3 templates with 30–95 % sequence identity (educational, not 100 %), validated locally with Ramachandran φ/ψ + DOPE per-residue + RMSD vs crystal, plus a manual SAVES-upload procedure.

- Templates: 3IHI (92.7 %), 6K7Q (76.0 %), 5H39 (72.3 %).
- Best by Cα RMSD: model 10 (0.367 Å vs 1HVY chain A).
- Best by DOPE: model 3 (DOPE = −35775.75).
- Ramachandran (local): favoured 83.5–85.3 %, allowed 12.6–14.4 %, outliers 1.4–2.8 % across the 10 models.
- SAVES web upload (PROCHECK, ERRAT, VERIFY3D, WHATCHECK) is documented as a manual step in `10_modeller/06_validation/SAVES_MANUAL.md` because SAVES has no programmatic API.

## Round 5 — v5 (FINAL for the docking phase)

**Verdict (4 reviewers):** Validator 10/10 PASS · Code review v4 punch list closed · Scientific officer SHIP · Structural bioinformatician CONDITIONAL PASS (3 reporting-only items, all baked into `report_FINAL.docx`).

**Hard fix in this round:** in-place reprotonation of the original 1HVY crystal cofactor coordinates (0.000 Å heavy-atom drift, 0 protein clashes). v4 had Kabsch-aligned the CCD-ideal D16 onto the bound conformer, producing 2.71 Å heavy-atom drift plus a 1.95 Å clash to PHE 80 CD2 — that placement artefact drove the v4 "cofactor expels dUMP" interpretation, which v5 dissolved. WT holo recovered to top affinity −8.25 kcal/mol with named-RMSD 0.999 Å vs crystal dUMP.

**Headline:** Rigid-receptor AutoDock Vina with AD4 partial charges and the physically correct (net −2) raltitrexed cofactor cannot resolve TYMS active-site point mutants at the kcal/mol scale. Largest holo Δ Vina = +0.77 kcal/mol (R215A_N226A) — well below Vina's ±0.85 noise floor. Null-result methodology paper.

## Round 4 — v4

**Critical issue surfaced:** v3's cofactor "pH 7.4 reprotonation" was a no-op (output file byte-identical to v2). Holo electrostatics still wrong → C195A negative Δ was an artefact, not biology.

**Fix in v4:** RDKit reprotonation from the PDB Chemical Component Dictionary's ideal SDF, bond-order assignment, deprotonation of both carboxylates, MD5-inequality assertion, plus atom-name preservation through the obabel↔Vina round-trip. v4 closed the v3 punch list — and exposed the deeper placement-artefact issue that round 5 fixed.

## Round 3 — v3

**Critical issues surfaced:**
- Receptor PDBQT had all-zero partial charges (silent meeko CLI fallback writing no charges).
- WT holo dock barely converged (3 poses, RMSD 4.32 Å).
- PyMOL `get_strain()` returned 1e30 in headless mode, so v2's rotamer minimisation was a no-op (every mutant used rotamer 0).
- Sign convention was backwards.

**Fix in v3:** charge waterfall (obabel Gasteiger → meeko → pdb2pqr) with `max|q| > 0.05` gate, multi-seed WT holo with affinity-based selection, sculpt rotamer minimisation, positive=destabilising sign convention with `mean_topk = mean(top min(3,n))`.

## Round 2 — v2

**Critical issues surfaced (the v1 audit):**
- 5 of 9 ortholog UniProt accessions point to wrong proteins (`P0CG53` is yeast polyubiquitin, `P11849` is a T4-phage photosystem-II-family protein, `P04996` is an L. casei membrane protein, etc.). Only human + mouse were real TYMS. JSD scores were noise.
- TYMS is an obligate homodimer; chain B was discarded though catalytic Arg175′/Arg176′ come from the partner subunit.
- The G217W mutant had a 0.98 Å Trp clash (sterically impossible).
- CME43 (covalently-modified Cys43) was silently dropped, leaving a 42→44 backbone gap.

**Fix in v2:** real TYMS panel of 10 verified orthologs across 5+ kingdoms; A+B dimer kept; CME43 re-mutated to native CYS in place; G217W dropped on the heavy-atom clash check; both apo and holo dockings produced.

## Round 1 — v1 (initial)

Initial pipeline. Stages 1–9 ran end-to-end against human TYMS (P04818) + dUMP from PDB 1HVY. The 4 reviewers immediately surfaced the critical MSA + dimer + G217W + CME43 issues that v2 had to fix.

## Reviewer reports

| Round | Folder |
| --- | --- |
| 1 (v1) | [`reviews/`](reviews/) |
| 2 (v2) | [`reviews_v2/`](reviews_v2/) |
| 3 (v3) | [`reviews_v3/`](reviews_v3/) |
| 4 (v4) | [`reviews_v4/`](reviews_v4/) |
| 5 (v5) | [`reviews_v5/`](reviews_v5/) |
| 6 (Phase 6, in progress) | [`reviews_phase6/`](reviews_phase6/) (will be added when reviewers complete) |
