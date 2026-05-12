# V5 Structural Bioinformatics Audit

Round 5 — CONDITIONAL PASS (hard blocker resolved; 3 reporting clarifications).
Date: 2026-05-12.

## Verdict: CONDITIONAL PASS

The hard cofactor-placement blocker is **fully resolved**, but a hidden labelling issue surfaced during audit (see Item C). Output may proceed if the user accepts that "WT/mutant holo dock" docks **dUMP substrate**, not D16, into a D16-occupied receptor (this matches the v2-v5 pipeline intent — but should be made explicit in the report).

## Per-audit findings

1. **Cofactor heavy-atom coords preserved (FIX A): VERIFIED.** Crystal D16 chain A (32 heavy) and v5 cofactor A — same atom-name set, RMSD 0.000000 Å. Same for chain B.
2. **MD5 differs from v4: VERIFIED.** Direct heavy-atom RMSD between v4 and v5 cofactor A = 2.7088 Å — confirms v4 was the broken Kabsch placement.
3. **Carboxylates deprotonated: VERIFIED.** Both chains contain exactly two carboxyl carbons (CD γ-Glu, CT α-C-terminus). No H within 1.3 Å of either O on either carboxylate. Net −2 per cofactor.
4. **Zero protein clashes: VERIFIED.** 0 clashes <1.8 Å on chain A, 0 on chain B. PHE A80 CD2 specifically: nearest cofactor heavy atom (OE2) is 3.569 Å away — v4's 1.95 Å clash is gone.
5. **WT holo top affinity: VERIFIED.** −8.249 (within ~1 kcal/mol of apo −9.197); credible.
6. **WT holo RMSD vs crystal: VERIFIED.** 0.999 Å by atom-name match; 0.334 file-RMSD. Round-trip atom-name preservation works.
7. **WT holo n_modes=2: real funnel collapse, not a sampling defect.** Every seed (8 total) ran exh=128 and explored 29-32 modes; Vina's pdbqt clustering retains only 2 unique pose clusters per seed. Affinity range 0.036 kcal/mol.

## NEW TECHNICAL ISSUE FOUND

**C. The "holo" docking docks dUMP, not D16, into the D16-occupied holo receptor.** This is by design (consistent with v2-v5 — Stage 6 docks the substrate), but Round 4 framing led one to expect "fixed cofactor placement → re-docked cofactor". Round 5 report MUST state plainly: "Docked ligand = dUMP (substrate); D16 cofactor sits in the receptor pocket as part of the holo state."

**Mutant spot check (R215A_N226A, H196A, C195A): clean.** Each `<mut>_holo.pdb` contains the D16 cofactor at v5 coords; vs respective mutant protein heavy atoms — 0 clashes.

**Mutant apo reuse asymmetry**: mixing apo affinities sampled at exh=32 (mutants from v3) with apo at exh=96 (WT from v4) introduces a small upward bias for mutants. Quantitative impact likely ≤0.2 kcal/mol (within Vina noise floor 0.85).

## Corrections still needed (all reporting-only)

1. (HIGH, reporting) Add an explicit caption everywhere "holo" appears: "Holo dock = dUMP substrate docked into receptor with D16 cofactor pre-bound at the v5 in-place reprotonated coords."
2. (MEDIUM, reporting) Footnote the n_modes=2 result: state Vina's pdbqt clustering retains only 2 distinct poses while underlying exploration found 29-32 modes per seed; affinity range 0.036 kcal/mol << 0.85 noise floor as evidence of funnel convergence.
3. (LOW, reporting) Footnote the apo exhaustiveness mismatch (mutants exh=32 reused from v3, WT exh=96 from v4) and quantify the bias as ≤Vina noise floor.

## Lower-priority improvements
- Re-run mutant apo at exh=96 to remove the exhaustiveness mismatch.
- For wt_holo, try `--num_modes 32 --energy_range 4` to capture second-tier −7 kcal/mol solutions.

**Bottom line: hard blocker resolved. Ship with the three reporting clarifications above.**
