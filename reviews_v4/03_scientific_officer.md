# V4 Scientific Officer Review

Round 4. Date: 2026-05-12.

## Verdict
v4 is methodologically honest and closes most v3 deficiencies. The headline finding — **no mutant signal exceeds Vina's noise floor under the (correctly ionised) holo condition** — is defensible and should be reported as the principal scientific result. The 12.95 Å WT-holo RMSD is real and explainable. Recommendation: **publish as a null-result methodology paper**.

## Items v4 closed
1. WT holo n_modes = 30 (spec ≥ 10).
2. Cofactor reprotonation real (MD5 differs; both carboxylates [O-]).
3. Holo signal: 0 mutants exceed Vina ±0.85 noise floor (top destab R176E +0.204).
4. Polyglutamylation statement present verbatim.
5. Vina noise floor citation present.
6. Atom-name preservation: WT apo `rmsd_top_to_native` = 0.9116 Å exactly matches `rmsd_top_named` = 0.9116 Å.
7. Mis-dock annotation: holo redefined as `|RMSD − WT_holo_RMSD| > 3 Å`.

## Items v4 did NOT close
1. **C195A holo Δ direction** — moved from −2.29 → −1.75. Still negative. Acceptable because filtered out of headline, but the residual −1.75 should not appear without `mis_docked` flag adjacent.
2. **Apo–holo correlation interpretation** — Pearson r = −0.30 (p=0.40), Spearman ρ = −0.45 (p=0.19), filtered ρ = +0.60 (n=4). Sign flip not discussed; with n=4 the +0.60 is meaningless.
3. **Mechanistic explanation of WT-holo 12.95 Å** — report attributes to "narrow funnel" but never explicitly states the cofactor-expulsion hypothesis.

## New issues introduced
1. Reference-frame ambiguity: holo "RMSD vs crystal" is now meaningless because WT itself is at 12.95 Å.
2. Selection-rule change: WT-holo "best" pose may be local minimum; affinity tie-breaker by n_modes does not address this.
3. Spearman flip not flagged — invites cherry-picking accusations.

## Required additions before final report
1. **One paragraph in §1** explaining why WT-holo RMSD is 12.95 Å (cofactor electrostatics expel dUMP from canonical pocket into remote sub-pocket).
2. **Headline sentence**: "Rigid-receptor AutoDock Vina with AD4 partial charges and the physically correct (net −2) raltitrexed cofactor cannot resolve TYMS active-site point mutants at the kcal/mol scale."
3. **Drop or annotate** the unfiltered apo–holo Spearman/Pearson; remove filtered n=4 ρ from prose.
4. **Suppress the C195A −1.75** from any summary surface that does not also display `mis_docked`.
5. **Flag the dual-RMSD-reference issue** in methods.
