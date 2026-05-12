# V3 Scientific Officer Review

Agent role: Scientific officer (round 3).
Date: 2026-05-12.

## Verdict
**Ship with caveats.** v3 closes the four highest-priority v2 findings (sign convention documented, WT holo dock recovered to a credible -8.314 / 2.08 Å, terminology cleaned to "Δ Vina score", explicit Limitations section added). However, two scientifically important problems persist:
(i) **WT holo dock returns only 3 poses**, and the **C195A holo Δ = -2.29** and **C195A_H196A Δ = -1.84** strongly suggest the WT holo run is bottoming out on a poorly-sampled local minimum that mutants escape — i.e. the WT reference is still soft.
(ii) The top-destabiliser ranking is **chemically reasonable but not strong**: best Δ Vina is +0.73 kcal/mol, well within Vina's empirical noise (~1 kcal/mol).

The pipeline can ship as a *qualitative* ranking exercise; it cannot ship as a quantitative ΔΔG study, and the report must say so more loudly.

## Items v3 closed
1. Sign convention documented; T170A control near zero.
2. Terminology cleaned (Δ Vina score, not ΔΔG).
3. Limitations section present, covers Vina ≠ ΔΔG, polyglutamylation, apo as negative control, T170A control behaviour, rigid receptor/ligand, RMSD definition.
4. Apo column explicitly framed as negative-control.
5. T170A behaves as null control.
6. Mis-dock filter applied (RMSD > 3 Å excluded from rankings).

## Items v3 did NOT close
1. **WT holo pose count** = 3 — top affinity gap (-8.3 → -6.1) means the WT reference is on a narrow funnel. C195A and C195A_H196A recover *more* poses with *better* affinities — same pathology as v2, less severe.
2. **C195A holo Δ = -2.29** is **not** scientifically defensible. C195 is the catalytic nucleophile that forms the covalent Michael adduct with dUMP; ablating its thiol cannot increase non-covalent affinity by 2.3 kcal/mol. Real explanation: WT holo undersampled and C195A frees up sampling. Must be flagged.
3. **C195A_H196A Δ = -1.84 holo** has the same problem.

## New issues introduced
1. **apo–holo r = 0.80 is misleading.** Dominated by mutants clustering near Δ ≈ 0 ± 0.5 in both conditions. Both axes are noise-dominated.
2. Holo top-5 ranking magnitudes (≤0.73 kcal/mol) are below Vina's discrimination threshold.
3. n_modes = 2–4 for many holo mutants — comparing top affinities across runs with such different sampling depths is statistically suspect.
4. R215A_N226A apo is mis-docked (RMSD 5.34 Å, Δ=+1.08) but holo is well-docked — interpret pair-wise effects with caution.

## Required additions before final report
1. **Re-dock WT holo with seed sweep ≥ 5** until n_modes ≥ 10, then recompute all `delta_vina_vs_wt`.
2. **Add explicit caveat next to the C195A row**: "negative Δ for the catalytic nucleophile is inconsistent with TYMS mechanism; attribute to under-sampled WT reference, not to true tighter binding."
3. **Replace Pearson r with Spearman rank correlation restricted to |Δ| > 0.3 kcal/mol**.
4. **Quote Vina noise floor (~0.7–1.0 kcal/mol)** in Limitations and state no mutant exceeds it on holo.
5. **Drop or footnote holo rows with n_modes < 5** from the destabiliser table.
6. **Polyglutamylation scope**: state that only mono-glutamate raltitrexed is probed.
