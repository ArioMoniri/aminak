# V5 Scientific Officer Review

Round 5 — SHIP as null-result methodology paper.
Date: 2026-05-12.

## Verdict
**(a) v5 SHIPS as a null-result methodology paper.** All five round-4 prose fixes are baked in, the cofactor placement artefact is fully resolved (independently verified: 0.0000 Å heavy-atom displacement vs 1HVY chain-A D16; zero protein-cofactor pairs <2.0 Å), the WT-holo top pose now overlays the crystal dUMP at 0.999 Å named-RMSD, and the headline 0/20-mutants-above-noise-floor result is preserved. One residual labelling weakness (see below) needs a one-line caption fix; otherwise this is publishable.

## Items v5 closed
1. **Cofactor placement artefact** — fully resolved. v5 receptor `protein_dimer_holo.pdb` retains the crystal D16 coordinates exactly. PHE-80 1.95 Å clash from v4 is gone (nearest is now 3.569 Å).
2. **WT holo top affinity = −8.249 kcal/mol**, tightly clustered across 8 seeds (range 0.036 kcal/mol).
3. **WT holo RMSD vs crystal**: file-RMSD 0.334 Å, named-RMSD 0.999 Å vs chain-A dUMP. Native re-docking. v4's 12.95 Å is dead.
4. **All 5 sci-off prose fixes baked in**:
   - (a) Headline sentence verbatim;
   - (b) Mechanistic paragraph naming the 2.71 Å v4 RMSD and the cofactor-A O1 ↔ PHE 80 CD2 1.95 Å clash;
   - (c) "No statistically significant apo–holo correlation" wording;
   - (d) C195A row visually flagged (pink `.row-c195a` class);
   - (e) Dual-RMSD-reference note in Methods §7a.
5. **Headline survives**: `n_above_noise_holo = 0`; top holo Δ = +0.772 (R215A_N226A), well below ±0.85.

## Items v5 did NOT close
1. **C195A holo Δ went MORE negative, not less** (−2.25 vs v4's −1.75 vs v3's −2.29). Same n_modes=2 sampling pathology as WT holo. Marked `low_confidence=True`, excluded from rankings, pink-flagged. Defensible only because the headline owns it.
2. **WT-holo `low_confidence=True` itself** (n_modes=2). The WT reference against which every mutant Δ is computed is itself flagged low-confidence under v5's own n_modes<5 rule.

## New issues
1. **"Funnel collapse" explanation for WT-holo n_modes=2 is partially defensible**. Strong evidence: 8 seeds converge to top affinity within 0.036 kcal/mol; RMSD 0.32–0.34 Å across seeds; named-RMSD 0.999 Å. Weak evidence: same n_modes≤2 hits 8 of 20 holo mutants, suggesting a search-space / clustering-radius artefact. Acceptable as working hypothesis if hedged.
2. **Apo–holo correlation values** (Pearson r=−0.740, p=0.470; Spearman ρ=−0.500, p=0.667; both n=3) reported with explicit non-significance.
3. **Mutant-apo / WT-apo protocol asymmetry** disclosed in Methods §7b. Quantitative impact ≤Vina noise floor.

## Required additions before final report
1. **One-line caption fix on `top5_destab_holo_clean`**: WT holo itself is `low_confidence=True` under the same rule that excludes 14/20 mutants. Add: *"Δ values are computed against a WT-holo reference that itself sits at n_modes=2; this is a property of the binding-site funnel under rigid-receptor Vina, not a confidence weakening of mutant Δ values."*
2. **Hedge the "funnel collapse" sentence**: replace "consistent with high-confidence binding, not poor sampling" with "consistent with convergent sampling to the crystallographic pose; we cannot exclude that Vina's clustering radius truncates additional near-native modes."
3. **R215A_N226A apo is mis_docked** (RMSD 5.34 Å) — confirm pink-flagged in apo top-5.

The top-3 holo destabilisers are chemically plausible (R215 contacts dUMP phosphate; H196 is conserved catalytic dyad partner; the double removes both an orienting Arg and an H-bond donor) AND all sit below ±0.85 noise floor. **Right shape for a null-result methodology paper. Ship it with the three additions above.**
