# v2 Scientific Officer Review

Agent role: Scientific Officer — round 2 (peer-review-grade defensibility).
Date: 2026-05-12.

## Verdict
**Ship with caveats — but only after the sign / labelling problem is fixed.** v2 closes the structural and panel-construction concerns from v1 cleanly: real TYMS orthologs across 5+ kingdoms, gap-corrected JSD, full A+B dimer, no force-augmentation, G217W properly skipped. However the docking analysis itself contains an internal inconsistency that must be addressed: nearly every "destabilising" mutation has a NEGATIVE ΔΔG, which contradicts the report's own legend ("Positive = destabilising"). The root cause is a near-failed WT holo dock (`top_affinity = -3.215` with only 3 poses, `mean_top3 = -1.69`, RMSD 4.32 Å) — the WT reference is itself an outlier, so almost any mutant looks "better." This needs framing, not erasure, before the report ships.

## Items v2 successfully closed
- **Real ortholog panel.** `01b_msa_v2/input.fa` lists 10 bona-fide TYMS accessions (P04818 human TYMS, P0A884 *E. coli* ThyA, P00469 *L. casei* TYMS, P00471 T4 phage td, etc.) covering Metazoa (4), Fungi, Plantae, Bacteria (2), Protozoa, and Bacteriophage. No DHFR/polyubiquitin contamination. C. elegans rejected for length (1059 aa).
- **Conservation now meaningful.** `01b_msa_v2/conservation_scores.csv`: C195 100%, H196 96.1%, R175 89.5%, R176 90.8%, R215 94.4%, N226 97.1% — all in or at the edge of the top decile. Long-runner claim verified.
- **Homodimer handled.** `protein_dimer_h.pdb` contains 4626 ATOM records each in chains A and B (symmetric).
- **Apo vs holo behave differently.** Apo poses cluster at top_affinity ≈ −7.0 to −7.9 with RMSD 5–8 Å (mis-docked, insensitive). Holo poses range −4.4 to −8.3 with RMSD ~2.1–2.2 Å (in-pocket, discriminating).
- **G217W dropped.** No G217W rows in `mutant_results_v2.csv`.
- **No augmentation.** `selected_meta.json: force_augmented: []`.
- **Arg-clamp panel ran.** R50A, R175A, R176A, R215A all present.

## Items v2 did NOT close
- **Sign of ΔΔG is backwards relative to the legend (HIGH severity).** Report Figure 6 says "Positive = destabilising" but the "Top destabilising mutations (holo)" table lists H196A at ΔΔG = −1.17, R50A at −2.79, etc. Either the WT reference is broken or the sign convention is flipped.
- **Δ Vina vs ΔΔG terminology (MEDIUM).** Report calls the Vina-affinity differences "ΔΔG" throughout. Vina is an empirical scoring function, not a free energy.
- **Polyglutamylation scope limitation (MEDIUM).** Not mentioned in v2 report.
- **Apo "insensitivity" is partly mis-docking (MEDIUM).** The apo column is largely uninformative.
- **WT holo dock is unreliable (HIGH).** Only 3 poses survived for WT holo with `mean_top3 = -1.69` and RMSD 4.3 Å. The whole holo-ΔΔG analysis rides on this single weak number.

## New issues introduced by v2
- `mean_top3` is `nan` for 11/20 holo mutants — analysis silently uses `top_affinity` only.
- **`T170A` ("control_surface") shows holo ΔΔG = −1.27**, similar magnitude to H196A (−1.17). A surface control behaving like a catalytic-dyad mutant is a red flag for the methodology, not a positive control passing.

## Required additions before final report
1. **Re-dock WT holo** with larger exhaustiveness / different seed until ≥10 poses survive and RMSD < 3 Å, then recompute all holo ΔΔGs against this corrected reference.
2. **Fix the sign/legend mismatch.** Either flip the convention so destabilising = positive ΔΔG, or relabel figures and tables.
3. **Add a Limitations section** covering: (a) cofactor polyglutamylation not modelled, (b) Vina score is not a true ΔΔG of binding, (c) apo runs are dominated by mis-docking and should be treated as a negative control, (d) T170A surface-control behaving like an active-site mutant indicates residual scoring noise.
4. **Annotate mis-docked rows** (RMSD > 3 Å) in the mutant tables; exclude them from "top destabilisers" rankings or move to a separate "pose lost" table.
5. **Restore mean_top3 for all mutants** (fall back to mean of all surviving poses when fewer than 3) so the column is not silently NaN.
