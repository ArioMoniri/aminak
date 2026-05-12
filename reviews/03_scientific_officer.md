# Scientific Officer Review

Agent role: Strategic / interpretive review (peer-reviewer perspective).
Date: 2026-05-12.

## Verdict

**Needs major revision before any external dissemination.** The structural/docking half (target choice, WT redocking RMSD 1.08 Å, mutant panel design, control) is competently executed and the limitations section is unusually candid. However, **the cross-species conservation panel is broken** — most of the eight "orthologs" are not TYMS proteins at all, which invalidates the JSD scores, the "top-25%" selection, and the framing of the entire analysis as a *conservation*-driven study. Until the MSA is rebuilt, the docking results stand on their own merits but the central claim ("conserved active-site residues drive binding") is unsupported.

## Strengths

- Defensible target/CRC framing: TYMS is the validated 5-FU target; 1HVY (1.9 Å, dUMP + raltitrexed in chain A) is a reasonable starting structure.
- The natural substrate dUMP was actually docked (not an analog); WT re-dock RMSD of 1.08 Å with top affinity −7.73 kcal/mol is a credible positive control.
- Mutational panel logic is thoughtful: Ala-scan + chemically-opposite singles + mechanism-motivated doubles (catalytic dyad, Arg175E/R176E phosphate-clamp swap, charge swap, aromatic swap, polar-neutral compensator).
- Distant-surface negative control **T170A (Δ = +0.04 kcal/mol)** is appropriate — chemically inert, ~18 Å away, behaves as expected (∼0 effect), validating the box/protocol.
- The C195 anomaly is *acknowledged* in the analysis: rigid-receptor docking of substrate cannot capture loss of the Cys195 nucleophilic attack.

## Weaknesses (severity-ordered)

1. **CRITICAL — the ortholog panel is wrong.** Of the 9 UniProt accessions used in `scripts/stage1_msa.py`:
   - `P0CG53` (S. cerevisiae) = **polyubiquitin**
   - `P11849` (T4 phage) = photosystem II protein D2 family, not phage TS
   - `P04996` (L. casei) = ~360-aa membrane protein (correct: `P00469`)
   - `P04394` (E. coli) = NOT E. coli ThyA (correct: `P0A884`)
   - `O44019` (C. elegans) and `P21520` (P. falciparum) — pairwise identity vs human is **15–20%** with only 80–260 aligned columns out of 313.
   Only Homo + Mus are real TYMS. JSD scores are noise. This is exactly why the pipeline had to *force-augment* the top-25% set with Cys195/His196.

2. **MAJOR — phylogenetic spread is borderline even with correct IDs.** Eight sequences with no archaeal representative. A defensible TYMS conservation panel needs ≥15 sequences spanning archaea/bacteria/fungi/plants/metazoa with the 30–80% identity sweet spot enforced.

3. **MAJOR — "force-augment Cys195/His196 into the top-25%" is not honestly disclosed in the report.** With a correct alignment these residues *are* invariant; the augmentation became necessary only because of weakness #1. As-written this looks like cherry-picking.

4. **MAJOR — homodimer not handled.** TS is an obligate homodimer with a composite active site spanning chains A and B. Stage-3 strips to chain A only, deletes the partner subunit, and docks into the resulting half-pocket. This invalidates contributions from across-the-interface contacts (Arg175/176) and biases all Δ-affinity values.

5. **MODERATE — raltitrexed/D16 was removed but the cofactor pocket is left empty.** Methylene-THF cofactor pocket physically constrains the dUMP site; docking dUMP into an apo (cofactor-free) pocket overestimates pocket flexibility and can artificially improve scores — exactly the "C195S/C195A score better than WT" pattern observed.

6. **MODERATE — Arg phosphate clamp coverage is thin.** R175E_R176E is included as a double, but no single-mutant breakdown (R175A, R176A, R215A) — yet Arg215 is annotated as a binding-site residue and was *not* mutated despite being in the conserved set.

7. **MINOR — overstatement risk.** "Δ affinity" is repeatedly used as if it measured ΔΔG of binding. Vina score is a heuristic; phrase as "Vina score change" throughout.

## Required additions before final report

1. **Rebuild the MSA with verified TYMS sequences** (correct UniProt IDs: P00469 L. casei, P0A884 E. coli ThyA, P07807 S. cerevisiae CDC21, Q23381 C. elegans, Q08758 P. falciparum, P04019 T4 td, plus Arabidopsis Q05762, an archaeon, etc.); ≥12–15 sequences, document % identity matrix, target 30–80% range. Recompute JSD and re-run the active-site overlap. Drop the manual augmentation step.
2. **Re-prepare the receptor as the biological homodimer** (chains A+B); re-run WT redocking and the full mutant panel.
3. **Run with the methylene-THF cofactor (or raltitrexed) retained** as a second docking condition; compare Δ-affinity rank order.
4. **Add R175A, R176A, R215A, R50A single Ala-scans** so the Arg phosphate-clamp claim is actually probed.
5. **Per-pose clustering and pose-energy analysis** (cluster the nine Vina modes per mutant; report cluster populations and whether the top score corresponds to the productive pose).

## Optional improvements

1. MM/GBSA or single-point rescoring of top poses.
2. Brief Rosetta/FoldX ΔΔG_fold for each mutant — separates "destabilises pocket" from "removes a contact."
3. Explicit paragraph on **polyglutamylation** of the physiological cofactor.
4. Reframe the headline from "conserved residues drive binding" to "active-site annotated residues, intersected with conservation, drive binding."
5. Replace "Δ affinity" with "ΔVina score (kcal/mol-equivalent)" throughout.
