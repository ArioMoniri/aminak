# Phase-6 Scientific Officer Review

Round 6. Date: 2026-05-12.

## Verdict: PUBLISH-READY (educational tier) — with two recommended additions

The Phase 6 sub-pipeline is internally consistent, scientifically defensible at its declared educational scope, and the numbers check out. Two narrative gaps prevent a flawless rating.

## Strengths
1. **Templates are real TYMS, organism-diverse**: 3IHI = mouse TYMS, 6K7Q = shrimp TYMS, 5H39 = KSHV ORF70 thymidylate synthase. Three independent eukaryotic/viral lineages spanning 92→72% identity.
2. **Best-by-DOPE selection is correct for blind prediction.** In the absence of a crystal one cannot know RMSD; DOPE/molpdf/GA341 is the only legitimate criterion.
3. **DOPE/molpdf magnitudes sensible**: DOPE ≈ −35,400 to −35,775 for a 313-aa enzyme (~−113/residue, in the well-folded globular band). GA341 = 1.0 across all 10 models.
4. **0.37 Å Cα RMSD to 1HVY is excellent** for templates 70-90% identity.
5. **Educational caveat clear** in both README and DOCX.
6. **SAVES manual instructions adequate**.

## Weaknesses
1. **Ramachandran statistics need a reference frame.** Models 83.5-85.3% favoured / 12.6-14.4% allowed / 1.4-2.8% outlier. MolProbity norm is ≥90% favoured. **However**, the 1HVY crystal itself under the same hand-drawn polygon scheme scores 82.3% favoured / 15.6% allowed / 2.1% outlier — i.e. **models match or beat the crystal under this scheme**. The polygons are conservative, not the models.
2. **No bridge to Phase 5 docking finding**: would re-docking against a Modeller receptor change the v5 conclusion?
3. **AutoModel run at default refinement** (no `library_schedule`, no `md_level`).
4. **10-model ensemble is small** (Modeller best-practice 20-100).
5. **Dimer interface missing**: TYMS is obligate homodimer; monomeric model is mechanistically incomplete.
6. **PSI-domain / cofactor sensitivity**: TYMS undergoes apo↔ternary loop conformational change.

## Required additions before final report
1. **Ramachandran calibration paragraph**: state that the local validator's hand-drawn polygons are more conservative than MolProbity; the 1HVY crystal scores 82.3% favoured under the same scheme.
2. **Docking bridge paragraph**: substituting a homology receptor would add template-biased rotamer noise; v5 conclusion would not change.
3. **Limitations augmentation**: AutoModel was at default refinement; production work needs `md_level=refine.slow` + 20+ models. TYMS undergoes cofactor/substrate-induced loop change; model inherits an undefined intermediate state.
