# Phase-6 Structural Bioinformatics Audit

Round 6. Date: 2026-05-12.

## Verdict: CONDITIONAL PASS

Pipeline is technically sound and reproducible, but contains:
1. one silent sequence bug (Cys43 missing from local FASTA);
2. a methodology mis-statement on RMSD (Cα-refined, presented as raw);
3. an under-honest framing of the exercise (3IHI dominance not stated).

## Per-audit findings

1. **Sequence-identity check vs P04818** — PARTIAL FAIL. Local FASTA: 287 aa; UniProt: 313 aa. **Cys43 is missing entirely** (1-residue deletion). The catalytic Cys becomes "Cys198" in the model coordinate system, not Cys195. This is cosmetic for modelling but will confuse anyone cross-referencing the literature.

2. **BLAST templates** — PASS. All three confirmed real TYMS:
   - 3IHI = TYMS, *Mus musculus*, 1.94 Å.
   - 6K7Q = TYMS, *Penaeus vannamei* (white shrimp), 2.27 Å.
   - 5H39 = ORF70 = TYMS of *Human gammaherpesvirus 8* (KSHV), 2.0 Å.
   Choices look intentional: one mammalian close homolog + two distant TYMS to span 70-95%.

3. **Alignment quality** — PASS. All four target Cys residues (UniProt 179, 194, 198, 209) align to Cys in all 3 templates. **Caveat**: numbering shifted by one from residue 43 onward.

4. **Heavy-atom clashes (models 1, 3, 10)** — PASS. 0 inter-residue contacts <1.8 Å.

5. **RMSD methodology** — FAIL (statement) / PASS (mechanics). The 0.367 Å is **CA-only**, after outlier rejection (n_atoms_aligned 246-267, not full 287). With 3IHI at 92.7% identity, AutoModel essentially copied 3IHI's backbone with a few side-chain swaps — the RMSD mostly measures "how close 3IHI's backbone is to 1HVY's backbone", not predictive power.

6. **Ramachandran φ/ψ regions** — CONDITIONAL FAIL. Hand-drawn polygons; no separate Gly/Pro maps. Result: 84-85% favoured is *under-estimated* vs PROCHECK.

7. **Per-residue DOPE profile** — PASS. Correct API; profile normalised and smoothed (window 15).

8. **Best-model selection (model 3 vs model 10)** — KEY FINDING. CA-only RMSD m3↔m10 = 0.37 Å, but all-heavy-atom RMSD = **1.01 Å**. Top-divergent CA positions: residues 93-101 and 119-121 — **loop regions where templates were uninformative**. Models differ in how this loop was rebuilt — exactly where one would expect real homology-modelling uncertainty.

## Required corrections
1. **`step1_clean_pdb.py`**: investigate Cys43 missing from FASTA. Restore from SEQRES OR document in FASTA header.
2. **`step5_pymol_compare.py`**: replace `cmd.align(... and name CA)` with both `cmd.super` and full `cmd.align` so all-atom RMSD also reported. Annotate CSV columns as "Cα-only, PyMOL-align outlier-rejected".
3. **`step6_validate.py`**: add separate Gly and Pro Ramachandran polygons (Lovell et al. 2003).
4. **`README_PHASE6.md`**: explicit paragraph that with 3IHI at 92.7% identity the modelling reduces to side-chain optimisation + minor loop refinement.
5. **`summary.json`**: rename `best_by_rmsd_to_crystal` → `best_by_ca_rmsd_to_crystal_pymol_align`.

## Lower-priority improvements
- Compute template-vs-crystal Cα RMSD (3IHI vs 1HVY directly) as the "free lunch" baseline.
- Run MolProbity instead of in-house Ramachandran approximation.
- Add `n_atoms_aligned / n_total_ca` ratio to RMSD CSV.
- Include a more distant template (~50% identity, e.g. *E. coli* TYMS).
