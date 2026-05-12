# Structural Bioinformatics Methods Audit

Agent role: Professional structural bioinformatician (deep technical methods review).
Date: 2026-05-12.

## Methods rigor: **FAIL**

The pipeline runs end-to-end and produces numbers, but the underlying inputs are corrupted in two load-bearing ways: (1) **the MSA is built from non-orthologous sequences** (≥5 of 9 input proteins are NOT thymidylate synthase), so JSD "conservation" is meaningless and the active-site selection inherits that error; and (2) **chain B was discarded** even though the documented dUMP-binding Arg175/Arg176 come from the dimer partner. On top of these, **one mutant (G217W) has heavy-atom overlap of 0.98 Å** (chemically impossible), and the structure cleanup silently dropped a covalently-modified residue (CME43), introducing a 42→44 backbone gap. Docking parameters and Vina invocation are mostly fine in isolation, but garbage-in dominates.

## Detailed findings

### MSA & conservation
- **Wrong orthologs.** Only `P04818` (human) and `P07607` (mouse) are TYMS. `P04394` is *E. coli* ATP synthase, `P04996` is *L. casei* PSII D1, `P0CG53` is yeast polyubiquitin, `O44019` is a *C. elegans* protein, `P11849` is T4 phage PSII-like, `P21520` is a 799-aa Plasmodium protein NOT the bifunctional DHFR-TS (`P13922`). Verifiable from per-row gap-fraction (66–74% for "orthologs", ~17% for the Plasmodium entry — opposite of expected).
- **Catalytic-Cys window is NOT aligned across orthologs.** At alignment column 369 (human C195), only mouse aligns the canonical `LPLMALPPCHALCQFYVV` block. Every other sequence shows unrelated residues or all-gap. JSD scores are dominated by 2-of-9 identity.
- **Plasmodium "DHFR-TS" trimming.** `input.fa` is 799 aa for `P21520`, identical to UniProt; no trimming was performed. Comments claim trimming would happen "after alignment" but no trimming is implemented.
- **Gap-column policy is downweighting, not exclusion.** 575/970 columns (59.3%) are >50% gap — these are not excluded from percentile ranking, inflating apparent significance of low-gap columns.

### Structure prep
- **Chain B removed; dimer partner contacts lost.** Direct check on `1hvy.pdb` shows chain-B Arg175 and Arg176 sit within 5 Å of the chain-A dUMP phosphate. Active-site CSV explicitly annotates R50 and C195 as "in other chain". Mutating residues 175/176 in chain A (as in `R175E_R176E`) does NOT mutate the dUMP-clamping arginines.
- **CME43 silently dropped → backbone gap.** `accept_residue: r.id[0]==' '` rejects HETATMs. CME43 is a covalently-modified Cys; after stripping, residue numbering jumps 42→44.
- **No metals to handle** in 1HVY (verified). BME and CME are S-modifications, not cofactors.
- **dUMP atom inventory is correct.** PDB residue name "UMP" with no O2' atom — matches deoxy-UMP geometry.
- **His tautomer not assigned.** `obabel -h` invoked without `-p`; `protein_h.pdb` has only generic "H" labels (no HD1/HE2 distinction). His196 tautomer is a default guess, not pH-aware.

### Docking parameters (Vina)
- **Box 22³ Å (10 648 Å³).** dUMP itself spans ~12 × 9 × 5 Å — box ~3× larger than ligand. Tolerable for a defined site, but loose; off-site minima are scoreable (modes 11–20 cluster 7–8 Å away from mode 1).
- **Exhaustiveness 16, num_modes 20.** Adequate. Convergence is poor: top mode at −7.73 isolated, modes 2–10 RMSD lower-bound 2.4–4.2 Å vs mode 1, indicating multiple competing minima rather than a clean funnel. Recommend exhaustiveness ≥32 and tighter 18³ box.
- **Seed 42 confirmed.** Both `stage5_6_dock_wt.py:188` and `06_docking_wt/vina_wt.log` ("random seed: 42").
- **Gasteiger charges on the phosphate.** Gasteiger systematically under-polarizes phosphate oxygens; AM1-BCC (Meeko) or RESP would be more appropriate.
- **Receptor side-chain placement** is uncontrolled — no Reduce/PROPKA, no relaxation.

### Mutagenesis quality
Heavy-atom clash check (script in `/tmp/clash_check.py`):
- **Y258A**: 0 clashes (<1.8 Å). PASS.
- **D218K**: 0 clashes. PASS.
- **G217W**: **9 clashes**, including TRP CD1 ↔ VAL223 CG2 at **0.98 Å**, CE2 ↔ VAL223 CG2 at 1.33 Å, CE2 ↔ CB at 1.69 Å. PyMOL Mutagenesis Wizard placed the only acceptable rotamer with no pocket relaxation, giving a sterically impossible structure. Any docking against `G217W_h.pdb` is invalid.

### RMSD vs native dUMP
- Native reference is `03_structure/ligand.pdb` (UMP from chain A, 20 heavy atoms) — correct source.
- `rmsd_to_native()` matches by atom name first, falls back to greedy nearest-pair. Top WT pose RMSD = 1.08 Å, plausible.
- **No Kabsch superposition.** RMSD computed in absolute crystal frame — appropriate for "docked pose vs crystal in same coordinate system" but undocumented.

## Specific corrections needed

1. **Replace `01_msa/input.fa`.** Real TYMS orthologs only: human P04818, mouse P07607, rat P45352, *E. coli* P0A884, *L. casei* P00469, *S. cerevisiae* P07807 (CDC21), *C. elegans*, Drosophila Q9V3K2 (verify), *P. falciparum* P13922 trimmed to TS domain (~residues 280–608). Drop the polyubiquitin/photosystem entries. Re-run Stage 1.
2. **Trim PfDHFR-TS to TS domain BEFORE alignment.**
3. **Exclude columns with >50% gap from percentile ranking.**
4. **Keep both chains** in `scripts/stage3_structure.py` (`accept_chain: chain.id in ('A','B')`). Re-prepare receptor for both WT and all 175/176 mutants.
5. **Preserve CME43** by re-mutating it back to CYS programmatically (rename, strip the 2-hydroxyethyl atoms) instead of dropping it.
6. **Drop or rebuild G217W** with side-chain relaxation (PyRosetta `pack_rotamers` or FoldX `RepairPDB`) before re-docking.
7. **Switch ligand charges to AM1-BCC** via Meeko `mk_prepare_ligand` as the primary path.
8. **Add His tautomer / PROPKA assignment** before adding hydrogens (e.g. `pdb2pqr30 --ff=AMBER --with-ph=7.4`).

## Lower-priority improvements

- Tighten Vina box to 18³ Å and raise exhaustiveness to 32; report cluster populations not just affinities.
- Add Kabsch-aligned RMSD as a second column for comparison.
- Validate atom-name match success rate in `rmsd_to_native()`.
- For double mutants, dock with explicit `flex` side-chains at the mutated residues to absorb local rearrangement.
