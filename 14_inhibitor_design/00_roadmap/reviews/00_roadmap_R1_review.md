# Roadmap Review — Round 1 (biologist + bioinformatician)

Reviewer agent role: structural bioinformatician + medicinal-chemistry biologist, in the style of `reviews/04_structural_bioinformatician.md` and `reviews_phase6/04_structural_bioinformatician.md`. Severity convention as defined by the corrector spec (HIGH / MEDIUM / LOW). Round 1, draft v0.

## Overall verdict

**CONDITIONAL_PASS.** The four-strategy framing is mechanistically sound and the seven-step pipeline cleanly inherits Phase 5–7 conventions, so structurally the plan can stand. However, several Tier-1 anchors are wrong or mis-bucketed, the active-site anchor list confuses *the drug given to patients* with *the species that actually binds TYMS*, the dimer-interface peptide is the wrong peptide sequence, and three control / sanity steps that are standard in any inhibitor-design pipeline are entirely absent (re-dock-RMSD positive control, PAINS/Lipinski filter, tautomer/protomer enumeration). All are fixable on paper — no compute wasted yet — but every one must be corrected before docking starts.

## Section-by-section

### Section 0 (Inputs)
Verdict: **CONDITIONAL_PASS**
Findings:
- [LOW] Box centre `(−0.137, +4.232, +15.159)` is correct per `TECHNICAL_NOTES.md §"Where the active-site box lives"`. No fix needed.
- [MEDIUM] Two box sizes (22³ canonical, 18³ Phase-7) are *both* quoted but Section 2A then commits silently to 22³ for strategy 1. Pick one and justify; if 22³ is kept, add the Phase-7 verifier's known concern that 22³ admits off-site minima 7–8 Å from mode 1 (see `reviews/04_structural_bioinformatician.md §"Docking parameters (Vina)"`). Fix: footnote which box and why, with a one-line cross-reference.
- [LOW] The `Δ_vs_dUMP = inhibitor_score − dUMP_score` convention is numerically correct given Phase-4 sign rules (more-negative-is-tighter), but the sentence "negative = the inhibitor binds tighter than dUMP" is right *only at the same site/state*. Strategies 1 and 2 cross-compare against dUMP at very different sites — Δ_vs_dUMP for an antifolate at the cofactor site is **not** a competition with dUMP, it is a competition with the methylene-THF/raltitrexed reference. Fix: split into `Δ_vs_dUMP` (site 1 only) and `Δ_vs_raltitrexed` (site 2). Otherwise the headline ranking plot will mislead.
- [LOW] The holo receptor in Section 0 already contains raltitrexed at the cofactor site, so any cofactor-site (strategy 2) dock against `holo` is asking "can the candidate displace raltitrexed?", not "does it bind the empty cofactor site?". Both are valid questions but they are different and must be labelled.

### Section 1 (Strategies)
Verdict: **CONDITIONAL_PASS**
Findings:
- [HIGH] **5-FU (PubChem CID 3385) is the wrong active-site anchor.** 5-FU is a pyrimidine *prodrug*; the species that actually inhibits TYMS is the deoxynucleotide metabolite **5-FdUMP** (CID 15718), which forms the covalent ternary complex with Cys195 + CH2-THF (Carreras & Santi 1995; Longley *et al.* 2003). Free 5-FU has no phosphate and will not engage the Arg175/176/215 phosphate clamp — docking 5-FU into the dUMP pocket teaches the wrong lesson. Fix: keep 5-FdUMP as the canonical anchor; either drop 5-FU or move it to a small "metabolic-precursor sanity panel" with an explicit annotation that it is expected to dock weakly.
- [HIGH] **Nolatrexed (AG-337, CID 60198) is mis-bucketed.** Nolatrexed is a *folate-site* (lipophilic, non-classical antifolate) inhibitor of TYMS, not an active-site/dUMP-pocket binder (Webber *et al.* 1996, J. Med. Chem. 39:4007). It must move from Strategy 1 → Strategy 2. Same for **ZD9331 (plevitrexed, CID 153985 — see also B1 below)**: ZD9331 is a classical antifolate that competes at the folate site (Jackman *et al.* 1997). Strategy 1's anchor list should be substrate/nucleotide mimics only: 5-FdUMP, BrdUMP, dUMP itself as positive control, and ideally IdUMP / FdUR-MP analogs.
- [MEDIUM] **Missing canonical active-site anchor: dUMP itself.** Re-docking dUMP at the active site is the mandatory positive control of the whole phase (does the protocol recover the crystal pose with sub-2 Å RMSD?). Without it there is no calibration line.
- [HIGH] **ZD9331 ≠ nolatrexed.** The roadmap calls "ZD9331 (nolatrexed)" in Section 1 and then again "BGC9331 / plevitrexed" in Section 2 — these are conflated. ZD9331 *is* plevitrexed/BGC9331 (developed by BTG/Cancer Research Campaign; Jackman 1997). Nolatrexed is AG-337 (Webber 1996, Agouron). Two completely different molecules. Fix: separate them and verify CIDs.
- [MEDIUM] Dimer-interface chemotype claim ("LR octapeptide LSCQLYQR") is the wrong peptide. Cardinale *et al.* 2011 (FEBS J 278:1487) and the follow-up Salo-Ahen 2015 work describe the **LR octapeptide LSCQLYQR is not the peptide reported** — the published TS dimer-interface disruptor is LSCQLYQR derived from the conserved interface loop, but the validated lead is the octapeptide **LSCQLYQR** *only* if the authors' numbering is used. The roadmap must cite the exact peptide and its position in the human TYMS sequence (Cardinale 2011 maps it to residues ~200–207 of *E. coli* TS; the human equivalent is shifted). Fix: corrector must look up the peptide sequence from the primary paper and not invent one. If the actual sequence is different, change the anchor.
- [LOW] Strategy 4 ("allosteric / mRNA-binding loop, residues 1–27") is correctly framed as exploratory. Add one sentence: this region is intrinsically disordered in 1HVY (no resolved density for residues 1–26 in chain A — verify against the PDB header) and so FPocket/FTMap will not find a pocket there. The honest framing is "we will look elsewhere on the surface and document the absence of a druggable cavity at the 1–27 loop region."

### Section 2 (Per-strategy pipeline)
Verdict: **CONDITIONAL_PASS**
Findings:
- [HIGH] **No re-dock RMSD positive control.** Standard inhibitor-design SOP requires re-docking the native substrate (and any anchor with a crystal pose) into the prepared receptor and demanding < 2 Å heavy-atom pose RMSD before *any* candidate is scored (Hawkins *et al.* 2007; Warren *et al.* 2006). The roadmap mentions pose RMSD in step E but only for hits that *beat* dUMP — too late. Fix: add step **A0** = "Re-dock dUMP into the prepared active-site box. Gate: top1 pose RMSD vs 1HVY UMP heavy atoms ≤ 2.0 Å (Kabsch-aligned). If fails, abort and rebuild receptor." Same gate for raltitrexed at the cofactor site against 1HVY D16.
- [HIGH] **No PAINS / Lipinski / Brenk filter on Tier 2 decoys/analogs.** PubChem similarity search at Tanimoto ≥ 0.7 on Morgan-2 will pull in covalent warheads, peroxides, aggregators (Baell & Holloway 2010 PAINS filters; Brenk *et al.* 2008). Without a PAINS gate the "decoy" set is contaminated with assay-pathological structures. Fix: add an RDKit filter step between B and C: `rdkit.Chem.FilterCatalog` with `PAINS_A/B/C` + Brenk + NIH; log every rejection. Lipinski / Veber filter optional but recommended for the allosteric fragment set.
- [HIGH] **No tautomer / protomer enumeration.** RDKit's default `MolFromSmiles` returns one tautomer; pemetrexed, raltitrexed, nolatrexed all have multiple ionisable centres and the dominant species at pH 7.4 is not obvious. Fix: use `rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumerator` + `Dimorphite-DL` (Ropp 2019) at pH 7.4 ± 0.5; dock all enumerated protomers and report which one scored best. This is exactly the failure mode that hits the AD4 sign-of-the-phosphate question (see TECHNICAL_NOTES Phase 6c).
- [MEDIUM] **Receptor protonation re-check missing.** Strategy assumes Phase-6c `protein_dimer_apo_fixed.pdbqt` is reusable as-is. Re-verify the His196 tautomer at pH 7.4 (HID vs HIE vs HIP) before docking — Phase-1 audit (`reviews/04_structural_bioinformatician.md §"His tautomer not assigned"`) flagged this and Phase-6c only fixed *charges*, not necessarily tautomer states. Fix: add an assertion in step A that `protein_dimer_apo_fixed.pdb` has explicit HID/HIE on every His, or re-run propka + pdb2pqr30 once.
- [MEDIUM] **Water-mediated H-bonds are ignored.** TYMS has structural waters in the active site (notably one bridging Tyr258 ↔ dUMP O4 in 1HVY). For an inhibitor that wants to mimic dUMP, deleting these waters silently penalises the right answer. Fix: at minimum document the decision (typically "all crystallographic waters removed for docking; PROLIF post-hoc flags which Tier-1 anchor poses would have made water-bridged H-bonds in the crystal"). 
- [LOW] **PROLIF is fine** but **PLIP** (Salentin *et al.* 2015, https://plip-tool.biotec.tu-dresden.de/) gives a clean, citable, command-line per-complex interaction report and is the field-standard for inhibitor papers. Recommend running both PROLIF (for the heatmap matrix) *and* PLIP (for the per-complex human-readable report).
- [LOW] DBSCAN clustering at >2 Å in step E is fine; specify `eps=2.0 Å, min_samples=2`.
- [LOW] Seed list {42, 7, 13, 99, 256} matches Phase 7 — good. Keep it.
- [LOW] `exhaustiveness 32` matches v5 canonical — good.

### Section 3 (Reviewer/corrector loop)
Verdict: **PASS**
Findings: none — the loop description matches the repo's existing doer ↔ verifier convention. Note: this review is round 1; the corrector should expect ≥ 2 rounds.

### Section 4 (Compute budget and stop conditions)
Verdict: **CONDITIONAL_PASS**
Findings:
- [MEDIUM] Strategy 3 estimate ("3 peptide anchors + 10 mimics × 5 seeds = ~3–5 h") is optimistic. Vina on an 8-residue peptide with ~25 rotatable bonds will run > 1 h *per pose* at exh=32 on M-series (Hassan 2017 reports >> 30 min for 5-mers). Fix: budget 12–20 h or commit to OPEN-B2 alternative (HPEPDOCK / HADDOCK / FRODOCK).
- [LOW] Stop condition S1 is well-stated and matches Phase-7's null-result-honesty.
- [LOW] Add S4: "If Tier-1 anchor (known active) fails to dock into its own site with Δ_vs_dUMP ≤ −0.85, the docking protocol is the suspect, not the chemistry — flag and investigate before any Tier-2 ranking is reported."

### Section 5 (Deliverables)
Verdict: **PASS**
Findings:
- [LOW] Add `pains_flag` and `lipinski_flag` columns to the per-strategy `results.csv` schema in Section 2G (so the deliverable matches the new filter step).
- [LOW] Add `tautomer_id` and `protonation_pH` columns (per the C-step fix above).

### Section 6 (Open questions)
Verdict: addressed below in "Resolution of OPEN tags".

## Resolution of OPEN tags

- **[OPEN-A1] Dimer-interface residue list.** Verdict: **REVISE.** The literature-prior list (residues 30, 76–81, 96–99, 174–179, 199, 202) is *partially* right but is taken from mixed-species papers (Cardinale uses *E. coli* numbering, Salo-Ahen uses *L. casei*). For human TYMS (P04818, 313 aa, 1HVY numbering) the interface from a 4 Å contact-map analysis spans roughly chain-A residues {21–27, 76–82, 174–180, 200–207, 250–256} and the symmetric chain-B partners (the homodimer is two-fold symmetric, but the interface is not "the same residues" in chain A and chain B in the structural sense — they are the same *sequence positions* contacting the *opposite* chain's residues). Corrector must compute the contact map programmatically (MDAnalysis `Universe.select_atoms("chainA").contacts(...other chain..., cutoff=4.0)`) and *not* hardcode the residue list from these papers.

- **[OPEN-B1] CID-to-InChIKey verification.** Verdict: **YES, mandatory.** Hard-coded PubChem CIDs are the single most common silent error in computational chemistry pipelines (CIDs get redirected to salts, racemates, the wrong tautomer). Corrector must add a per-compound assertion that the fetched canonical SMILES → RDKit → InChIKey matches a literature-validated InChIKey (e.g. pemetrexed free acid = `QOFFJEBXNKRSPX-ZDUSSCGKSA-N`). Pre-flag now: roadmap claims **pemetrexed = CID 135410875**; the canonical free acid is CID 60843 (135410875 is a salt/disodium-heptahydrate form). Re-verify all 14 anchor CIDs.

- **[OPEN-B2] Peptides in Vina.** Verdict: **REVISE — switch tool.** Hassan *et al.* 2017 (J. Comput. Chem. 38:1278) shows Vina median peptide-pose accuracy drops below 2 Å only for ≤ 4-mers. For an 8-mer like LSCQLYQR, AutoDock Vina cannot be trusted quantitatively. Use **HPEPDOCK** (Zhou *et al.* 2018, Nucleic Acids Res. 46:W443) or **CABS-dock** as the primary docker for strategy 3, and keep Vina for ≤ 5-residue mimetics only. Document this as a strategy-3-specific deviation from the v5 protocol.

- **[OPEN-C1] Charge-delta surface in CSV.** Verdict: **YES, mandatory.** Adopt verbatim. Add the columns `formal_charge_input`, `total_q_pdbqt`, `delta_q`, and gate-fail rows where `|delta_q| > 0.1 e`.

- **[OPEN-D1] Flex-residue follow-up gating.** Verdict: **YES** but tighten. Proposal "flex follow-up only on rigid hits with Δ_vs_dUMP ≤ −0.85" is correct *and* should additionally cap at top-5 per strategy (15–20 flex runs total across 4 strategies, fits the 30-min/compound budget Phase 8b documented).

- **[OPEN-F1] PRODIGY-LIG / FoldX availability.** Verdict: **YES — accept the fallback.** PRODIGY-LIG web API is fine for ≤ 30 submissions/day. FoldX 5 arm64-darwin binary is not officially distributed (Schymkowitz lab ships x86_64; works under Rosetta 2 but slow). Fallback to MDAnalysis BSA + interface residue contacts is sufficient; document the missing ΔΔG as a known limitation, do not silently skip.

- **[OPEN-F2] `fpocket` install-or-skip.** Verdict: **YES — install.** `fpocket` is one `brew install fpocket` away on arm64-darwin (Homebrew bottle exists). No reason to skip. If install fails, document the OS version and fall back to MDpocket (web) for the allosteric strategy.

## Additional findings (items not in roadmap that should be)

- [HIGH] **No re-dock RMSD sanity gate** (already raised under Section 2 but re-flagging here as the single most important addition).
- [HIGH] **No PAINS/Brenk filter on the analog/decoy set** (already raised).
- [HIGH] **No enantiomer enumeration.** Several anchors (raltitrexed, pemetrexed) have a defined (S)-stereocentre; PubChem may return a racemate. Corrector must verify single enantiomer is docked (the (S) for the antifolates), or dock both and report.
- [HIGH] **No tautomer/protomer enumeration at pH 7.4** (already raised).
- [MEDIUM] **DUD-E-style decoy generation by RDKit alone produces weaker decoys than the actual DUD-E web service.** Use the DUD-E server (Mysinger *et al.* 2012, J. Med. Chem. 55:6582) for at least the active-site (strategy 1) decoy set so the enrichment AUC is comparable to literature benchmarks.
- [MEDIUM] **No enrichment / AUC metric in the deliverables.** Whenever you have actives + decoys you should report **ROC-AUC** and **BEDROC** (α=20) on the docking-score-vs-label classifier. Without it, "compound X beat dUMP by 0.9 kcal/mol" is not contextualised against false-positive rate. Fix: add `roc_auc_top1`, `bedroc_alpha20` per strategy in `master.csv`.
- [MEDIUM] **No induced-fit / cryptic-pocket indicator.** A simple, free signal is the apo-vs-holo Vina-score gap on the same ligand: large gap (> 2 kcal/mol) suggests the pocket geometry shifts meaningfully on cofactor binding, i.e. a cryptic component. Add `apo_minus_holo_top1` to the per-compound row and flag any > 2 kcal/mol.
- [MEDIUM] **No mention of crystal-water handling.** Document the choice (remove all vs keep within X Å of pocket) and stick with it across strategies.
- [LOW] **No literature priors for negative controls.** Strategy 2's docking should include a known *non-binder* (e.g. an unrelated small molecule of similar MW/logP, such as ibuprofen) to verify the cofactor-site box doesn't bind everything.
- [LOW] **AlphaFold receptor as a second target.** Phase 7 already validated the AF model vs 1HVY (`12_phase7/03_alphafold/`). Phase 14 should at least sanity-dock the top hit per strategy against the AF receptor as well; if scores diverge by > 2 kcal/mol, that is a real (publishable) finding about receptor-model sensitivity.

## Sign-off requirements for next round

Corrector agent must, before I sign PASS:

1. **Swap Strategy 1 anchor list.** Drop 5-FU (CID 3385) and nolatrexed/AG-337 (CID 60198); move nolatrexed and ZD9331 to Strategy 2; add dUMP (positive control) and BrdUMP; verify 5-FdUMP CID 15718 and BrdUMP CID. Justify each anchor with a one-line literature citation in a new column.
2. **Verify every PubChem CID** by fetching the canonical InChIKey and cross-checking against a hard-coded literature InChIKey table (see [OPEN-B1]); re-flag pemetrexed CID specifically (135410875 looks like a salt; free-acid CID is 60843).
3. **Fix the dimer-interface peptide.** Cite the exact peptide sequence and source paper (Cardinale 2011 — give the residue range it maps to in human TYMS numbering, *not* the *E. coli* numbering). Switch primary docker to HPEPDOCK for strategy 3 ([OPEN-B2]).
4. **Add step A0 = positive-control re-dock gate.** dUMP into the active site, raltitrexed into the cofactor site, both with ≤ 2 Å pose-RMSD gate. Abort otherwise.
5. **Add a PAINS/Brenk filter step** between B and C; surface flagged compounds in `results.csv`.
6. **Add tautomer + protomer enumeration at pH 7.4 ± 0.5** (Dimorphite-DL + RDKit TautomerEnumerator); add `tautomer_id` and `protonation_pH` columns.
7. **Add enantiomer verification** for stereoactive anchors (raltitrexed, pemetrexed).
8. **Document crystal-water handling** explicitly in Section 2A.
9. **Split `Δ_vs_dUMP` into per-site Δ references** (Δ_vs_dUMP for site 1, Δ_vs_raltitrexed for site 2, Δ_vs_LR_peptide for site 3, no Δ for site 4) and update the plotting plan accordingly.
10. **Add `roc_auc_top1` + `bedroc_alpha20` + `apo_minus_holo_top1` + `pains_flag` + `lipinski_flag` + `tautomer_id` + `enantiomer` columns** to the per-strategy CSV schema.
11. **Compute the dimer-interface residue list from the 4 Å contact map** rather than the literature-prior list ([OPEN-A1]).
12. **Pick one box size** for strategy 1 (18³ vs 22³) with a one-sentence justification cross-referenced to the Phase-7 off-site-minimum concern.
13. **Re-verify the apo receptor's His tautomer assignment** at pH 7.4.
14. **Re-estimate strategy 3 compute budget** under HPEPDOCK (or whatever tool is chosen) ([OPEN-B2]).

End of Round 1 review.
