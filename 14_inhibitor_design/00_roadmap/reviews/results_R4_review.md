# Phase 14 Results Review — Round 4 (post-execution)

## Overall verdict
**CONDITIONAL_PASS** — the pipeline executed end-to-end, the four CSVs and master union are consistent, the headline numbers are chemically defensible, and the README's tone is teaching-oriented. Three concrete issues block an unconditional PASS: (i) the `sasa_buried_pct` column in `01_active_site/results_analysed.csv` returns values 1.05–1.46 (i.e., > 1, which the README itself flags as a bug indicator) — the metric is broken and should either be fixed or removed from the analysed CSV before the figures are claimed to depend on it; (ii) the Strategy 1 A0 re-dock gate failed numerically (RMSD 5.83 Å vs 2.0 Å) and the "frame-mismatch" defence, while plausible, is not actually demonstrated anywhere — a one-line `analysis_post.py` extension extracting dUMP coordinates from the Phase-6c-frame receptor and re-running the RMSD would prove the point; (iii) `results_analysed.csv` in S1 carries duplicate rows for three failed decoys (`decoy_CID49846` appears twice with `ok=False`) — harmless but noisy, and indicates a non-deduplicated build step.

## Per-strategy verdict + biology defence

### Strategy 1 (active-site)
**verdict:** CONDITIONAL_PASS.

**biology:** Numbers are chemically sensible. The Tier-1 nucleotide cluster (dUMP −8.78, 5-FdUMP −9.04, BrdUMP −8.88) sits within 0.27 kcal/mol of each other — entirely consistent with the published TYMS literature: the 5-halogen at C5 of the uracil ring doesn't change the binding pose, it changes the *covalent* step (Cys195 Michael addition) that Vina does not model. The ~1.5 kcal/mol drop from 5-FdUMP (−9.04) to floxuridine (−7.48) cleanly attributes phosphate-clamp engagement (Arg23, Arg175', Arg176', Arg215, Ser216) at ~1.5 kcal/mol of binding energy — defensible against the Phan/Stout 2001 1HVY arginine-clamp characterization. The 5-FU prodrug (−4.95) at the floor is correct: 5-FU has no sugar, no phosphate, nothing for the clamp to grab — it's a *prodrug* whose active species is 5-FdUMP (only formed after thymidine kinase + thymidylate phosphorylase). The README correctly labels it as "precursor sanity panel" — exactly the right framing. The ~3 kcal/mol Tier-1-vs-prodrug gap matches the empirical "Vina can rank active vs decoy but not active vs active" pattern Phase 7 already documented.

**methods:** The A0 gate is the weakest part. RMSD 5.83 Å vs a 2.0 Å gate is a *failed* gate per the roadmap's own pass criterion, and the defence ("frame mismatch between pre-prep crystal coords and Phase-6c-frame receptor") is a hypothesis, not a proof. The score-vs-Phase-7 match (−8.78 vs −8.785, Δ 0.005 kcal/mol) is *strong* circumstantial evidence the docked pose is the canonical pose, but a 5.8 Å RMSD against a correctly-frame-aligned calibrant would falsify that. **Required fix:** add `A0_redock_gate/A0_frame_check.py` that extracts dUMP from `06f_receptor_fixed/protein_dimer_holo_fixed.pdbqt` (or equivalent Phase-6c-frame source) and re-computes RMSD; if it then passes <2 Å, the README can keep its current narrative; if it still fails, the narrative needs to soften to "score-equivalent, pose-uncertified."

**caveats correctly disclosed?:** Mostly. The README's Limitation #1 (Vina ±0.85 noise floor) and #5 (waters removed → E1b annotation) are appropriately surfaced. The SASA-buried-% bug is not mentioned in the README at all — values of 1.0–1.5 are mathematically impossible for a fraction, and the README text claims to use this column. Either the column is a per-pose ΔSASA in Å² mis-scaled by area normalization, or it's a freesasa output mis-divided. **Required fix:** either correct the SASA calc or drop the column from the analysed CSV and the README's claim.

### Strategy 2 (cofactor-site)
**verdict:** PASS.

**biology:** Plevitrexed (ZD9331) at −10.01 *is* defensible. ZD9331 was specifically designed as a non-polyglutamatable TYMS inhibitor (replaces the glutamate γ-carboxylate of raltitrexed with a quinazoline-tethered propargyl chain) and its published Ki against human TYMS is sub-nanomolar (Jackman et al. 1997 — actually tighter than raltitrexed). A 0.88 kcal/mol Vina advantage over raltitrexed is therefore directionally right, even if "just above noise" at the rigid-Vina scale. The pemetrexed S vs R observation (−9.72 vs −9.63) is the *correct* null finding for a rigid-receptor docker without a chiral scoring term — Vina's force-field-based scoring is achiral by construction; the README's framing ("Vina rigid-receptor cannot distinguish enantiomers here, which is a known limitation") is technically accurate. Methotrexate at −9.59 is reasonable: MTX *is* a weak TYMS binder (its primary target is DHFR) but it's structurally a 2,4-diaminopteridine antifolate and engages the same folate-binding residues; the Vina score reflects shape complementarity, not selectivity.

**methods:** Box centre (0.401, 12.392, 17.766) is correct — the README documents it as the D16 heavy-atom centroid from the holo cofactor PDBQT, which is exactly where raltitrexed sits in 1HVY chain A (Phan 2001 reports raltitrexed contacting Asp218, Phe225, Leu221, Met311, Tyr135 — a pocket centred ~12 Å from active-site Cys195, consistent with the box centre offset from the S1 box centre of about 10 Å). Exhaustiveness 16 (vs canonical 32) is documented in TECHNICAL_NOTES with user authorisation and the raltitrexed self-control still passes (−9.13 apo, score reproducible across seeds within 0.22 kcal/mol). This is sound.

**caveats correctly disclosed?:** Yes — the "holo state penalty is brutal" teaching point (drop of ~3 kcal/mol for every antifolate) is the correct interpretation of displacement-vs-empty-pocket dynamics and the README states this plainly.

### Strategy 3 (dimer-interface)
**verdict:** PASS as a documented null result.

**biology:** The technical interpretation of the +86 / +85 positive scores is correct: positive Vina scores mean the search algorithm could not find any conformation with favourable van-der-Waals + H-bond + electrostatic terms inside the search box, which for a 938 Da, 70-heavy-atom 8-mer in a 26×22×22 Å box is *exactly* the expected failure mode for rigid-receptor small-molecule Vina (Hassan 2017 documents this for any peptide ≥ 6 residues). The "specificity vs scrambled = +1.48 kcal/mol" framing in the README is honest: the canonical peptide actually scored *slightly worse* than the scrambled control, which by the strategy's own pre-registered logic (canonical should beat scrambled if there's a real interface preference) is a clean null. **The README is correct to call this a "documented null result, not a finding."**

The contact-map result (46 chain-A interface residues, 42 chain-B) is in the right order of magnitude for the TYMS dimer interface — Cardinale et al. 2011 (PMID 21282636) and Salo-Ahen et al. 2015 (PMID 26390032) work in E. coli numbering on a smaller homodimer, but the cross-monomer arginine residues (Arg175', Arg176' in human, R166'/R167' equivalents in E. coli) reaching into the partner active site is preserved chemistry — the ~40-residue contact zone on each monomer is consistent with the published interface footprint.

**methods:** Box centre and size are reasonable. The 4-mer fragment scores (−4.1 to −4.7 kcal/mol, near-degenerate across positions) consistent with shallow PPI binding — exactly what one would predict for a flat interface with no pre-formed druggable hot spot.

**caveats correctly disclosed?:** Yes — both the README and TECHNICAL_NOTES surface HPEPDOCK/CABS-dock outage and the "right tool for the question" disclaimer. `unreliable_flag = True` is set in the per-row CSV for the 8-mer rows. This is appropriate disclosure.

### Strategy 4 (allosteric)
**verdict:** CONDITIONAL_PASS — the framing is honest but understated.

**biology:** The −4 to −5.5 kcal/mol band for all 60 fragment runs is the correct interpretation: this is the absolute scale of "drug-like fragment bouncing off a non-pocketed surface patch" — *not* the −8 to −10 kcal/mol scale of true pocket binding. The README's conclusion that "TYMS does not present an obvious druggable allosteric pocket on its resolved chain-A surface outside the substrate and cofactor sites" is supported by the data.

However, the freesasa-fallback site-picking method is **structurally biased** in two ways the README does not flag: (i) picking the *highest-SASA Cα* selects for solvent-exposed loop residues, which are by construction the *worst* candidates for druggability (druggable pockets are concave, partially buried, with mixed hydrophobic/H-bond character — not high-SASA bumps); (ii) the 15 Å exclusion zone around the active and cofactor sites guarantees the picks miss any *cryptic* pocket adjacent to but not overlapping the known sites, which is precisely where the literature (e.g., the Phan 2001 "open-loop" conformation, Anderson 2009 conformational ensembles) suggests TYMS allosteric modulation may live. The negative result is therefore weaker evidence of "no druggable allosteric site" than the README implies — it is evidence of "no druggable allosteric site *at three high-SASA surface centroids ≥ 15 Å from substrate/cofactor*."

**methods:** The cavity centres ((−13, +21, −8.5), (−21, +5, +9), (+18, +27, −3.5)) span the chain-A surface reasonably but, by selecting `res26`, `res42`, `res284`, hit the N-terminal arm, a peripheral loop, and a C-terminal region — none of which are mechanistically privileged for TYMS regulation. The negative result is real for *these three patches*, not for the protein at large.

**caveats correctly disclosed?:** Partially. The README says "freesasa-fallback method tells us *where the spatial candidates are* but not *whether they're druggable*" — good, but it should also explicitly add: *"Strategy 4 should be re-run on x86-64 with working FPocket (or with PrankWeb / DoGSiteScorer) before any 'no druggable allosteric pocket' claim can be made strongly; the present result speaks only to the three surface centroids tested."*

## Cross-strategy consistency

- **Noise floor consistency:** The ±0.85 kcal/mol Vina noise floor is invoked consistently in S1 (dUMP / 5-FdUMP separation declared sub-noise), S2 (plevitrexed declared "just above noise"), and S4 (entire band declared "well below the Tier-1 threshold"). S3 correctly does not apply it (peptide engine scale is non-comparable). Consistent.
- **Δ-reference choices:** S1 → dUMP (correct: substrate is the natural reference at the substrate pocket); S2 → raltitrexed (correct: the only co-crystallised cofactor-site antifolate in 1HVY); S3 → scrambled-sequence permutation control (correct for PPI null hypothesis); S4 → absolute Vina score with `delta_vs_reference` left null (correct for an exploratory screen with no Tier-1). All four references are the *right* ones for their respective questions.
- **Master CSV union:** Row count is 86 data rows = S1 (16, after dropping 6 failed decoys) + S2 (18) + S3 (7) + S4 (45). Column union is wide (38 columns) and per-strategy columns are appropriately NaN-filled outside their strategy. No rows are silently dropped. Faithful union.
- **One small issue:** S1's `results_analysed.csv` has three failed decoys (`decoy_CID49846`, `decoy_CID6149`, `decoy_CID65349`) each appearing in duplicate or quadruplicate `ok=False` rows. They are correctly excluded from the master (master only carries `ok=True` rows for S1), but the duplication suggests `analysis_post.py` is being re-run without deduplication. Cosmetic.

## README educational tone

The separation between README (teaching) and TECHNICAL_NOTES (audit/build log) is **clean and appropriate**. The README leads each strategy with "What we did → What we found → Teaching points" — a textbook pedagogical structure. The Teaching Points are genuinely instructional ("the canonical 5-fluoro substitution is barely visible at the Vina rigid-receptor scale", "the decoy/weak-binder separation IS clean — that ~3.5 kcal/mol active-vs-prodrug gap tracks the chemistry", "holo state always weaker because the cofactor occupies geometry the substrate would otherwise use") and would be defensible in a medicinal-chemistry course handout.

The Honest Limitations section at the end is the right length, lists the six material caveats (Vina noise floor, DUD-E outage, HPEPDOCK outage, FPocket crash, waters removed, no GNINA), and explicitly defers the "agent-grade caveats" to TECHNICAL_NOTES — exactly the separation the user requested. TECHNICAL_NOTES contains the wrong-CID audit table, the per-strategy execution caveats, and the cross-strategy ranking limitations — which is the right *build-log* content and would clutter the README if mixed in.

The README does *not* read as a build log. It reads as a methods-and-results write-up for a docking screen. Pass on tone.

## Sign-off requirements for R5

1. **Fix or drop `sasa_buried_pct`** in `01_active_site/results_analysed.csv`. Values are 1.05–1.46 which are mathematically impossible as fractions. Either rename to `sasa_buried_A2` (and re-document units), divide by total ligand SASA to get a true fraction, or remove the column. The README's Strategy 1 narrative does not currently quote a SASA-buried number, so dropping the column is the lowest-risk fix.
2. **Demonstrate the A0 frame-mismatch defence.** Add `A0_redock_gate/A0_frame_check.py` that extracts dUMP coordinates from a Phase-6c-frame source (e.g., `06f_receptor_fixed/cofactor_A.pdbqt`-style file containing the bound dUMP if one exists, or `08_*/holo_apo_extract/`) and re-runs RMSD vs the docked pose. If RMSD now passes < 2 Å, log "A0 gate passes after frame alignment." If it still fails, soften the README's "matches the Phase-7 canonical exactly" language to "score-equivalent; pose un-validated due to RMSD calibrant frame mismatch."
3. **Deduplicate `01_active_site/results_analysed.csv`.** Three failed decoy rows appear 2–4 times each. Add a `drop_duplicates(subset=["compound","state","seed"])` at the end of `analysis_post.py`.
4. **Strengthen Strategy 4 caveat in README.** Add one sentence after the Strategy 4 teaching point: *"This negative result speaks only to the three high-SASA surface centroids tested; a definitive allosteric-pocket survey requires FPocket (broken on this host), PrankWeb, or DoGSiteScorer on x86-64 hardware."*
5. **Optional:** consider re-rendering `fig1_distributions.png` and `fig4_tier_separation.png` after the dedupe in (3) so the violins/boxes are not subtly weighted by the duplicated failed-decoy rows. Likely cosmetic but worth checking.

With items 1–4 addressed, this is a publishable result. Item 5 is housekeeping.

Files reviewed:
- `/Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/README.md`
- `/Users/ario/Downloads/aminak-inhibitor/TECHNICAL_NOTES.md` (Phase 14 section, lines 148–199)
- `/Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/01_active_site/results_summary.csv`
- `/Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/01_active_site/results_analysed.csv`
- `/Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/02_cofactor_site/results_summary.csv`
- `/Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/03_dimer_interface/results_summary.csv`
- `/Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/04_allosteric/results_summary.csv`
- `/Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/05_aggregate/master.csv`
