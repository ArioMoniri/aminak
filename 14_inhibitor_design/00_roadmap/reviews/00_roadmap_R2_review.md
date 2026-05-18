# Roadmap Review — Round 2 (biologist + bioinformatician)

Reviewer agent role: structural bioinformatician + medicinal-chemistry biologist, in the style of `reviews/04_structural_bioinformatician.md`. Severity: HIGH / MEDIUM / LOW. Round 2, against `ROADMAP.md` v1 (post-R1 corrector pass).

## Overall verdict

**CONDITIONAL_PASS.** The corrector made an honest, point-by-point pass and 12 of 14 R1 items are genuinely CLOSED. v1 is markedly better than v0 — the strategy framing is now chemically correct, Δ references are per-site, the A0/A1/A2/A3 prefix gates are well-structured, and the budget revision (20–35 h → 50–80 h) is honest. However, two R1 items are only *partially* closed (CID verification is deferred to runtime with placeholder text that includes the guess word "likely", and PROLIF is asked to do something it physically cannot do), and three NEW v1 additions introduce defects that would bite a fresh agent at runtime (HPEPDOCK has no fallback / no per-job timeout, the cofactor-site box centre is ambiguous for the apo case, PLIP install needs an explicit OpenBabel note for arm64-darwin Python 3.14). Fixable on paper. No compute wasted.

## Status of R1 sign-off items

1. **Swap Strategy 1 anchor list** — CLOSED — 5-FU moved to precursor sanity, nolatrexed and ZD9331 relocated to S2, dUMP + BrdUMP added; per-anchor literature citation present in §1.1 column.
2. **Verify every PubChem CID by InChIKey** — OPEN — §A2 is structured correctly but four of nine rows in the verification table contain the literal string `*(to verify A0.2)*` and the BrdUMP row says `*(to verify A0.2; likely 167671 for free acid)*`. "Likely" is not an InChIKey. The whole point of A2 is that the corrector pre-verifies the literature InChIKey table *now*, on paper — runtime can only check whether the CID maps to the literature InChIKey if the literature InChIKey is hard-coded. Deferring to runtime defeats the gate. (See sign-off #1 below.)
3. **Fix dimer-interface peptide + switch to HPEPDOCK** — CLOSED with caveat — §1.1 row 11 + §A3 + §D Strategy-3 deviation are well-specified. Caveat: the actual peptide sequence is still not written down anywhere; A3 promises to compute it from a Biopython pairwise alignment to Cardinale 2011. That is acceptable *only if* the Cardinale 2011 *E. coli* peptide is itself written down in the roadmap so the alignment is reproducible. It is not. (See sign-off #6.)
4. **A0 positive-control re-dock gate, ≤ 2 Å RMSD** — CLOSED — well-defined, per-strategy native references listed, abort-on-fail explicit, diagnostic plot committed.
5. **PAINS/Brenk filter between B and C** — CLOSED — §1.3, flagged-not-dropped policy is correct.
6. **Tautomer + protomer enumeration at pH 7.4 ± 0.5** — CLOSED — Dimorphite-DL + RDKit TautomerEnumerator, all states docked, best reported.
7. **Enantiomer verification for stereoactive anchors** — CLOSED — §1.4 explicit (S)-form check for raltitrexed and pemetrexed.
8. **Crystal-water handling explicit** — CLOSED — §0 documents removal + PROLIF post-hoc flag … but the PROLIF flag itself is broken; see #10 below and sign-off #2.
9. **Per-site Δ references** — CLOSED — §0 table is clean, plot 2 splits scales properly.
10. **`roc_auc_top1`, `bedroc_alpha20`, `apo_minus_holo_top1`, `pains_flag`, `lipinski_flag`, `tautomer_id`, `enantiomer` columns added** — CLOSED — §G schema lists ~30 columns; `enrichment.csv` per strategy is a clean addition.
11. **Dimer-interface residue list from 4 Å contact map** — CLOSED — §A3 script provided; runtime computation is the right call here (unlike for CIDs).
12. **Pick one box size for Strategy 1** — CLOSED — 22³ chosen, off-site-minimum risk surfaced via `pose_distance_from_box_centre` and `off_site_minimum` columns. Good mitigation.
13. **Re-verify His tautomer at pH 7.4** — CLOSED — §A1 assertion script + pdb2pqr30 fallback.
14. **Re-estimate Strategy 3 budget under HPEPDOCK** — CLOSED-but-incomplete — 12–20 h is honest in expectation but unbounded in worst case; see sign-off #3 (HPEPDOCK timeout).

## New v1 additions — audit

- **Dimorphite-DL (Ropp 2019)** — ACCEPT — correct choice for pH-dependent protomer enumeration; combines cleanly with RDKit TautomerEnumerator.
- **PLIP alongside PROLIF + `prolif_plip_agreement` column** — REVISE — scientifically sound (PLIP is the field-standard per-complex reporter, PROLIF is the heatmap matrix tool) but the install plan is one line and may not work as-written. PLIP on Python 3.14 / arm64-darwin requires OpenBabel as a transitive dependency; openbabel wheels for Python 3.14 are not on PyPI as of writing. PLIP does *not* require Java (only the public web service does); the CLI is pure Python. The roadmap should either pin Python ≤ 3.12 for the PLIP step, install OpenBabel via Homebrew (`brew install open-babel`) before `pip install plip`, or fall back to running PLIP via Docker.
- **HPEPDOCK web-service hard dependency for ≥ 6-mer peptides** — REVISE — see sign-off #3. No fallback. No per-job timeout. The Zhou 2018 service has been intermittently down in 2024–2025 and the queue can exceed 24 h. If the service is offline for the whole Strategy-3 window, the strategy stalls indefinitely and S1 ("null result") cannot fire because the issue is infrastructure, not chemistry.
- **ROC-AUC + BEDROC (α=20) + `enrichment.csv` + `fig_enrichment.png`** — ACCEPT — Truchon & Bayly 2007 α=20 is the field-standard early-enrichment metric, correctly excluded from Strategy 4 (no actives).
- **`apo_minus_holo_top1` + `cryptic_pocket_flag`** — ACCEPT — cheapest possible signal for induced-fit / cryptic-pocket effects, well-targeted.
- **`pose_distance_from_box_centre` + `off_site_minimum` columns** — ACCEPT — clean way to surface the 22³ off-site-minimum risk per row rather than forcing a re-audit of the whole panel.
- **~30-column CSV schema in §G** — ACCEPT in principle, but the `state` column has only two documented values (`apo` / `holo_raltitrexed_bound`) and Strategy 3/4 dock against `apo` only — the schema should explicitly say which `state` values are valid per strategy so a downstream pivot doesn't accidentally aggregate over an empty cell.
- **Stop Condition S4 (anchor failure = protocol suspect)** — ACCEPT — exactly the right policy; auto-diagnostic on firing is a nice touch.
- **AlphaFold-receptor sanity dock for top hit per strategy (§A.x)** — ACCEPT — adds < 2 h compute for a potentially publishable model-vs-crystal sensitivity finding.
- **Ibuprofen (CID 3672) negative control for Strategy 2** — ACCEPT.

## Section-by-section (only sections with issues)

- **§0 Crystal-water handling.** The sentence "PROLIF post-analysis (step E.1) explicitly flags whether the top pose of each Tier-1 active-site anchor would have made a water-bridged H-bond in the crystal" is **incorrect as written**. PROLIF analyses ligand–protein interactions in the *docked* complex, which by construction has *no waters* (the corrector chose "all crystallographic waters removed prior to docking"). PROLIF cannot detect a "would-have-made" water bridge; that is a *crystal-overlay annotation* task, not an interaction-fingerprint task. Concrete fix below (sign-off #2).

- **§A Pocket definition row 2 (cofactor-site).** "centroid of raltitrexed D16 (chain A, recomputed from `cofactor_A.pdbqt`)" — but Strategy 2 docks against both `apo` and `holo_raltitrexed_bound` receptors (per §D). When docking against the apo receptor, D16 is by definition absent from the receptor. The box centre must be computed from the *holo* receptor's D16 coordinates and then applied verbatim to the apo receptor's docking run. The roadmap does not say this. A fresh agent will read row 2, see "apo receptor for Strategy 2", and either (a) fail because `cofactor_A.pdbqt` is not loaded for apo, or (b) compute the centre from something else and get a different box. One-sentence fix.

- **§A0 Strategy 3 gate.** "HPEPDOCK score Δ vs scrambled-seq control ≥ 1 σ" is well-specified, but the scrambled-sequence control itself isn't defined. Is the alphabet preserved? Are cysteines kept fixed (relevant if the peptide is LSCQLYQR — randomising the Cys position changes the chemistry meaningfully)? Specify.

- **§D Strategy 3 deviation paragraph.** The peptide-vs-Vina scale-mismatch tradeoff is correctly acknowledged ("Strategy 3 is not numerically compared to Strategies 1/2/4"), but the headline plot then says "Strategy 3 on a separate subpanel (HPEPDOCK-scored, different scale)" — this is right. However, the per-strategy `enrichment.csv` for Strategy 3 uses HPEPDOCK scores as the ranking input; BEDROC and ROC-AUC are scale-invariant (rank-based), so this is fine, but the roadmap should say so explicitly to pre-empt the question.

- **§4 compute budget.** "~12–20 h" for HPEPDOCK is honest in queue-normal conditions, but the unbounded tail risk needs an explicit per-job timeout. See sign-off #3.

## New findings (issues introduced in v1)

- [HIGH] **CID verification deferred to runtime with placeholder text (BrdUMP "likely 167671", four `*(to verify A0.2)*` rows).** The R1 fix demanded a *hard-coded* literature InChIKey table; v1 has half a table. At runtime, if BrdUMP CID 135398598 returns InChIKey X and there is no literature InChIKey to compare to (because the corrector said "likely 167671"), A2 cannot fail — there is nothing to compare against. The gate becomes a no-op for those four rows.

- [HIGH] **PROLIF cannot flag a *missing* water bridge** — see §0 above. The roadmap conflates interaction-fingerprint analysis (what PROLIF does) with crystal-overlay annotation (what is actually needed). Proposed concrete fix: write a separate post-analysis script (call it `E1b_water_bridge_check.py`) that overlays the top docked pose on 1HVY (Kabsch-aligned on receptor Cα), and for every Tier-1 active-site anchor flags `water_bridge_lost = True` iff (i) at least one heavy atom of the docked pose lies within 3.5 Å of crystal Tyr258 OH, AND (ii) the crystal water O at the bridge position lies within 3.5 Å of both Tyr258 OH and the *dUMP O4 equivalent atom* of the docked pose. Only then is the bridge truly "lost" — i.e. only then did the pose displace a water that was making a real H-bond. This is a 30-line MDAnalysis script, not a PROLIF column.

- [HIGH] **HPEPDOCK has no fallback, no per-job timeout, no offline detection.** If the Zhou 2018 web service queues or returns 503 for the full Strategy-3 window, the strategy stalls and the operator cannot tell whether to wait, retry, or abort. The 50–80 h budget becomes unbounded.

- [MEDIUM] **Cofactor-site box centre ambiguity for apo docking** (§A row 2; see above).

- [MEDIUM] **PLIP install plan does not address the OpenBabel transitive dependency** on arm64-darwin Python 3.14. A fresh agent following `PIP_BREAK_SYSTEM_PACKAGES=1 pip install plip` may get a wheel-build failure and no clear fallback.

- [LOW] **Scrambled-sequence control for the HPEPDOCK A0 gate is not defined** (alphabet preservation, Cys handling).

- [LOW] **`state` column valid-values per strategy not enumerated** in §G schema.

- [LOW] **Strategy 3 enrichment.csv** should explicitly say "rank-based metrics (ROC-AUC, BEDROC) are scale-invariant, so HPEPDOCK scores can be used directly" — pre-empts the obvious reader concern.

## Sign-off requirements for R2 → R3

Corrector agent must, before I sign PASS:

1. **Pre-verify all Tier-1 CIDs *on paper now*.** Every row in the §A2 table must contain a hard-coded literature InChIKey. Specifically: BrdUMP (replace "likely 167671" with the actual free-acid CID + its actual InChIKey looked up from Santi & McHenry 1972 or PubChem), FdUR (CID 5790), methotrexate (CID 126941), raltitrexed (CID 104758), pemetrexed (CID 60843), nolatrexed (CID 60198), plevitrexed/ZD9331/BGC9331 (CID 153985 or the corrected one). Remove every `*(to verify A0.2)*` literal from the table; A0.2 is then a true cross-check, not a placeholder-resolver.

2. **Rewrite the water-bridge step.** Delete the claim that PROLIF flags missing water bridges. Replace with a new step E.1b = `water_bridge_check.py` per the script sketch above (overlay top pose on 1HVY, flag `water_bridge_lost` iff the pose displaces a water that was making a real Tyr258 ↔ dUMP O4 bridge by the dual-3.5 Å criterion).

3. **Add an HPEPDOCK fallback + per-job timeout.** Concrete: (i) per-submission timeout 4 h, on expiry mark the row `null_reason = "HPEPDOCK timeout"` and continue; (ii) if > 50 % of Strategy-3 peptide submissions time out, abort Strategy 3 with `null_reason = "HPEPDOCK service unavailable"` per Stop Condition S3 (extend S3 to cover service-down as well as not-installed); (iii) name a secondary docker fallback (CABS-dock or FRODOCK) so a fresh agent has a deterministic next step rather than silent wait.

4. **One sentence in §A row 2** clarifying that the cofactor-site box centre is computed from the *holo* receptor's D16 chain-A coordinates and applied identically to the apo and holo docking runs.

5. **Pin PLIP install path** for arm64-darwin Python 3.14: either `brew install open-babel && PIP_BREAK_SYSTEM_PACKAGES=1 pip install plip`, or pin a Python ≤ 3.12 environment for PLIP only, or use the official PLIP Docker image. Pick one and write it in §F tooling.

6. **Write the Cardinale 2011 *E. coli* / *L. casei* source peptide sequence** into §A3 so the Biopython pairwise alignment is reproducible. One line.

7. **Define the scrambled-sequence control for HPEPDOCK A0 gate** (alphabet preserved, Cys position preserved or shuffled — pick one and justify).

8. **Enumerate `state` valid values per strategy** in the §G schema (S1/S2: `apo`, `holo_raltitrexed_bound`; S3/S4: `apo`).

9. **Add one sentence to §G or §D** noting that ROC-AUC and BEDROC are rank-based and therefore scale-invariant, so HPEPDOCK kcal/mol can be used directly as the ranking score for Strategy-3 enrichment without re-scaling against Vina.

If items 1–3 land cleanly, items 4–9 are LOW-cost rewordings and R3 should be PASS.

End of Round 2 review.
