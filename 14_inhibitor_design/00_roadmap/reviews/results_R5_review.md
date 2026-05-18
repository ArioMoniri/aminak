# Phase 14 Results Review — Round 5 (verification of R4 fixes)

## Overall verdict
**PASS** (subject to two ongoing parallel fixes, see below)

All four R4 blockers are closed by inspection of the current artefacts. The corrector took the lowest-risk path on each (drop the broken column rather than reinterpret it; soften A0 prose rather than fabricate a pose; dedupe rather than re-run; add explicit scoping disclaimer rather than retract S4). No new issues introduced.

## R4 blocker status

1. **SASA column bug — CLOSED.** `head -1` on all four analysed CSVs confirms `sasa_buried_pct` is absent everywhere:
   - `01_active_site/results_analysed.csv` columns end at `water_bridge_residues` (no SASA col).
   - `02_cofactor_site`, `03_dimer_interface`, `04_allosteric` also free of the column.
   - `05_aggregate/master.csv` likewise has no `sasa_buried_pct` (verified via pandas).
   - `grep -rn "sasa_buried_pct"` across `14_inhibitor_design/` returns hits *only* in (a) historical roadmap docs (`ROADMAP*.md`, expected — spec history), and (b) the R4 review file itself. No Python script and no live README narrative reference the dropped column, so no figure/aggregator cascade.

2. **A0 frame-mismatch defence — CLOSED (prose softening; pose-correctness still open on parallel track).** README line 38 (Strategy 1, step 5) now reads: "top1 = −8.78 kcal/mol (**score-equivalent** to the Phase-7 canonical −8.785 to within 0.01 kcal/mol). The heavy-atom pose RMSD vs the crystal `ligand_h.pdb` calibrant is 5.83 Å, *above* the roadmap's 2.0 Å gate — but the calibrant file lives in pre-Phase-3 coordinates while the Phase-6c receptor underwent pdb2pqr30 rigid-body reformatting, so the RMSD compares poses in different frames. Honest framing: the *Vina score* is reproduced exactly; the *pose orientation* is uncertified by A0 because the calibrant is in the wrong frame. Fixing the calibrant … is on the corrector backlog." This is the requested "score-equivalent, pose-uncertified" language. The original overclaim ("matches Phase-7 canonical exactly") is gone.

3. **Duplicate decoy rows — CLOSED.** `01_active_site/results_analysed.csv` is 44 data rows; `drop_duplicates(subset=['compound','state','seed'])` returns 44 — zero duplicates. Each of the three failed decoys (`decoy_CID65349`, `decoy_CID6149`, `decoy_CID49846`) now appears exactly 4× (apo/holo × seeds 7/42), the expected design grid, with no triplication. Master.csv also clean (0 dups in S1 slice).

4. **S4 overstated conclusion — CLOSED.** README line 179 carries the new heading "**Honest scoping (R4 reviewer correction).**" with the exact requested points: (a) result *only* speaks to the three tested centroids, (b) freesasa fallback biases toward convex high-SASA loops vs concave druggable pockets, (c) a real survey needs FPocket on x86-64, or PrankWeb/DoGSiteScorer, plus dimer interior / C-term helix tip / disordered N-term loop, (d) explicit "We do not claim TYMS has no allosteric pocket." The conclusion is appropriately scoped without retracting the empirical observation.

## New issues introduced by R4 fixes
None observed.
- No Python aggregator/figure script references the dropped column (grep -rln `sasa_buried` over `*.py` returns nothing).
- Master.csv aggregates cleanly (86 rows, expected strategy counts: S1=16, S2=18, S3=7, S4=45).
- README internal consistency preserved — the dropped column was never quoted in narrative tables.
- All four `figN_*.png` files still present in `figures/` (not regenerated, but no dependency on dropped column means staleness risk is nil).

## Subject to ongoing fixes
- **A0 frame-aligned RMSD recomputation** (parallel track): once the calibrant is re-extracted in Phase-6c coordinates and a real frame-aligned RMSD lands, the softened "pose-uncertified" language on README line 38 should be tightened back to a concrete "frame-aligned pose RMSD = X Å" statement, with the A0 gate either passing (<2 Å) or formally documented as a real failure (not a frame artefact). Today's PASS is conditional on that distinction not flipping the verdict.
- **Strategy 4 re-run with working FPocket** (parallel track): if real druggability-ranked cavities are found and re-screened, the Strategy 4 "null" headline and the §179 disclaimer will both need to be rewritten against the new cavity set. The current scoped null is defensible *as scoped*.

## Sign-off
**PASS.** All four R4 blockers are concretely closed in the current artefacts, with no cascade damage from the fixes. Phase 14 results are publishable in their current form, with the two parallel tracks (A0 frame-aligned RMSD, S4 FPocket re-run) treated as planned upgrades rather than open blockers. If either parallel track returns a result that contradicts the current README narrative (real A0 pose still > 2 Å, or a real druggable allosteric pocket exists), re-open R6 to update the affected sections — but no action required from R5 reviewers in the meantime.
