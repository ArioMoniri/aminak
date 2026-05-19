# PPTX v3 + README review

## Overall verdict
CONDITIONAL_PASS

Most claims check out against the underlying CSV/JSON, the BSA and Plasmodium-21 numbers are real, and the WT −8.55 / −8.78 reconciliation is sound. Two issues block ship: (1) the Modeller-vs-AlphaFold figure is **broken** (every bar reads 0.00 Å and a placeholder warns "Lovell %favoured not in Phase-7c summary"), and (2) the labels on the three new PyMOL renders overlap to the point of being unreadable. Both are fixable without redesign.

## Findings (severity HIGH / MEDIUH / LOW)

1. **[HIGH] `modeller_vs_alphafold.png` is non-functional.** All three rows (Modeller B99990010, B99990003, AlphaFold v6) show "0.00 Å" bars of zero length, and the right margin text says "Lovell %favoured not in Phase-7c summary". The plotting script did not load the actual columns from `modeller_vs_alphafold.csv` (which holds RMSD ≈ 0.37–0.39 Å and %favoured 94.53–95.44). **Fix:** re-render reading `rmsd_vs_1HVY_super_A` and `pct_favoured`; plot Lovell as a second panel. The slide-10 caption claim "2 of 3 beat the crystal on Lovell %favoured" is actually *under-counted* — the CSV shows **3 of 3** (AlphaFold 94.53, Modeller B99990003 95.44, Modeller B99990010 95.09, all > 1HVY 92.2). Update caption to "all 3 beat the crystal" or "2 of 3 Modeller models + AlphaFold" depending on which framing is intended.

2. **[HIGH] Slide-5 PyMOL renders have overlapping residue labels.**
   - `dimer_overview_clean.png`: central labels collide into an illegible blob ("ArgAA-AlAA295"-shaped overstrike near the active site).
   - `dimer_activesite_clean.png`: R175 / R176 fused on the right; multiple labels stacked at left and center; can't tell which sidechain belongs to which name.
   - `cavity18_carve_clean.png`: `GLN203` / `PHE201` / `ASN201` overlap inside the carve-out.
   **Fix:** in PyMOL set `set label_position, [3,3,3]`, use `label_size 14`, drop redundant duplicated labels (currently the same residue is labelled on both protomers), or anchor labels with `label_outline_color white` for readability.

3. **[MEDIUM] R215 is not a "top dimer hot-spot".** Slide 11 says "R175 / R176 / R215 are also the phosphate-clamp residues — substrate + dimer share the same hot-spots". `ppi_summary.json` confirms R175 #1 (184.5 Å²) and R176 #3 (137.9 Å²), but **R215 is #10** (66.9 Å²) — mid-pack, ~36% of R175's BSA. The narrative is defensible but the word "top" overstates. **Fix:** soften to "R175 / R176 are top-3 dimer hot-spots, and R215 also contributes (66.9 Å²)", or list P59 (#2, 156.7 Å²) and R202 (#4, 116.3 Å²) which actually rank above R215.

4. **[MEDIUM] Slide 7 "+0.77 largest holo Δ Vina" is the v5 panel number, not Phase 7 multi-replica.** Slide header says "PHASE 7 · MUTAGENESIS" but the +0.77 figure comes from the 20-mutant v5 panel (README line 27). In the Phase 7 multi-replica CSV, R215A_N226A_holo = −7.4886 vs WT_holo = −7.4864, giving Δ = +0.002. Both panels make the same point (below noise) but the slide attributes the v5 number to Phase 7. **Fix:** caption "Phase 5–7" or annotate "v5 panel (canonical seed 42)" so a reader cross-referencing the multi-replica CSV doesn't see a 100× mismatch.

5. **[LOW] Slide 5 deep-link `06f_receptor_fixed/protein_dimer_apo_fixed.pdbqt`.** The file exists at `/Users/ario/Downloads/aminak-inhibitor/06f_receptor_fixed/protein_dimer_apo_fixed.pdbqt` (repo root), NOT under `14_inhibitor_design/`. If the PPTX hyperlink is relative to the presentation folder it will 404; if it is relative to repo root it resolves. **Fix:** confirm hyperlink anchoring; safest is to make it `../../06f_receptor_fixed/...` from the slide.

6. **[LOW] Slide 9 says "1HVY crystal: 92.2 %"** which is consistent with the CSV-implied baseline; no action needed, just flag that the figure on slide 10 needs to surface this comparison line.

7. **[LOW] Slide 6 prose "Holo state (cofactor present): −7.49 kcal/mol".** Matches CSV (WT_holo = −7.4864). Good.

## Spot-check verdicts

- **WT −8.55 / −8.78 consistency: PASS.** `multi_replica_aggregate.csv` row WT_apo: mean = −8.5512, n_seeds_ok = 5, spread 0.113, SD 0.050 — the slide-6 number "−8.55 ± 0.05" is exact. The seed-42 row inside multi-replica is −8.586 (18 Å box), distinct from the Phase-6 canonical −8.785 (22 Å box, README line 658) — the README explains the box-size difference, slide 6 carries the −8.78 footnote correctly.
- **PPI BSA = 2079 Å² reasonable: PASS.** Within the typical 1500–3000 Å² band for a stable obligate homodimer; consistent with TYMS literature (~2000 Å² per side). The 4158 Å² total = 2 × 2079 arithmetic is correct.
- **R175 / R176 / R215 top hot-spots claim: PARTIAL.** R175 #1, R176 #3 confirmed; R215 ranks #10 (see finding #3). Not a top hot-spot.
- **Plasmodium 21 / 36 divergence: PASS.** Re-counted `cavity18_residues.csv`: 36 rows total, `aa_Plasmodium_falciparum != human_aa` (excluding "−" gaps) = **21**. Slide 20 and README line 912 are consistent. One row (B25) has Pf as "−" (gap); the 21 count correctly excludes it.
- **Modeller-vs-AlphaFold "2 of 3 beat crystal": FAIL on figure, INCORRECT on claim.** Figure broken (see finding #1). True count is 3 of 3 beat 1HVY's 92.2 % (AF 94.53, M003 95.44, M010 95.09).
- **Hyperlink targets exist: MOSTLY PASS.** Verified: `04_allosteric/cavity18_evidence/viewers/cavity18_indazole.html` ✓, `…ibuprofen.html` ✓, `06_smina_rescore/` ✓, `07_advanced_methods/haddock3/` ✓ (with README + config + active-residues), `07_advanced_methods/mmgbsa/` ✓ (with README + tleap.in + mmpbsa.in). The `06f_receptor_fixed/protein_dimer_apo_fixed.pdbqt` target lives at repo root, not under 14_inhibitor_design — relative-path risk (finding #5).
- **New images readable: FAIL on three of six.** `ppi_dimer_interface.png` (clean bar chart) and `conservation_multipanel.png` (3-panel JS / tree / clamp table) are publication-quality. `modeller_vs_alphafold.png` is empty (finding #1). `dimer_overview_clean.png`, `dimer_activesite_clean.png`, and `cavity18_carve_clean.png` have colliding residue labels (finding #2).

## Ready-to-ship verdict
NO — fix `modeller_vs_alphafold.png` (regenerate from the real CSV columns, and update the caption to "3 of 3 beat the crystal") and re-render the three PyMOL panels with non-overlapping labels; the other findings are wording tweaks that can ride along.
