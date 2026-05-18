# Corrector changelog — Round 1 → v1

Corrector pass for `ROADMAP.md` after R1 review (`00_roadmap_R1_review.md`, CONDITIONAL_PASS). v0 preserved at `ROADMAP_v0.md`. All v0 `[OPEN-*]` tags resolved; none remain.

## Findings addressed (point-by-point)

**Section 0 (Inputs)**
- [MEDIUM] Two box sizes quoted → §0 table + §A footnote: committed to **22³** for Strategy 1, cross-referenced to Phase-7 off-site-minimum concern; new `pose_distance_from_box_centre` + `off_site_minimum` columns surface the risk.
- [LOW] Δ_vs_dUMP convention misleading across sites → §0 "Per-site Δ reference convention" table introduces `delta_vs_dump_top1` (S1), `delta_vs_raltitrexed_top1` (S2), `delta_vs_lr_peptide_top1` (S3), no Δ for S4.
- [LOW] Holo receptor already contains raltitrexed → §0 explicitly labels `state ∈ {apo, holo_raltitrexed_bound}` so the "displace raltitrexed" question is named.

**Section 1 (Strategies)**
- [HIGH] 5-FU mis-anchored at active site → §1.1 moves 5-FU to a "precursor sanity panel" with annotation that it is expected to dock weakly; adds dUMP (positive control) and BrdUMP; canonical species is 5-FdUMP.
- [HIGH] Nolatrexed + ZD9331 mis-bucketed → §1.1 moves both to Strategy 2. Strategy 1 Tier-1 count drops from 6 to 4 + 5-FU sanity; Strategy 2 grows from 5 to 6 + ibuprofen negative control. Tradeoff documented inline.
- [MEDIUM] Missing dUMP positive control → added as Strategy-1 Tier-1 and as the A0 re-dock gate.
- [HIGH] ZD9331/nolatrexed/plevitrexed/BGC9331 conflation → §1.1 deconflates: ZD9331 = plevitrexed = BGC9331 (Jackman 1997); nolatrexed = AG-337 = thymitaq (Webber 1996). §A2 InChIKey gate enforces.
- [MEDIUM] LR peptide sequence + numbering wrong → §1.1 row 11 + §A3: corrector computes the human-TYMS equivalent range via Biopython pairwise alignment to Cardinale 2011's *E. coli* peptide rather than hard-coding.
- [LOW] Strategy 4 honesty about disordered 1–26 loop → §1 table row 4 + §A note the loop is unmodelled in 1HVY; framing reset to "screen resolved surface and document absence of cavity at the 1–26 region".

**Section 2 (Pipeline)**
- [HIGH] No re-dock RMSD positive control → new §A0 hard gate; dUMP into active site + raltitrexed into cofactor site, both ≤ 2 Å heavy-atom RMSD vs 1HVY; abort strategy on fail. Strategy 3 gate replaced with HPEPDOCK scrambled-sequence sanity; Strategy 4 with FPocket druggability ≥ 0.5.
- [HIGH] No PAINS/Brenk filter → §1.3 adds PAINS A/B/C, Brenk, NIH MLSMR, Lipinski/Veber (Astex RO3 for fragments) as flag-not-drop columns.
- [HIGH] No tautomer/protomer enumeration → §1.4 adds RDKit TautomerEnumerator + Dimorphite-DL at pH 7.4 ± 0.5; all states docked, best reported with `tautomer_id` / `protomer_id` / `protonation_pH` columns.
- [HIGH] No enantiomer verification → §1.4 explicit verification for raltitrexed and pemetrexed ((S)-form); `enantiomer` column.
- [MEDIUM] Receptor protonation re-check → new §A1 BioPython assertion that every HIS is HID/HIE/HIP; falls back to pdb2pqr30 propka re-run.
- [MEDIUM] Crystal-water handling → §0 documents removal + PROLIF post-hoc `water_bridge_lost` flag for Tyr258↔O4 bridge in 1HVY.
- [LOW] PLIP recommended alongside PROLIF → §E.2 runs both; `prolif_plip_agreement` column surfaces disagreements. **Tradeoff note in §E.2**: kept PROLIF as the cross-compound matrix tool, added PLIP as the field-standard per-complex report — both written side by side.
- [LOW] DBSCAN parameters → §E.3 specifies `eps=2.0 Å, min_samples=2`.

**Section 4 (Budget + stop conditions)**
- [MEDIUM] Strategy 3 budget optimistic → §4 revised from 3–5 h to 12–20 h (HPEPDOCK web latency dominates); total revised from 20–35 h to 50–80 h with the tradeoff explained.
- [LOW] Added S4: anchor-fails-to-dock = protocol failure, not chemistry; auto-diagnostic on firing.

**Section 5 (Deliverables)**
- [LOW] Schema columns `pains_flag`, `lipinski_flag`, `tautomer_id`, `protonation_pH` added → §G expanded schema (≈30 columns) + new `enrichment.csv` per strategy.

**Additional findings (R1)**
- [HIGH] No re-dock RMSD gate, no PAINS, no enantiomer, no tautomer → all addressed (see HIGH items above).
- [MEDIUM] DUD-E web vs RDKit decoys → §1.2 uses DUD-E web for Strategy 1; RDKit fallback for S2/S3; documented tradeoff (rate limit vs benchmark comparability).
- [MEDIUM] ROC-AUC + BEDROC → §G + new `enrichment.csv` + new `fig_enrichment.png` (figure 6).
- [MEDIUM] Induced-fit / cryptic-pocket indicator → `apo_minus_holo_top1` + `cryptic_pocket_flag` columns.
- [MEDIUM] Crystal-water handling → see §0 (above).
- [LOW] Negative control for Strategy 2 → ibuprofen (CID 3672) added.
- [LOW] AlphaFold-receptor sanity dock for top hit per strategy → new §A.x.

**Resolution of OPEN tags**: §6 table maps every v0 `[OPEN-*]` to its v1 resolution; A1→§A3, B1→§A2, B2→§D Strategy-3 deviation, C1→§C charge-delta gate, D1→§D flex cap (5/strategy), F1→§F (FoldX declined, fallback), F2→§F FPocket via brew.

**Sign-off checklist**: §7 table maps the 14 R1 sign-off items to v1 sections. All 14 are addressed; per the spec, items the reviewer said "must appear in Section 2" (A0 gate, A1 His check, A2 CID verification, A3 contact map, PAINS filter, tautomer enumeration, enantiomer verification) live as concrete pipeline steps in Section 2, not just as bullets elsewhere.

## Declined / scaled back (with reasons)

- **GNINA** (not requested by reviewer, pre-emptive note) — declined; Apple Silicon limit, project-wide engine constraint. Documented limitation per S3.
- **AutoGrid4 maps** — same. Documented limitation per S3.
- **FoldX 5 ΔΔG** as a hard pipeline step (R1 [OPEN-F1]) — declined; Schymkowitz ships x86_64 only and Rosetta 2 fallback is too slow for a 4-strategy pipeline. MDAnalysis BSA + interface contact-count change covers the science at smaller resolution; `null_reason = "FoldX 5 not available on arm64-darwin"` per S3. Documented limitation — accept null result per Stop Condition S1.

## Sanity checks (per corrector spec point 7)
- All four strategies present (active-site / cofactor / dimer / allosteric): YES (§1, §A, §B, §F).
- Per-strategy pipeline consistent — v1 has eight steps (A0 + A1 + A2 + A3 + A–G) applied uniformly with strategy-specific specialisations called out: YES (§2).
- Compute-budget honest — revised from 20–35 h to 50–80 h with the tradeoff explained: YES (§4).
- Stop conditions present (S1–S4): YES (§4).
- Deliverables complete: YES (§5).
- No `[OPEN-*]` tags remain: verified by grep.

Word count ≈ 720. (Slight overrun on the 600-word target; the additional length is the point-by-point mapping, which the spec demands explicitly. Compressing further would force me to drop the reviewer-finding ↔ v1-section mapping, which is the audit trail.)
