# Phase 14 — Inhibitor Design Roadmap (DRAFT v0)

> **Scope.** Phases 1–13 of this repo characterised TYMS (P04818, PDB 1HVY) by docking the **natural substrate** dUMP and probing it with point mutants. Phase 14 inverts the question: **what small molecules will *out-compete* dUMP at the active site, and what other ligandable surfaces does TYMS expose?** We enumerate four orthogonal inhibitor strategies, dock each against the same canonical receptor used in Phases 5–13, and apply the same noise-floor honesty as the existing Phase 7 work (Δ Vina ≥ 0.85 kcal/mol to call a difference; Trott & Olson 2010).
>
> **This roadmap is a draft.** It is intended to be torn apart by a biologist+bioinformatician reviewer agent and patched by a corrector agent before any compute is spent. Open questions are tagged `[OPEN]`.

---

## 0. Inputs we inherit from earlier phases

| Artefact | Path | Provenance |
| --- | --- | --- |
| Apo dimer receptor (PDBQT, AMBER ff14SB charges, Phase-6c hardened) | `06f_receptor_fixed/protein_dimer_apo_fixed.pdbqt` | Phase 6c — see TECHNICAL_NOTES §"Phase 6c" |
| Holo cofactor A (raltitrexed, −1.82 e per copy) | `06f_receptor_fixed/cofactor_A.pdbqt` | Phase 6c |
| Holo cofactor B | `06f_receptor_fixed/cofactor_B.pdbqt` | Phase 6c |
| Reference substrate dUMP (PDBQT) | `05d_ligand_v4/dump.pdbqt` | Phase 5d |
| Active-site box centre | (−0.137, +4.232, +15.159) Å | TECHNICAL_NOTES §"Where the active-site box lives" |
| Active-site box size | 22 × 22 × 22 Å (v5 canonical), 18³ (Phase 7) | same |
| Per-residue SASA (apo & holo) | `12_phase7/04_sasa/` | Phase 7 |
| Per-residue conservation (JS) | `01b_msa_v2/conservation_scores.csv` | Phase 1b |
| Vina noise floor for Δ-call | ±0.85 kcal/mol | Trott & Olson 2010 |

**Convention reused (from Phase 4):** *more negative Vina = better binding*; for inhibitor *ranking* against dUMP we report `Δ_vs_dUMP = inhibitor_score − dUMP_score`; **negative = the inhibitor binds tighter than dUMP** at the same site.

---

## 1. Strategies (run all four, one folder each)

| # | Folder | Site | Primary literature anchor | Prior-art chemotype to mine |
| --- | --- | --- | --- | --- |
| 1 | `01_active_site/` | dUMP pocket (catalytic Cys195, His196, Arg175/176/215 phosphate clamp, Asn226, Tyr258) | Carreras & Santi 1995 (TYMS mechanism); Costi 2002 (TYMS inhibitor review) | 5-FdUMP, 5-FU, BrdUMP, ZD9331 (nolatrexed), AG-337 (thymitaq) — substrate-mimetic nucleobases & nucleotides |
| 2 | `02_cofactor_site/` | mTHF / raltitrexed pocket (D16 site, Phe80, Trp80-loop) | Costi 2002; Berger & Berger 2004 | Methotrexate, raltitrexed, pemetrexed, BGC9331, plevitrexed — antifolate scaffolds |
| 3 | `03_dimer_interface/` | A↔B subunit interface (residues identified from chain A↔B contact map) | Cardinale 2011 (TYMS dimer-disruptor LR peptide); Salo-Ahen 2015 | Octapeptide LR; small-molecule mimics of LR (none clinical) |
| 4 | `04_allosteric/` | mRNA-binding loop (residues 1–27 of TYMS; autoregulatory) | Chu 1991 (TYMS autoregulation); Brunn 2014 | No clinical ligand exists — this strategy is **exploratory**, framed as a fragment-screen with FTMap-style hotspot detection rather than scaffold-based design. |

**Why these four and not others?** All four are *mechanistically distinct* (orthogonal site usage), and the first three each have validated chemotype precedent so docking-rank-against-known-actives is a meaningful internal control. The fourth is signposted as a fragment exploration so a flat Vina landscape is the correct teaching result rather than a failure.

---

## 2. Per-strategy pipeline (the seven steps each strategy runs)

Every strategy follows the **same seven-step pipeline** so cross-strategy comparison is apples-to-apples.

```
A) Pocket definition             — coordinates + residue list + box
B) Compound set assembly         — known actives + matched-decoy expansion
C) Ligand prep                   — protonation, charges, PDBQT
D) Docking (Vina, exh=32, multi-seed)
E) Pose analysis                 — interaction profile, pose RMSD vs known crystal where available
F) SASA / pocket-occupancy / PPI metrics (strategy-dependent)
G) Aggregate scoring + Δ_vs_dUMP table
```

### A — Pocket definition

| Strategy | Box centre (Å) | Box size (Å) | Reference atoms used |
| --- | --- | --- | --- |
| 1. Active-site | (−0.137, +4.232, +15.159) | 22 × 22 × 22 | Cα centroid of residues [80, 87, 109, 135, 175, 176, 195, 196, 214, 215, 217, 218, 221, 225, 226, 258] (chain A — same as Phase 5–7 canonical) |
| 2. Cofactor-site | centroid of raltitrexed D16 (chain A, recompute from `cofactor_A.pdbqt`) | 22 × 22 × 22 | All non-H atoms of D16:A |
| 3. Dimer interface | midpoint between Cα centroids of chain A and chain B interface residues | 26 × 22 × 22 (longer along the A↔B axis) | Interface residues identified via 4 Å contact map between chains (PRODY / MDAnalysis) |
| 4. Allosteric (TS-mRNA loop hotspot) | FTMap-style consensus hotspot on chain A surface, **excluding** the active-site and cofactor-site shells (≥ 8 Å) | 20 × 20 × 20 | Top FPocket / FTMap probe cluster, ranked by SASA-weighted druggability |

`[OPEN-A1]` Reviewer: confirm the dimer-interface residue list once we compute the 4 Å contact map; the residue list below is a literature-prior, not a structure-derived list.

Literature-prior interface residues (Cardinale 2011, Salo-Ahen 2015): chain A side {30, 76, 77, 78, 81, 96, 99, 174, 175, 176, 178, 179, 199, 202}; chain B side: the same residues (homodimer, symmetric).

### B — Compound set assembly

Each strategy gets a **two-tier** compound set:

| Tier | What | Size | Source |
| --- | --- | --- | --- |
| Tier 1 (known actives) | Clinically/literature-validated binders at that site | 5–10 per strategy | PubChem CIDs hard-coded in the strategy script; SDF fetched via PubChem REST |
| Tier 2 (matched decoys / analogs) | Topology-matched property-matched decoys for ranking-test, plus near-neighbour analogs from PubChem similarity search (Tanimoto ≥ 0.7 on Morgan-2 of a Tier-1 anchor) | 20–40 per strategy | PubChem similarity REST + DUD-E-style decoy generation (RDKit; matched MW, logP, HBA, HBD, RotB) |

Tier-1 anchors (proposal, subject to reviewer veto):

| Strategy | Tier-1 anchors (PubChem CIDs) |
| --- | --- |
| 1. Active-site | 5-FU (3385), 5-FdUMP (15718), 5-FdUR (5790), BrdUMP (135398598), AG-337 / nolatrexed (60198), ZD9331 (153985) — six anchors |
| 2. Cofactor-site | Methotrexate (126941), raltitrexed (104758), pemetrexed (135410875), BGC9331 / plevitrexed (122478), CB300638 (468595) — five anchors |
| 3. Dimer-interface | LR octapeptide (LSCQLYQR; build from sequence via RDKit `Chem.MolFromSequence`), peptidomimetic analogs from Cardinale 2011 — three peptide anchors |
| 4. Allosteric | No literature anchors. Use the **ZINC15 fragment library** (a 200-compound MW ≤ 250 Da, logP ≤ 3.5 random subset) as Tier-2; Tier-1 is empty by design and that absence is the teaching point. |

`[OPEN-B1]` Reviewer: PubChem CIDs above need verification — agent should not assume hard-coded CIDs are correct. The corrector must add an assertion step that fetches the SDF, parses the InChIKey, and cross-checks against the literature name.

`[OPEN-B2]` Reviewer: are LR peptide and peptidomimetics really tractable in Vina? AutoDock Vina's accuracy on peptides ≥ 6 residues is poor (Hassan 2017). Should we either (a) restrict to peptide fragments ≤ 5 residues, (b) use HPEPDOCK / HADDOCK instead, or (c) accept that PPI-disruptor Vina scores are *qualitative*?

### C — Ligand prep

```
for each compound (SDF from PubChem):
  RDKit  — sanitize, embed with ETKDGv3, optimise with MMFF94s (max 1000 iter)
  obabel — protonate at pH 7.4, assign Gasteiger charges
  meeko  — mk_prepare_ligand.py → PDBQT with explicit rotatable bonds
  Assertion: |total charge − formal charge| < 0.05 e  (gate-fail else)
```

`[OPEN-C1]` Reviewer: meeko's Gasteiger sometimes silently re-protonates differently from the obabel step upstream. Add a separate column in the results CSV: `formal_charge_input`, `total_q_pdbqt`, `delta`, and surface every row with `delta > 0.1`.

### D — Docking

**Same canonical Vina invocation as Phase 5–7**, so results are scale-comparable:

```bash
vina \
  --receptor 06f_receptor_fixed/protein_dimer_{apo,holo}_fixed.pdbqt \
  --ligand   14_inhibitor_design/<strategy>/ligands/<compound>.pdbqt \
  --center_x <cx> --center_y <cy> --center_z <cz> \
  --size_x   <sx> --size_y   <sy> --size_z   <sz> \
  --exhaustiveness 32 \
  --num_modes 20 \
  --seed <s> \
  --cpu 4 \
  --out  14_inhibitor_design/<strategy>/docked/<compound>_<state>_seed<s>.pdbqt \
  > 14_inhibitor_design/<strategy>/logs/<compound>_<state>_seed<s>.log
```

- **Seeds:** {42, 7, 13, 99, 256} — five replicates per compound × state (apo/holo) — matches Phase 7 multi-replica convention.
- **State:** strategies 1 & 2 dock against both **apo** and **holo** receptors (so dUMP-displacement vs cofactor-displacement can be inferred from the Δ between states); strategies 3 & 4 dock against **apo only** (cofactor is irrelevant for dimer-interface and an allosteric surface).
- **Aggregation:** report `top1`, `top3_mean`, `top5_mean`, `mean_over_seeds(top1)`, `sd_over_seeds(top1)` per (compound, state).

`[OPEN-D1]` Reviewer: Phase 8b already showed that flex-residue Vina costs >30 min/compound. Should Phase 14 stay rigid-receptor (and inherit the Phase 7 null-result honesty), or run flex on the top-5-by-rigid for each strategy? Proposal: rigid-first for the full panel; flex-residue follow-up only on compounds with rigid `top1 < dUMP_top1 − 0.85`.

### E — Pose analysis

For every compound that beats dUMP at rigid Vina by Δ ≥ 0.85 kcal/mol:

1. **Interaction profile** — `prolif` (RDKit-based; already a `requirements.txt` dep upstream) computes per-residue interaction fingerprint (HBond donor/acceptor, hydrophobic, π-stacking, salt-bridge, halogen). Output: per-compound `prolif_fingerprint.json` + a heatmap of interaction-type × residue, with the dUMP fingerprint as the reference column.
2. **Pose RMSD vs the crystal pose** (where a crystal exists): use the chain-A copy of the bound ligand in 1HVY (dUMP), 1HW3 (raltitrexed), 1JU6 (5-FdUMP), etc. RMSD computed via RDKit's symmetry-corrected `GetBestRMS`. Report only on compounds whose Tier-1 anchor has a matched PDB.
3. **Pose-cluster diversity** — DBSCAN on heavy-atom RMSD across the 5 seeds × top-3 modes. Report number of distinct clusters (>2 Å); a *single* tight cluster across seeds = high-confidence pose, *many* loose clusters = Vina searching, not converging.

### F — Strategy-specific structural metrics

| Strategy | Metric | Tool | Why |
| --- | --- | --- | --- |
| 1. Active-site | Pocket-occupancy ratio (ligand SASA buried / ligand SASA free) | `freesasa` (Python) | A genuine substrate-mimetic should bury ≥ 80 % of its surface, as dUMP does. |
| 2. Cofactor-site | ΔSASA on chain A residues 80, 130–134, 218–225 (mTHF binding face) | `freesasa` Δ between apo and complex | Antifolate occupancy should specifically de-expose these residues. |
| 3. Dimer interface | **PPI metrics**: BSA (buried surface area) on dimer interface residues, ΔΔG_bind (PRODIGY-LIG) for dimer ± ligand, FoldX `AnalyseComplex` interface energy | `freesasa`, PRODIGY-LIG (online), FoldX 5 | A real PPI disruptor should *reduce* dimer interface BSA when bound. |
| 4. Allosteric | FTMap-style hotspot recapitulation: does the ligand land in any of the top 5 FPocket cavities? | `fpocket` (Homebrew) | An allosteric hit should land in a druggable pocket distinct from active/cofactor. Also report distance from active-site Cα centroid. |

`[OPEN-F1]` Reviewer: PRODIGY-LIG is a web service (no API key needed, but rate-limited). FoldX 5 requires a non-commercial licence and arm64-darwin binary availability is uncertain. Fallback: report BSA + interface residue-residue contacts from MDAnalysis, document the missing ΔΔG, do not silently skip.

`[OPEN-F2]` Reviewer: `fpocket` may not be installed. Add an install-or-skip gate.

### G — Aggregate scoring + Δ_vs_dUMP

Per-strategy `results.csv` schema:

```
compound_name, pubchem_cid, inchikey, mw, logp, hba, hbd, rotb, tpsa,
tier, anchor_for, state, seed,
top1_kcalmol, top3_mean, top5_mean, pose_cluster_count, pose_rmsd_vs_xtal,
sasa_buried, sasa_buried_pct, ppi_bsa, ppi_delta_bsa, prolif_overlap_with_dump,
delta_vs_dump_top1, delta_vs_dump_significant_p085,
notes
```

Final `14_inhibitor_design/05_aggregate/master.csv` is the union of the four per-strategy CSVs plus a `strategy` column.

Headline plots (saved under `figures/`):

1. Per-strategy distribution of `top1_kcalmol` (violin × 4 strategies), with dUMP reference line — `fig_distributions.png`
2. `delta_vs_dump_top1` ranked horizontal bar, colour-coded by Tier and significance gate — `fig_delta_ranking.png`
3. PROLIF interaction-fingerprint heatmaps (one per strategy's top hit + dUMP) — `fig_prolif_heatmaps.png`
4. SASA / BSA / pocket-occupancy panel (strategy-specific) — `fig_structural_metrics.png`
5. Per-strategy pose overlay PyMOL render (top hit vs dUMP vs Tier-1 anchor) — `fig_pose_overlays.png`

---

## 3. Reviewer / corrector loop

> Exactly per user instruction: one biologist+bioinformatician **reviewer agent**, one **corrector agent**, unbounded loop until reviewer signs `PASS` on every section.

Each cycle:

1. **Reviewer agent** receives the current roadmap and writes `reviews/00_roadmap_R<n>_review.md` with verdicts `PASS / CONDITIONAL_PASS / FAIL` per section and per `[OPEN-*]` tag.
2. **Corrector agent** receives the review and writes the next roadmap version (`ROADMAP.md` overwritten; previous saved as `ROADMAP_v<n-1>.md`).
3. Repeat until the reviewer's overall verdict is `PASS` (no FAIL, no CONDITIONAL_PASS, every `[OPEN-*]` tag resolved or explicitly converted to a "documented limitation").

Same loop runs **again** after each strategy is executed (review the *results*, not the plan): cycle continues until reviewer signs the analysis off.

---

## 4. Compute budget and stop conditions

| Stage | Estimated wall-time (single-Mac, M-series, exh=32) | Comments |
| --- | --- | --- |
| Compound assembly + prep (all 4 strategies) | ~20 min | PubChem REST + RDKit/meeko |
| Strategy 1 docking (5 anchors + 25 analogs) × 2 states × 5 seeds | ~5–8 h | 300 Vina runs |
| Strategy 2 docking (5 anchors + 25 analogs) × 2 states × 5 seeds | ~5–8 h | 300 Vina runs |
| Strategy 3 docking (3 peptide anchors + 10 mimics) × 1 state × 5 seeds | ~3–5 h | peptides slower in Vina |
| Strategy 4 docking (200 fragments) × 1 state × 1 seed (fragment screen) | ~6–10 h | small ligands fast |
| Pose analysis + figures | ~1–2 h | |
| **Total** | **~20–35 h sequential** | parallelisable across strategies if the user-budget permits a 4-way fork |

Stop conditions (we honour these *up-front* so the unbounded reviewer loop doesn't churn on irreducible compute issues):

- **S1.** Rigid-receptor Vina cannot resolve binding-mode differences below Vina's ±0.85 kcal/mol noise floor. If a strategy yields no compound with `Δ_vs_dUMP ≤ −0.85`, the reviewer must accept "null result, see Phase 14 Limitations" rather than loop indefinitely on the same compute.
- **S2.** Phase 14 inherits Phase 8b's documented flex-residue compute budget (30 min per compound). Flex-only-on-rigid-hits gate (D1 above) is mandatory.
- **S3.** Tools that aren't installed (FoldX, PRODIGY-LIG, GNINA, autogrid4) are *documented as missing* and the corresponding column in the results CSV is left blank with a `null_reason` value — not silently skipped.

---

## 5. Deliverables (what gets committed back)

- `14_inhibitor_design/00_roadmap/ROADMAP.md` — this file, final version after reviewer PASS.
- `14_inhibitor_design/00_roadmap/reviews/` — all review iterations verbatim.
- `14_inhibitor_design/{01..04}_<strategy>/`:
  - `ligands/` — all input SDF + PDBQT
  - `docked/` — all output PDBQT + logs
  - `analysis/` — PROLIF JSONs, SASA tables, PPI metrics, pose-RMSD tables
  - `results.csv` — strategy-level master
  - `README.md` — strategy-level educational write-up
- `14_inhibitor_design/05_aggregate/master.csv` — cross-strategy union
- `14_inhibitor_design/figures/` — the five headline plots
- **Root `README.md` update** — educational summary of Phase 14 + commands + figures, *not* a think-pad.
- **Root `TECHNICAL_NOTES.md` update** — Phase 14 caveats, what-broke-and-why, agent-grade detail.
- **Root `CHANGELOG.md` update** — commit-level summary.

---

## 6. Open questions for the reviewer (tag index)

- `[OPEN-A1]` Confirm or veto the literature-prior dimer-interface residue list.
- `[OPEN-B1]` Add a CID-to-InChIKey verification gate.
- `[OPEN-B2]` Peptides in Vina — restrict, switch tool, or document?
- `[OPEN-C1]` Charge-delta surface in CSV.
- `[OPEN-D1]` Flex-residue follow-up gating threshold.
- `[OPEN-F1]` PRODIGY-LIG / FoldX availability + fallback.
- `[OPEN-F2]` `fpocket` install-or-skip.

End of draft v0.
