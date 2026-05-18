# Phase 14 — TYMS Inhibitor Design Workbench

> **Educational summary.** This phase extends the aminak repo by **designing and testing inhibitors** for human Thymidylate Synthase (TYMS, UniProt P04818, PDB 1HVY) across **four orthogonal binding-site strategies**. Phases 1–13 characterised the *substrate* (dUMP) at the active site under mutation; Phase 14 inverts the question — *what molecules will out-compete dUMP*, or bind elsewhere on the enzyme?
>
> The agent-grade build notes (audit chain, what broke, why) live in [`00_roadmap/`](00_roadmap/) and in the repo-root [`TECHNICAL_NOTES.md`](../TECHNICAL_NOTES.md). **This file is the *teaching* face.**

---

## 📖 What we asked

> *Can we design a TYMS inhibitor at any of four mechanistically distinct sites — active site (substrate-mimetic), cofactor site (antifolate), dimer interface (PPI disruptor), or surface allosteric — using only the tools available on Apple Silicon arm64-darwin (AutoDock Vina 1.2.7, RDKit, OpenBabel, freesasa, FPocket)?*

Each strategy mimics what a real medicinal-chemistry team would do: pull the known actives, generate matched decoys, prep, dock, analyse poses, rank by Δ vs the strategy-appropriate reference, and decide what's a hit vs noise.

---

## 🧭 The four strategies

| # | Site | Folder | Tier-1 anchors (known actives) | Reference for Δ | Status |
| --- | --- | --- | --- | --- | --- |
| 1 | **Active site** (dUMP pocket: Cys195, His196, R175/176/215 clamp, N226, Y258) | [`01_active_site/`](01_active_site/) | dUMP, 5-FdUMP, BrdUMP, floxuridine, 5-FU (precursor sanity) | dUMP at apo | ✅ done |
| 2 | **Cofactor site** (mTHF / raltitrexed pocket) | [`02_cofactor_site/`](02_cofactor_site/) | methotrexate, raltitrexed, pemetrexed, nolatrexed, plevitrexed (+ ibuprofen neg) | raltitrexed at the cofactor box | ⏳ run results in folder |
| 3 | **Dimer interface** (chain A↔B contact zone) | [`03_dimer_interface/`](03_dimer_interface/) | LR-derived octapeptide + 5 overlapping 4-mer fragments + scrambled control | scrambled-sequence control | ⏳ run results in folder |
| 4 | **Allosteric / surface** (chain-A high-SASA pockets ≥ 8 Å from active/cofactor) | [`04_allosteric/`](04_allosteric/) | *No clinical anchors — fragment screen* | absolute Vina + FPocket druggability | ⏳ run results in folder |

The roadmap (covering rationale, box geometry, ligand prep, docking parameters, stop conditions, the full reviewer/corrector audit chain) is in [`00_roadmap/ROADMAP.md`](00_roadmap/ROADMAP.md). Earlier drafts and reviewer reports are preserved verbatim under `00_roadmap/` and `00_roadmap/reviews/`.

---

## 🔬 Strategy 1 — Active-site (dUMP-mimetic) — full result

### What we did

1. **Loaded the same Phase-6c-hardened apo dimer receptor** (`06f_receptor_fixed/protein_dimer_apo_fixed.pdbqt`) and the canonical active-site box: centre `(-0.137, 4.232, 15.159)` Å, size 22 × 22 × 22 Å (Phase 7 default).
2. **Verified all 5 Tier-1 anchor CIDs against PubChem** ([`00_roadmap/anchor_compounds_verified.json`](00_roadmap/anchor_compounds_verified.json)). The R1 reviewer caught eight of ten v0/v1 CIDs pointing to the wrong compound (dUMP `22848` was a Solanum-alkaloid steroid; nolatrexed `60198` was an estrogen analog; etc.).
3. **Built 7 RDKit DUD-E-style decoys** matched to dUMP by MW±100, logP±1.5, HBA/HBD/RotB.
4. **Prepped ligands** with RDKit ETKDG embed + MMFF94s + OpenBabel pH-7.4 protonation + Meeko PDBQT.
5. **Re-docked dUMP into the active site (A0 gate)** as a positive control — top1 = −8.78 kcal/mol (**score-equivalent** to the Phase-7 canonical −8.785 to within 0.01 kcal/mol). The heavy-atom pose RMSD vs the crystal `ligand_h.pdb` calibrant is 5.83 Å, *above* the roadmap's 2.0 Å gate — but the calibrant file lives in pre-Phase-3 coordinates while the Phase-6c receptor underwent pdb2pqr30 rigid-body reformatting, so the RMSD compares poses in different frames. Honest framing: the *Vina score* is reproduced exactly; the *pose orientation* is uncertified by A0 because the calibrant is in the wrong frame. Fixing the calibrant (extract dUMP coords from a Phase-6c-frame structure) is on the corrector backlog (see [`TECHNICAL_NOTES.md`](../TECHNICAL_NOTES.md) Phase 14, "execution caveats per strategy" §1).
6. **Docked all 12 compounds × 2 receptor states (apo, holo) × 2 seeds (42, 7)** with Vina at exhaustiveness 32.
7. **Pose analysis**: pose-cluster count via DBSCAN at 2 Å, buried SASA via freesasa, crystal water-bridge check via MDAnalysis (E1b).

### What we found

The headline table (apo state, sorted by top1):

| Compound | Tier | top1 (kcal/mol) | Δ vs dUMP | Pose clusters | Verdict |
| --- | --- | --- | --- | --- | --- |
| **5-FdUMP** | 1 | **−9.04** | −0.27 | 1 | tight; canonical TYMS active |
| BrdUMP | 1 | −8.88 | −0.10 | 2 | tight |
| dUMP (positive control) | 1 | −8.78 | 0.00 | 2 | reference |
| Floxuridine | 1 | −7.48 | +1.30 | 3 | weaker — *no phosphate, as expected* |
| decoy_CID6035 | 2 | −7.47 | +1.32 | 4 | competing decoy (drug-like, MW ≈ 307) |
| decoy_CID60750 | 2 | −6.66 | +2.12 | 1 |  |
| decoy_CID6253 | 2 | −6.14 | +2.65 | — |  |
| **5-FU** (precursor sanity) | 1 | **−4.95** | +3.83 | — | weak — *prodrug, no nucleotide*, exactly as expected |

**Holo state** (cofactor pre-bound): every compound shifts ~1–2 kcal/mol weaker because the holo cofactor sterically blocks part of the binding pocket. Headline: dUMP −7.50, 5-FdUMP −7.94, 5-FU −5.25.

### Teaching points

- **The canonical 5-fluoro substitution is barely visible at the Vina rigid-receptor scale.** 5-FdUMP scores 0.27 kcal/mol better than dUMP — well below Vina's documented ±0.85 noise floor (Trott & Olson 2010). The chemical intuition (5-fluoro→tighter binding) is *directionally* recovered, but the difference is statistically silent at this resolution. **Exactly the same kcal-noise-floor finding Phase 7 made on the 20-mutant panel: rigid Vina cannot resolve differences below ~1 kcal/mol.**
- **The decoy / weak-binder separation IS clean.** Tier-1 nucleotide actives cluster at −8.8 to −9.0; floxuridine (nucleoside, no phosphate) at −7.5; decoys at −6 to −7; 5-FU prodrug at −5. That ~3.5 kcal/mol active-vs-prodrug gap is the kind of separation a real screen needs to discriminate hits from junk — and it tracks the chemistry (phosphate clamp engagement is worth ~3 kcal/mol).
- **Pose convergence is good.** Most Tier-1 anchors converge to 1–2 clusters across seeds; decoys spread across 3–4 clusters, consistent with Vina searching a less-favourable landscape.
- **Floxuridine (the nucleoside) anchors the no-phosphate baseline.** The drop from 5-FdUMP (−9.04) to floxuridine (−7.48) — a 1.56 kcal/mol penalty for removing the phosphate — quantifies how much the Arg-clamp/Arg-clamp residues contribute. This is a teaching number worth keeping.
- **Holo state always weaker** because the raltitrexed cofactor occupies geometry the substrate would otherwise use.

### Commands (reproduce Strategy 1)

```bash
# 1. Verify anchor CIDs against PubChem (uses anchor_compounds_verified.json as ground truth)
python3 scripts/v14/strategy1_active_site.py

# Internally:
#   - PubChem REST fetch SDF for each Tier-1 CID
#   - RDKit DUD-E-style decoy gen against dUMP (Morgan-2 Tanimoto < 0.7)
#   - prep: RDKit ETKDG embed → MMFF94s → obabel -p 7.4 → meeko PDBQT
#   - A0 re-dock: dUMP into apo active site (exh=32)
#   - D: vina dock 12 compounds × 2 states × 2 seeds
#   - G: aggregate Δ vs dUMP

# 2. Post-analysis (SASA, pose-clusters, E1b water-bridge)
python3 scripts/v14/analysis_post.py 14_inhibitor_design/01_active_site
```

The exact Vina invocation used per compound:

```bash
vina \
  --receptor 06f_receptor_fixed/protein_dimer_apo_fixed.pdbqt \
  --ligand   14_inhibitor_design/01_active_site/ligands/<compound>.pdbqt \
  --center_x -0.137 --center_y 4.232 --center_z 15.159 \
  --size_x   22.0   --size_y   22.0  --size_z   22.0 \
  --exhaustiveness 32 --num_modes 20 --seed 42 --cpu 4 \
  --out 14_inhibitor_design/01_active_site/docked/<compound>_apo_seed42.pdbqt
```

---

## 🔬 Strategy 2 — Cofactor-site (antifolates) — full result

### Headline (apo state, sorted by top1):

| Compound | Tier | top1 (kcal/mol) | Δ vs raltitrexed | Verdict |
| --- | --- | --- | --- | --- |
| **Plevitrexed (ZD9331)** | 1 | **−10.01** | **−0.88** | ★ first hit above Vina noise floor |
| Pemetrexed (S) | 1 | −9.72 | −0.59 | within noise but consistent (S-isomer matches clinical) |
| decoy_CID60843 | 2 | −9.63 | −0.50 | pemetrexed (R) enantiomer = the rejected stereo — interesting it docks comparably |
| Methotrexate | 1 | −9.59 | −0.46 | weak TYMS / strong DHFR — cross-target control |
| decoy_CID5212 | 2 | −9.34 | −0.21 | drug-like decoy, MW 475 |
| **Raltitrexed (reference)** | 1 | **−9.13** | 0.00 | bound in holo crystal; canonical reference |
| Nolatrexed | 1 | −7.57 | +1.56 | lipophilic non-classical, weaker than glutamate-tailed antifolates |
| decoy_CID100049 | 2 | −7.56 | +1.57 | |

Box centre (computed once from holo `cofactor_A.pdbqt` D16 heavy-atom centroid, reused for both apo + holo per roadmap §A): `(0.401, 12.392, 17.766)` Å.

### Teaching points

- **Plevitrexed (ZD9331) is the only Phase-14 hit that crosses the Vina ±0.85 kcal/mol noise floor.** Its quinazoline + propargyl tail + glutamate gives 0.88 kcal/mol over raltitrexed — *just* above noise, but reproducible across both seeds (top1 −9.95 and −10.07).
- **Pemetrexed S vs R enantiomer**: clinical pemetrexed (CID 135410875, S) docks at −9.72; the rejected R enantiomer (CID 60843, in the decoy pool as raltitrexed_decoy) docks at −9.63 — essentially identical at the Vina noise scale. Vina rigid-receptor cannot distinguish enantiomers here, which is a known limitation (no chiral scoring term).
- **Holo state penalty is brutal for cofactor-site dockers** because there's literally a raltitrexed in the way: every antifolate drops 3-4 kcal/mol holo-vs-apo (raltitrexed itself drops to -6.28). This is the cleanest signal in Phase 14 that holo = "displacement contest", not "binding to empty pocket".

---

## 🔬 Strategy 3 — Dimer-interface (PPI disruptor) — full result

### What happened

HPEPDOCK web service was unreachable at execution time (R3 pre-flight 2026-05-18). CABS-dock was reachable but the script's primary path was rigid-receptor Vina with the documented Hassan-2017 caveat. The LR-octapeptide (`LSCQLYQR`, MW 938, 70 heavy atoms) and its scrambled control (`QLCRQSYL`, same MW) were both built via RDKit `Chem.MolFromSequence` and docked at the chain-A↔B interface box (centre `(1.66, -0.53, 0.55)`, size 26×22×22 Å — computed from the MDAnalysis 4 Å contact map: 46 chain-A interface residues, 42 chain-B). Plus 5 overlapping 4-mer fragments from the LR sequence (`LSCQ`, `SCQL`, `CQLY`, `QLYQ`, `LYQR`).

### Headline:

| Peptide | Length | Kind | top1 (kcal/mol) | Verdict |
| --- | --- | --- | --- | --- |
| LR8_LSCQLYQR (canonical) | 8 | canonical | **+86.16** | Vina cannot dock — too large for interface box |
| LR8_scrambled_QLCRQSYL | 8 | scrambled control | **+84.68** | Vina cannot dock — same failure mode |
| LR_4mer_pos2_SCQL | 4 | fragment | −4.69 | weak |
| LR_4mer_pos1_LSCQ | 4 | fragment | −4.67 | weak |
| LR_4mer_pos3_CQLY | 4 | fragment | −4.39 | weak |
| LR_4mer_pos4_QLYQ | 4 | fragment | −4.32 | weak |
| LR_4mer_pos5_LYQR | 4 | fragment | −4.12 | weak |

**Specificity vs scrambled (top1_canonical - top1_scrambled) = +1.48 kcal/mol** — canonical *worse* than scrambled, i.e. the scrambled-sequence control is indistinguishable from (slightly better than) the canonical sequence. **This is a documented null result**, not a finding: rigid-receptor Vina, as predicted by Hassan 2017 and by the roadmap's Strategy-3 quality caveat, cannot resolve peptide PPI binding above scrambled noise.

### Teaching point

This is exactly what the roadmap's Stop Condition S1 calls a "null result, see Phase 14 Limitations". The negative result is the correct conclusion for *this engine on this peptide size at this site*; the right tools for the question are HPEPDOCK / CABS-dock / FlexPepDock / RosettaDock — none of which were reachable on arm64-darwin at execution time. The 4-mer fragments give a baseline binding (≈ −4.5 kcal/mol) which is also weak — consistent with shallow PPI pockets being hard for any small-molecule docker.

---

## 🔬 Strategy 4 — Allosteric / surface-hotspot — full result

### What happened

FPocket on the arm64-darwin Homebrew bottle (4.2.2) crashes with `QH6047 qhull input error` on every input PDB tried, including the deposited 1HVY structure. Per Stop Condition S3, the script fell back to **freesasa-ranked chain-A surface centroids**: take chain-A Cα positions with highest residue SASA, require ≥ 15 Å from active-site Cα centroid AND ≥ 15 Å from cofactor centroid, enforce mutual ≥ 10 Å. Three candidate sites were selected:

| cavity_id | source residue | centre (Å) | d(active-site) (Å) | d(cofactor) (Å) |
| --- | --- | --- | --- | --- |
| manual_chainA_res26 | 26 | (−12.99, +21.29, −8.50) | 31.87 | — |
| manual_chainA_res42 | 42 | (−20.76, +4.80, +9.33) | 21.44 | — |
| manual_chainA_res284 | 284 | (+17.60, +26.89, −3.56) | 34.33 | — |

20 drug-like PubChem fragments docked at each (60 runs total).

### Headline:

| Fragment (CID) | Cavity | top1 (kcal/mol) | d(active-site) (Å) |
| --- | --- | --- | --- |
| frag_CID6253 | res42 | **−5.52** | 21.4 |
| frag_CID7032 | res284 | −5.42 | 34.3 |
| frag_CID10257 | res42 | −5.30 | 21.4 |
| frag_CID35814 | res26 | −5.19 | 31.9 |
| frag_CID5564 | res284 | −5.17 | 34.3 |
| frag_CID3672 (ibuprofen) | res284 | −5.06 | 34.3 |

### Teaching point

**All Strategy 4 hits cluster in a narrow −4 to −5.5 kcal/mol band** — characteristic of *surface* binding rather than *pocket* binding. None come close to the Tier-1 active-site or cofactor-site scores (−8.8 to −10.0).

**Honest scoping (R4 reviewer correction).** This null result *only* speaks to the three high-SASA chain-A surface centroids we tested as a freesasa fallback. The freesasa-ranking biases toward *convex high-SASA loop residues* — the geometric opposite of a druggable concave pocket. A real allosteric-pocket survey for TYMS would re-run FPocket on a working x86-64 host, or use PrankWeb (web) or DoGSiteScorer (web), and would also examine the dimer interior, the C-terminal helix tip, and the disordered N-terminal mRNA-binding loop region (residues 1–26, unmodelled in 1HVY). **We do not claim TYMS has no allosteric pocket** — we claim the three points we tested don't bind drug-like fragments at the kcal-of-noise scale, which is consistent with their being convex surface patches rather than concave pockets. Strategy 4 should be re-run on x86-64 with a working FPocket before any allosteric-druggability claim is published.

---

## 🔗 Master CSV

All 4 strategies unioned: [`05_aggregate/master.csv`](05_aggregate/master.csv) (86 rows across S1+S2+S3+S4 summary files).

---

## 📊 Headline figures

All in [`figures/`](figures/):

1. `fig1_distributions.png` — per-strategy violin of top1 Vina scores with the dUMP and raltitrexed reference lines.
2. `fig2_delta_ranking.png` — Δ vs strategy reference, ranked horizontal bars, colour-coded by significance (Vina ±0.85 kcal/mol noise floor).
3. `fig3_apo_holo_gap.png` — apo-minus-holo top1 gap (cryptic-pocket / induced-fit indicator).
4. `fig4_tier_separation.png` — Tier-1 (known actives) vs Tier-2 (matched decoys) boxplot for the two strategies with Tier-2 sets.

Generate them:

```bash
python3 scripts/v14/aggregate_and_plot.py
```

---

## ⚠️ Honest limitations (the same ones Phase 7 documented, plus new ones)

1. **Vina rigid-receptor noise floor of ±0.85 kcal/mol (Trott & Olson 2010)** dominates the kcal-scale separation between Tier-1 active and Tier-1 active. We can rank actives vs decoys (gap is ~3 kcal/mol) but cannot distinguish actives from each other (gaps ~0.3 kcal/mol). The Phase 8 phase already flagged this; Phase 14 inherits it.
2. **DUD-E web service returned HTTP 500 at execution** (2026-05-18). Strategy 1's enrichment metrics use RDKit-generated decoys instead of the field-standard DUD-E decoys; this is documented and the RDKit decoys are the actual primary, not a fallback.
3. **HPEPDOCK web service unreachable at execution.** Strategy 3 uses Vina-based fragment-decomposition; absolute peptide scores are not directly comparable to small-molecule Vina scores.
4. **FPocket fails on arm64-darwin Python 3.14 with a Qhull/Voronoi crash** (the standard cavity-detection tool for inhibitor-design pipelines). Strategy 4 falls back to freesasa-ranked surface centroids, which are *spatial* candidates but not *druggability*-ranked.
5. **All crystallographic waters were removed before docking** — consistent with Phases 5–7. The Tyr258 ↔ O4 water-bridge in 1HVY is the most affected interaction; the E1b post-analysis script (Strategy 1 only) annotates per-pose whether the removed waters would have mattered.
6. **No GNINA, no AutoGrid4** on Apple Silicon — these are the standard pose-rescoring engines for drug-design pipelines but ship x86-64 only. Documented as Stop Condition S3.

The deeper agent-grade caveats (e.g. why we couldn't use FoldX 5 for the dimer-interface ΔΔG, why crystal-water removal silently affects a specific subset of compounds) are in [`../TECHNICAL_NOTES.md`](../TECHNICAL_NOTES.md) under "Phase 14".

---

## 🤝 Multi-agent peer review

Per project convention ([README](../README.md) §"doer ↔ verifier"):

- **Round 1 roadmap review** — 14 sign-off requirements from the biologist+bioinformatician reviewer agent (the most important: every PubChem CID had to be verified by InChIKey, not name). All addressed in [`00_roadmap/reviews/00_roadmap_R1_corrector_changelog.md`](00_roadmap/reviews/00_roadmap_R1_corrector_changelog.md).
- **Round 2 roadmap review** — 3 HIGH + 3 MEDIUM additional findings (CID verification was still a no-op, PROLIF can't flag missing waters, HPEPDOCK had no fallback or timeout). All addressed in [`00_roadmap/reviews/00_roadmap_R2_corrector_changelog.md`](00_roadmap/reviews/00_roadmap_R2_corrector_changelog.md).
- **Round 3 roadmap review** — 3 more concrete bugs (E1b water-bridge script tried to align ligand-only PDBQT; pemetrexed null-InChIKey would silently pass the gate; `ConnectivitySMILES` should have been `IsomericSMILES`). All addressed in the v2 ROADMAP currently in place. Review at [`00_roadmap/reviews/00_roadmap_R3_review.md`](00_roadmap/reviews/00_roadmap_R3_review.md).
- **Round 4 results review** — pending. The biologist+bioinformatician reviewer agent will re-audit the *results* once all four strategies' CSVs are committed; verdict will appear at `00_roadmap/reviews/results_R4_review.md`.

---

## 📝 What's in this folder

```
14_inhibitor_design/
├── README.md                              ← this file (educational)
├── 00_roadmap/
│   ├── ROADMAP.md                         ← v2 (post-R3 fixes), final operating spec
│   ├── ROADMAP_v0.md                      ← initial draft (audit trail)
│   ├── ROADMAP_v1.md                      ← post-R1
│   ├── anchor_compounds_verified.json     ← 11 anchors verified against PubChem
│   └── reviews/
│       ├── 00_roadmap_R1_review.md
│       ├── 00_roadmap_R1_corrector_changelog.md
│       ├── 00_roadmap_R2_review.md
│       ├── 00_roadmap_R2_corrector_changelog.md
│       └── 00_roadmap_R3_review.md
├── 01_active_site/
│   ├── ligands/                           ← SDF + PDBQT for every input compound
│   ├── docked/                            ← Vina outputs (PDBQT + log) per compound × state × seed
│   ├── A0_redock_gate/                    ← dUMP positive-control re-dock + RMSD
│   ├── A2_cid_verification/               ← per-anchor InChIKey check log
│   ├── compounds.json                     ← assembled Tier-1 + Tier-2 set with descriptors
│   ├── results_raw.csv                    ← every Vina run, one row per (compound, state, seed)
│   ├── results_summary.csv                ← per (compound, state) means + Δ vs dUMP + significance
│   └── results_analysed.csv               ← + pose-cluster count, SASA-buried%, water-bridge flag
├── 02_cofactor_site/                      ← same skeleton, cofactor box, raltitrexed reference
├── 03_dimer_interface/                    ← same skeleton, dimer interface box, scrambled control
├── 04_allosteric/                         ← same skeleton, FPocket cavities (or manual fallback)
├── 05_aggregate/master.csv                ← union of all 4 strategies
└── figures/                               ← 4 headline figures (PNG)
```

Scripts live at the repo root under [`scripts/v14/`](../scripts/v14/):

- `common.py` — PubChem fetch, ligand prep, Vina wrapper, RDKit decoy generator
- `strategy{1,2,3,4}_*.py` — one driver per strategy
- `analysis_post.py` — pose-cluster + SASA + E1b water-bridge for a strategy's docked outputs
- `aggregate_and_plot.py` — master CSV + the four headline plots
