# Phase 14 — Inhibitor Design Roadmap (v2, post-R2 corrector pass)

> **Scope.** Phases 1–13 characterised TYMS (P04818, PDB 1HVY) by docking the **natural substrate** dUMP and probing it with point mutants. Phase 14 inverts the question: **what small molecules will *out-compete* dUMP at the active site, and what other ligandable surfaces does TYMS expose?** We enumerate four orthogonal inhibitor strategies, dock each against the same canonical receptor used in Phases 5–13, and apply the same noise-floor honesty as Phase 7 (Δ Vina ≥ 0.85 kcal/mol to call a difference; Trott & Olson 2010).
>
> **Status.** Round 2 reviewer audit returned CONDITIONAL_PASS with 3 HIGH + 3 MEDIUM findings; this v2 bakes all of them in. The ten Tier-1 anchor CIDs have been **directly verified** against PubChem (table baked in, audit JSON at `anchor_compounds_verified.json`). E1b crystal water-bridge check, HPEPDOCK timeouts + CABS-dock fallback + scrambled-control, and arm64-darwin PLIP install gate are all explicit. The reviewer/corrector loop continues until reviewer signs PASS.
>
> **Engine note.** Apple Silicon constrains us to **AutoDock Vina 1.2.7** as the only docking engine. The reviewer suggested HPEPDOCK for strategy 3 peptides; that is a web service (Zhou 2018) so it remains usable. GNINA / AutoGrid4 / FoldX-x86 are not first-class on arm64-darwin and we document, not silently skip, anywhere they are absent.

---

## 0. Inputs we inherit from earlier phases

| Artefact | Path | Provenance |
| --- | --- | --- |
| Apo dimer receptor (PDBQT, AMBER ff14SB charges, Phase-6c hardened) | `06f_receptor_fixed/protein_dimer_apo_fixed.pdbqt` | Phase 6c — see TECHNICAL_NOTES §"Phase 6c" |
| Holo cofactor A (raltitrexed, −1.82 e per copy) | `06f_receptor_fixed/cofactor_A.pdbqt` | Phase 6c |
| Holo cofactor B | `06f_receptor_fixed/cofactor_B.pdbqt` | Phase 6c |
| Reference substrate dUMP (PDBQT) | `05d_ligand_v4/dump.pdbqt` | Phase 5d |
| Active-site box centre | (−0.137, +4.232, +15.159) Å | TECHNICAL_NOTES §"Where the active-site box lives" |
| **Active-site box size — single canonical for Phase 14** | **22 × 22 × 22 Å** (v5 canonical) | see §2A footnote below |
| Per-residue SASA (apo & holo) | `12_phase7/04_sasa/` | Phase 7 |
| Per-residue conservation (JS) | `01b_msa_v2/conservation_scores.csv` | Phase 1b |
| Vina noise floor for Δ-call | ±0.85 kcal/mol | Trott & Olson 2010 |

**Per-site Δ reference convention (revised after R1).** The reviewer flagged that "Δ_vs_dUMP" across sites is misleading because dUMP is not the displaced species at the cofactor site or the dimer interface. We therefore report **site-appropriate Δ references**:

| Strategy | Δ reference | Variable | Sign convention |
| --- | --- | --- | --- |
| 1. Active-site | dUMP at the active-site box | `delta_vs_dump_top1` | negative = inhibitor tighter than dUMP at same site/state |
| 2. Cofactor-site | raltitrexed at the cofactor-site box | `delta_vs_raltitrexed_top1` | negative = candidate tighter than raltitrexed at the cofactor site |
| 3. Dimer-interface | LR-derived peptide reference at the interface box | `delta_vs_lr_peptide_top1` | negative = small molecule tighter than the peptide lead at the interface |
| 4. Allosteric (fragment screen) | **no Δ reference** — report absolute Vina score and FPocket druggability instead | n/a | fragment exploration, ranking by `top1_kcalmol` directly |

The Phase-4 sign convention is preserved (more negative Vina = better binding). The cross-strategy headline plot (Section 2G figure 2) now uses each strategy's own Δ reference, not a single dUMP column.

**Holo-vs-apo labelling (R1 §0).** Docking against the holo receptor at the cofactor site asks "can the candidate displace raltitrexed?", which is *different* from "does it bind the empty cofactor pocket?". We label both questions explicitly in the `state` column (`apo` vs `holo_raltitrexed_bound`) and report both for strategies 1 and 2.

**Crystal-water handling (R1 §2).** Decision: **all crystallographic waters are removed prior to docking**, consistent with Phase 5–7. The known Tyr258 ↔ dUMP O4 bridging water in 1HVY is documented as a limitation in `TECHNICAL_NOTES` Phase-14 section. PROLIF post-analysis (step E.1) explicitly flags whether the top pose of each Tier-1 active-site anchor would have made a water-bridged H-bond in the crystal; if it would have, the row gets `water_bridge_lost = True` so the bias is visible in the CSV, not hidden.

---

## 1. Strategies (run all four, one folder each)

| # | Folder | Site | Primary literature anchor | Prior-art chemotype to mine |
| --- | --- | --- | --- | --- |
| 1 | `01_active_site/` | dUMP pocket (catalytic Cys195, His196, Arg175/176/215 phosphate clamp, Asn226, Tyr258) | Carreras & Santi 1995 (TYMS mechanism); Longley *et al.* 2003 (5-FU mechanism review) | **Nucleotide / nucleobase mimics only** — see anchor list below |
| 2 | `02_cofactor_site/` | mTHF / raltitrexed pocket (D16 site, Phe80, Trp80-loop) | Costi 2002 (TYMS inhibitor review); Berger & Berger 2004 | Classical and non-classical antifolates — nolatrexed and ZD9331 **relocated here** from Strategy 1 per R1 |
| 3 | `03_dimer_interface/` | A↔B subunit interface (residues identified from the 4 Å contact map; see §2A) | Cardinale *et al.* 2011 (FEBS J 278:1487); Salo-Ahen 2015 | LR-derived octapeptide and short peptidomimetics. Primary docker for ≥ 6-mer peptides: **HPEPDOCK** (Zhou 2018); Vina retained for ≤ 5-residue mimetics only. |
| 4 | `04_allosteric/` | Surface hotspot detection on chain A, **excluding** active-site and cofactor-site shells (≥ 8 Å). Note: residues 1–26 of TYMS are intrinsically disordered in 1HVY (no resolved electron density per the PDB header for chain A); FPocket will not find a pocket on the unmodelled mRNA-binding loop. **The honest framing of Strategy 4 is therefore "we screen the resolved surface for cryptic druggable pockets and document the absence of one at the 1–26 region."** | Chu 1991; Brunn 2014; FPocket (Le Guilloux 2009) | Exploratory fragment screen — no clinical anchors. |

**Why these four?** All four are *mechanistically distinct* (orthogonal site usage). Strategies 1, 2, and 3 each have validated chemotype precedent so docking-rank-against-known-actives is a meaningful internal control. Strategy 4 is signposted as a fragment exploration so a flat Vina landscape is the correct teaching result rather than a failure.

### 1.1 Tier-1 anchors (v2 — directly verified from PubChem)

R1 caught three serious chemistry errors in v0; R2 demanded the v1 "to verify" placeholders be replaced by *actually verified* CIDs. v2 baked-in verification was performed by the corrector via `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/...` queries on 2026-05-18; the per-compound JSON record (CID, InChIKey, formula, MW, IUPACName, rejected-CID audit list) is at `00_roadmap/anchor_compounds_verified.json` and is the ground truth.

| Strategy | Tier-1 anchor | **Verified PubChem CID** | **Verified InChIKey** | Formula | MW | Role |
| --- | --- | --- | --- | --- | --- | --- |
| 1. Active-site | dUMP *(positive control, re-dock gate)* | **65063** | `JSRLJPSBLDHEIO-SHYZEUOFSA-N` | C9H13N2O8P | 308.18 | Substrate; mandatory re-dock RMSD calibrant (1HVY co-crystal, Phan 2001) |
| 1. Active-site | 5-FdUMP *(canonical active species)* | **8642** | `HFEKDTCAMMOLQP-RRKCRQDMSA-N` | C9H12FN2O8P | 326.17 | Covalent ternary-complex inhibitor with Cys195 + CH₂THF (Carreras & Santi 1995; Longley 2003) |
| 1. Active-site | BrdUMP | **93036** | `LHLHVDBXXZVYJT-RRKCRQDMSA-N` | C9H12BrN2O8P | 387.08 | Halogenated dUMP mimic (Santi & McHenry 1972) |
| 1. Active-site | Floxuridine (FdUR / 5-FdUrd) | **5790** | `ODKNJVUHOIMIIZ-RRKCRQDMSA-N` | C9H11FN2O5 | 246.19 | Nucleoside, no phosphate — expected to bind weakly, sanity-pair against 5-FdUMP |
| 1. Active-site (precursor sanity) | 5-FU | **3385** | `GHASVSINZRGABV-UHFFFAOYSA-N` | C4H3FN2O2 | 130.08 | Prodrug — expected to dock weakly, annotated as such |
| 2. Cofactor-site | Methotrexate (MTX) | **126941** | `FBOZXECLQNJBKD-ZDUSSCGKSA-N` | C20H22N8O5 | 454.4 | Classical antifolate; cross-target control (primarily DHFR) |
| 2. Cofactor-site | Raltitrexed (D16, Tomudex) | **104758** | `IVTVGDXNLFLDRM-HNNXBMFYSA-N` | C21H22N4O6S | 458.5 | Cofactor-site reference; **already bound in the holo receptor** |
| 2. Cofactor-site | Pemetrexed (S, Alimta) | **135410875** | *(PubChem returns no InChIKey for this record; v2 derives it via RDKit `MolToInchiKey` from the verified canonical SMILES — see §A2)* | C20H21N5O6 | 427.4 | Multi-target antifolate (TYMS + DHFR + GARFT). R2 footnote: CID 60843 is the (2R) enantiomer (clinical drug is S); CID 446556 is another S tautomer. |
| 2. Cofactor-site | Plevitrexed (ZD9331 / BGC9331) | **135430970** | *(derived via RDKit at A2)* | C26H25FN8O4 | 532.5 | TYMS-selective non-polyglutamatable antifolate (Jackman 1997) |
| 2. Cofactor-site | Nolatrexed (AG-337, Thymitaq) | **135400184** | *(derived via RDKit at A2)* | C14H12N4OS | 284.34 | Lipophilic non-classical antifolate (Webber 1996) |
| 2. Cofactor-site *(negative control)* | Ibuprofen | **3672** | `HEFNNWSXXWATRW-UHFFFAOYSA-N` | C13H18O2 | 206.3 | MW/logP-matched unrelated drug; should not bind cofactor site |
| 3. Dimer-interface | LR-derived peptide *(exact sequence resolved at A3)* | n/a (peptide) | computed at A3 | per A3 | per A3 | Cardinale 2011 lead, residue range mapped to human-TYMS 1HVY numbering via Biopython pairwise alignment in A3 |
| 3. Dimer-interface | Two short peptidomimetic analogs (≤ 5 residues) | n/a | computed at A3 | per A3 | per A3 | From Cardinale 2011 Table 1 + Salo-Ahen 2015 review |
| 4. Allosteric | *no Tier-1 anchors by design* | — | — | — | — | The absence is the teaching point (Brunn 2014) |

**v0/v1 wrong CIDs rejected (audit trail).** Eight of the ten v0/v1 CIDs were verified to point to *the wrong molecule*. The full rejected-CID list (with the actual compound each wrong CID returned) is in `anchor_compounds_verified.json` under `wrong_cids_rejected`. The egregious cases were:

- v0/v1 dUMP CID `22848` → returns a Solanum-alkaloid steroid C27H43NO8 (off by formula, off by MW, off by InChIKey).
- v0/v1 5-FdUMP CID `15718` → returns 2-(4-tert-butylphenoxy)acetic acid C12H16O3.
- v0/v1 BrdUMP CID `135398598` → returns GTP C10H12N5O13P3.
- v0/v1 nolatrexed CID `60198` → returns an estrogen-analog steroid C20H24O2.
- v0/v1 plevitrexed CID `122478` → returns an imidazole dimer C46H36Cl2N4O4.

Had the pipeline run on v1 CIDs the active-site "results" would have been ranking the *wrong molecules*. This is exactly the failure mode R2 demanded be closed before any compute.

**Tradeoff note (R1 [OPEN-B2] / sign-off #1, #3).** Moving nolatrexed and ZD9331 out of Strategy 1 shrinks the Strategy-1 anchor count from "six" (v0) to **four mechanistically-correct anchors + 5-FU as precursor sanity** (v1). Strategy 2 grows from five to six anchors (plus one negative control). This is the correct tradeoff: mis-bucketed anchors would have produced inflated active-site rankings and a misleading enrichment AUC. The total Tier-1 anchor count across the four strategies is **4 + 6 + 3 + 0 = 13** (down from 6 + 5 + 3 + 0 = 14 in v0); the matched-decoy expansion budget is unchanged.

### 1.2 Tier-2 (analogs + decoys)

| Source | Method | Target size | Endpoint |
| --- | --- | --- | --- |
| PubChem similarity neighbours of every Tier-1 anchor | Tanimoto ≥ 0.7 on Morgan-2 (radius 2, 2048-bit) via PubChem similarity REST | ~20 per strategy | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/<cid>/cids/JSON?Threshold=70` |
| DUD-E web service (Mysinger 2012, J. Med. Chem. 55:6582) — *attempted, fall back to RDKit on failure* | upload Tier-1 SMILES → DUD-E generates property-matched decoys | 50 decoys/anchor → cap at 200/strategy | `https://dude.docking.org/generate` |
| **RDKit DUD-E-style decoy generator** (the actual primary, per R3 pre-flight) | property-matched (MW±20, logP±0.5, HBA±2, HBD±2, RotB±2, formal-charge equal) random sampling from a 100k-compound PubChem subset | 30 decoys/anchor | local |
| ZINC15 fragment subset (Strategy 4 only) | MW ≤ 250 Da, logP ≤ 3.5, 200 random fragments | 200 | `https://zinc15.docking.org/substances/subsets/fragment/?count=200&format=smi` |

**R3 pre-flight result (2026-05-18): DUD-E web returns HTTP 500.** Per Stop Condition S3, v2 demotes DUD-E from "primary" (v1) to "attempted-with-fallback". The actual execution path is the RDKit DUD-E-style generator (third row), with DUD-E retried at the start of each strategy and used only if available. ROC-AUC and BEDROC numbers will therefore be RDKit-decoy-based; this is honest but slightly less benchmarkable to literature. Documented limitation.

### 1.3 Pre-docking filters (R1 sign-off #5) — applied between B and C

Every compound (Tier 1 and Tier 2) is passed through three structural-alert filters before ligand prep:

| Filter | Tool | Action |
| --- | --- | --- |
| PAINS A/B/C | `rdkit.Chem.FilterCatalog.FilterCatalogParams.PAINS_A/B/C` (Baell & Holloway 2010) | `pains_flag` column populated; compounds are *not* dropped — they are docked and ranked separately so the user sees how PAINS-flagged scaffolds behave |
| Brenk | `rdkit.Chem.FilterCatalog.FilterCatalogParams.BRENK` (Brenk 2008) | `brenk_flag` column populated, same handling |
| NIH MLSMR | `rdkit.Chem.FilterCatalog.FilterCatalogParams.NIH` | `nih_flag` |
| Lipinski / Veber (Strategy 4 fragments use Astex Rule-of-Three instead) | RDKit `Descriptors` | `lipinski_flag`, `veber_flag` (or `ro3_flag` for Strategy 4) |

### 1.4 Stereochemistry + tautomers + protomers (R1 sign-off #6, #7)

For every Tier-1 anchor with a defined stereocentre (notably raltitrexed and pemetrexed, both (S)), the corrector script **explicitly verifies the canonical SMILES enantiomer matches the marketed/literature stereochemistry** before docking. If PubChem returns a racemate, both enantiomers are docked and reported (`enantiomer` column = `S`, `R`, or `racemate`).

Tautomer + protomer enumeration at **pH 7.4 ± 0.5** uses:

```
rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumerator   # tautomers
Dimorphite-DL (Ropp 2019)                                       # protomers at pH 7.4 ± 0.5
```

All enumerated states are docked; the best-scoring state per compound is the row reported, and a `tautomer_id`, `protomer_id`, `n_states_docked`, `protonation_pH` set of columns captures the choice.

**Why this matters (R1 [OPEN-C1]).** Pemetrexed, raltitrexed, and nolatrexed all have multiple ionisable centres; the dominant species at physiological pH is not the default RDKit tautomer. This is the same failure mode that complicates the AD4 phosphate sign question documented in TECHNICAL_NOTES Phase-6c.

---

## 2. Per-strategy pipeline (the eight steps each strategy runs)

> v0 had a seven-step pipeline (A–G). v1 inserts step **A0 (positive-control re-dock gate)** in front, per R1 sign-off #4. Renumbering keeps the original A–G labels for the operating pipeline but A0 is a hard gate that must pass before A is even configured.

```
A0) POSITIVE-CONTROL RE-DOCK GATE — re-dock the native reference into its own site;
    require top1 heavy-atom pose RMSD ≤ 2.0 Å vs crystal pose; abort strategy on fail.
A1) Receptor protonation re-check — assert His196 + every His has explicit HID/HIE/HIP;
    re-run propka + pdb2pqr30 once if any His lacks tautomer assignment.
A2) CID → InChIKey verification — every Tier-1 PubChem CID is fetched, the canonical
    SMILES → RDKit → InChIKey is computed, and the InChIKey is cross-checked against
    a hard-coded literature reference table (pemetrexed free acid = QOFFJEBXNKRSPX-ZDUSSCGKSA-N, etc.).
    Gate-fail on mismatch.
A3) Dimer-interface contact-map computation — for Strategy 3 only, compute the 4 Å
    chain-A ↔ chain-B contact map via MDAnalysis and use the resulting residue list to
    define the interface box and the LR-derived peptide source range.
A)  Pocket definition (uses A3 output for Strategy 3)
B)  Compound set assembly (Tier 1 + Tier 2, with stereochem + tautomer + protomer
    enumeration per §1.4; then PAINS/Brenk/NIH/Lipinski/Veber filters per §1.3 — flagged,
    not dropped)
C)  Ligand prep (charge-delta gate per R1 [OPEN-C1])
D)  Docking (Vina, exh=32, multi-seed) — and HPEPDOCK for Strategy 3 ≥ 6-mer peptides
E)  Pose analysis — re-dock RMSD vs crystal, PROLIF + PLIP interaction profile,
    pose-cluster diversity (DBSCAN eps=2.0 Å, min_samples=2)
F)  Strategy-specific structural metrics (SASA / BSA / pocket-occupancy / FPocket druggability)
G)  Aggregate scoring + per-site Δ reference + ROC-AUC + BEDROC + apo_minus_holo gap
```

### A0 — Positive-control re-dock gate (R1 sign-off #4) — **hard gate**

The single most important addition from R1. Standard inhibitor-design SOP (Hawkins 2007; Warren 2006) demands that the native reference re-docks to within ≤ 2.0 Å RMSD of the crystal pose *before* any candidate is scored.

| Strategy | Native reference | Crystal source | Pass criterion |
| --- | --- | --- | --- |
| 1. Active-site | dUMP | 1HVY chain A UMP ligand (heavy atoms; UMP is the dehydroxylated form of dUMP in the deposited model — Kabsch-aligned heavy-atom RMSD on the matched atom set) | top1 pose RMSD ≤ 2.0 Å |
| 2. Cofactor-site | Raltitrexed (D16) | 1HVY chain A D16 ligand | top1 pose RMSD ≤ 2.0 Å |
| 3. Dimer-interface | LR-derived peptide ≥ 6 residues docked via HPEPDOCK | no crystal PPI-disruptor co-structure exists; gate replaced by **HPEPDOCK score sanity** (peptide must rank above a length-matched scrambled-sequence control by ≥ 1 standard deviation across 10 scrambles) | HPEPDOCK Δ vs scrambled-seq control ≥ 1 σ |
| 4. Allosteric | n/a (fragment screen, no native ligand) | gate replaced by **FPocket druggability sanity**: top-5 FPocket cavities on the resolved-chain-A surface must have druggability score ≥ 0.5 (Schmidtke & Barril 2010) | ≥ 1 FPocket cavity with druggability ≥ 0.5 outside the active-site/cofactor shells |

**On fail:** the strategy is aborted; the corrector or human investigator inspects the receptor (re-build / re-protonate / re-partition charges) before retry. Failure is *not* silently absorbed — the failure mode and the diagnostic plot are committed to `<strategy>/A0_redock_gate/`.

### A1 — Receptor protonation re-check (R1 §2 MEDIUM)

Phase-6c hardened **charges**, not necessarily His tautomers. Assertion script:

```python
# A1_check_his_tautomers.py — gate-fail on any HIS without HID/HIE/HIP
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)
pdb = parser.get_structure("apo", "06f_receptor_fixed/protein_dimer_apo_fixed.pdb")
his_residues = [r for r in pdb.get_residues() if r.resname == "HIS"]
assert len(his_residues) == 0, (
    f"Found {len(his_residues)} ambiguous HIS residues; expected only HID/HIE/HIP. "
    "Re-run propka + pdb2pqr30 at pH 7.4."
)
```

If the assertion fires, re-run `pdb2pqr30 --ph 7.4 --titration-state-method propka --with-ph 7.4 protein_dimer_apo.pdb protein_dimer_apo_fixed.pdb` once and re-export the PDBQT.

### A2 — PubChem CID → InChIKey verification gate (R2-closed)

Hard-coded CIDs are the most common silent error in computational chemistry pipelines (CIDs redirect to salts, racemates, the wrong tautomer). v2 ground-truth table (verified 2026-05-18; primary source = PubChem PUG REST cross-checked against IUPAC name + molecular formula + literature):

| Compound | Verified CID | Verified InChIKey (PubChem-returned, RDKit-derived where PubChem omits) |
| --- | --- | --- |
| dUMP | **65063** | `JSRLJPSBLDHEIO-SHYZEUOFSA-N` |
| 5-FdUMP | **8642** | `HFEKDTCAMMOLQP-RRKCRQDMSA-N` |
| BrdUMP | **93036** | `LHLHVDBXXZVYJT-RRKCRQDMSA-N` |
| Floxuridine | **5790** | `ODKNJVUHOIMIIZ-RRKCRQDMSA-N` |
| 5-FU | **3385** | `GHASVSINZRGABV-UHFFFAOYSA-N` |
| Methotrexate ((S)) | **126941** | `FBOZXECLQNJBKD-ZDUSSCGKSA-N` |
| Raltitrexed (free acid) | **104758** | `IVTVGDXNLFLDRM-HNNXBMFYSA-N` |
| Pemetrexed ((S)) | **135410875** | derived at runtime via `Chem.MolToInchiKey(Chem.MolFromSmiles(canonical_smiles))`; canonical SMILES baked in `anchor_compounds_verified.json` (CID 60843 is the (R) enantiomer — rejected) |
| Plevitrexed (ZD9331) | **135430970** | derived at runtime via RDKit (PubChem returns no InChIKey for this record) |
| Nolatrexed (AG-337) | **135400184** | derived at runtime via RDKit (PubChem returns no InChIKey for this record) |
| Ibuprofen | **3672** | `HEFNNWSXXWATRW-UHFFFAOYSA-N` |

Runtime gate behaviour (`scripts/v14/A2_verify_cids.py`):

```python
# R3-fix: use IsomericSMILES (stereo-preserving), not the obsolete ConnectivitySMILES.
# R3-fix: every anchor MUST carry a verified InChIKey before runtime; null is a gate-fail at load time.
import json, urllib.request
from rdkit import Chem

ground_truth = json.load(open("14_inhibitor_design/00_roadmap/anchor_compounds_verified.json"))

# Load-time gate (fail fast — before any PubChem call) — R3 HIGH-2 fix
for anchor in ground_truth["anchors"]:
    if anchor.get("cid") is None: continue          # peptides are exempt
    if anchor.get("inchikey") is None:
        raise RuntimeError(
            f"Ground-truth InChIKey missing for {anchor['display_name']} (CID {anchor['cid']}). "
            f"Run scripts/v14/A2_populate_inchikeys.py first to derive it from canonical SMILES."
        )

# Runtime gate — PubChem cross-check
for anchor in ground_truth["anchors"]:
    cid = anchor.get("cid")
    if cid is None: continue
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/InChIKey,IsomericSMILES,MolecularFormula/JSON"
    pubchem = json.loads(urllib.request.urlopen(url, timeout=15).read())["PropertyTable"]["Properties"][0]
    expected_key = anchor["inchikey"]   # mandatory after load-time gate
    # PubChem may omit InChIKey on stereo-specific records — fall back to RDKit-from-IsomericSMILES
    actual_key = pubchem.get("InChIKey") or Chem.MolToInchiKey(Chem.MolFromSmiles(pubchem["IsomericSMILES"]))
    if actual_key != expected_key:
        raise RuntimeError(f"InChIKey MISMATCH for {anchor['display_name']} (CID {cid}): expected {expected_key}, got {actual_key}")
    if pubchem["MolecularFormula"] != anchor["formula"]:
        raise RuntimeError(f"Formula MISMATCH for {anchor['display_name']}: expected {anchor['formula']}, got {pubchem['MolecularFormula']}")
```

Gate-fail action: abort the strategy, dump the mismatch table to `<strategy>/A2_cid_verification/MISMATCH.md`, and notify the corrector. The InChIKey is the ground truth; the CID may have been redirected upstream.

**Bootstrap step `scripts/v14/A2_populate_inchikeys.py` (R3 HIGH-2 fix).** Before A2 runs, this one-time script populates the null-InChIKey rows in `anchor_compounds_verified.json` by fetching their `IsomericSMILES` from PubChem and converting via `Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))`. The result is committed to the JSON so the load-time gate above can succeed. This bootstrap runs as part of A0 setup, before any docking is attempted.

### A3 — Dimer-interface contact map (R1 sign-off #11, [OPEN-A1])

For Strategy 3 *only*. v0 hard-coded the interface residue list from mixed-species papers; the reviewer correctly identified that Cardinale 2011 uses *E. coli* / *L. casei* numbering. v1 computes the interface programmatically against human-TYMS 1HVY numbering:

```python
# A3_dimer_interface.py
import MDAnalysis as mda
u = mda.Universe("06f_receptor_fixed/protein_dimer_apo_fixed.pdb")
chainA = u.select_atoms("protein and segid A")
chainB = u.select_atoms("protein and segid B")
# 4 Å contact map heavy-atom – heavy-atom
contacts = mda.analysis.contacts.Contacts(
    u, select=("protein and segid A and not name H*", "protein and segid B and not name H*"),
    refgroup=(chainA, chainB), radius=4.0
)
# residues from chain A whose any heavy atom is within 4 Å of any heavy atom of chain B
interface_A = sorted(set(a.resid for a in chainA if any(
    (a.position - b.position).dot(a.position - b.position) ** 0.5 <= 4.0
    for b in chainB
)))
# box centre = midpoint of chain-A interface residue Cα centroid and chain-B equivalent
```

The resulting `interface_residues_chainA.csv` and `interface_residues_chainB.csv` go into `03_dimer_interface/A3_contact_map/` and the box centre + dimensions in §2A are populated from this output (not the literature list).

**Expected output (post-R1 prediction; final list populated by A3 at runtime).** Approximate residue ranges per the reviewer's expert estimate: chain-A residues {21–27, 76–82, 174–180, 200–207, 250–256} and symmetric chain-B partners. This expectation is documented so the A3 output can be visually checked, *not* hard-coded into the pipeline.

The LR-derived peptide source range (Strategy-3 Tier-1) is then defined as the chain-A subrange that aligns to the Cardinale 2011 *E. coli* peptide via a pairwise human-vs-*E. coli* TYMS sequence alignment (Biopython `pairwise2.align.globalxx`), with the human-numbering residue range cited explicitly in the strategy README.

### A — Pocket definition (v1, after R1)

| Strategy | Box centre (Å) | Box size (Å) | Reference atoms used |
| --- | --- | --- | --- |
| 1. Active-site | (−0.137, +4.232, +15.159) | **22 × 22 × 22** (single canonical, see footnote) | Cα centroid of residues [80, 87, 109, 135, 175, 176, 195, 196, 214, 215, 217, 218, 221, 225, 226, 258] (chain A — Phase 5–7 canonical) |
| 2. Cofactor-site | centroid of raltitrexed D16 (chain A, **computed once from the holo receptor** `06f_receptor_fixed/cofactor_A.pdbqt` and reused unchanged for both apo and holo dockings — the apo receptor has no D16 by definition) | 22 × 22 × 22 | All non-H atoms of D16:A |
| 3. Dimer interface | midpoint between Cα centroids of chain-A and chain-B interface residues *as computed by A3* | 26 × 22 × 22 (longer along A↔B axis) | A3 output |
| 4. Allosteric | Top FPocket cavity centre on chain A surface, **excluding** the active-site and cofactor-site 8 Å shells. Note: residues 1–26 disordered in 1HVY → no FPocket cavity expected on the mRNA-binding loop. | 20 × 20 × 20 | FPocket cavity centroid, ranked by druggability score |

**Footnote — box size choice for Strategy 1 (R1 §0 MEDIUM, sign-off #12).** v0 quoted both 22³ (v5 canonical) and 18³ (Phase 7) as "still in play". v1 commits to **22³**. Justification: 22³ is the canonical size used in Phases 5, 6, 6c, 6d, 6e and 7; cross-strategy comparison demands a single box. The Phase-7 structural-bioinformatician verifier (`reviews/04_structural_bioinformatician.md §"Docking parameters (Vina)"`) flagged that 22³ *admits* off-site minima 7–8 Å from mode 1; we mitigate this by reporting the per-pose distance from the active-site Cα centroid in column `pose_distance_from_box_centre`, and any pose > 7 Å from the centroid is flagged `off_site_minimum = True`. The reviewer can re-audit those rows specifically rather than the whole panel.

### A.x — Optional second receptor: AlphaFold (R1 §"Additional findings" LOW)

For every strategy's top hit (after step G), also dock against the AlphaFold receptor used in `12_phase7/03_alphafold/`. If the score diverges from the 1HVY-receptor score by > 2 kcal/mol, the divergence is logged as a `receptor_model_sensitivity` column and surfaced in the strategy README — that is a real, publishable finding about model-vs-crystal sensitivity, not a bug.

### B — Compound set assembly

See §1.1 (Tier-1, R1-corrected) + §1.2 (Tier-2, DUD-E for Strategy 1 + RDKit fallback) + §1.3 (PAINS/Brenk/NIH/Lipinski filters, flagged not dropped) + §1.4 (stereochem + tautomer + protomer enumeration at pH 7.4 ± 0.5).

### C — Ligand prep (with R1 [OPEN-C1] charge-delta gate)

```
for each (compound × tautomer × protomer × enantiomer):
  RDKit  — sanitize, embed with ETKDGv3, optimise with MMFF94s (max 1000 iter)
  obabel — protonate at pH 7.4, assign Gasteiger charges (if Dimorphite-DL already did
           the protonation, obabel only charges — do not re-protonate)
  meeko  — mk_prepare_ligand.py → PDBQT with explicit rotatable bonds
  Assertions:
    (i)  |total charge − formal charge| < 0.05 e  (gate-fail else)
    (ii) record formal_charge_input, total_q_pdbqt, delta_q in CSV;
         surface any row with |delta_q| > 0.10 e in a separate review-me file.
```

### D — Docking

**Same canonical Vina invocation as Phase 5–7** (so results are scale-comparable):

```bash
vina \
  --receptor 06f_receptor_fixed/protein_dimer_{apo,holo}_fixed.pdbqt \
  --ligand   14_inhibitor_design/<strategy>/ligands/<compound>_<state>.pdbqt \
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
- **State:** strategies 1 & 2 dock against both `apo` and `holo_raltitrexed_bound` receptors (so dUMP-displacement vs cofactor-displacement vs raltitrexed-displacement can be inferred from the per-state Δ). Strategies 3 & 4 dock against `apo` only (cofactor irrelevant for dimer-interface and allosteric surface). Strategy 4 fragments: **1 seed only** (fragment screen scale), the multi-seed convergence test is restricted to the top-20 fragments by absolute Vina score in a second pass.
- **Aggregation:** report `top1`, `top3_mean`, `top5_mean`, `mean_over_seeds(top1)`, `sd_over_seeds(top1)` per (compound × tautomer × protomer × enantiomer × state).

**Flex-residue follow-up (R1 [OPEN-D1] resolved).** Rigid-receptor Vina is the primary protocol. Flex-residue follow-up is run *only* on compounds with rigid `Δ vs strategy reference ≤ −0.85 kcal/mol` (per the reviewer's tightened threshold), and is **capped at top 5 per strategy** — i.e. at most 5 × 4 = 20 flex runs across the phase. Flex side chains: catalytic Cys195 + His196 + Arg175/176/215 for Strategy 1; D16-contact residues Phe80, Trp80-loop, Asp225 for Strategy 2; no flex for Strategies 3 and 4 (peptide / fragment scale makes flex meaningless).

**Strategy 3 docking deviation (R1 [OPEN-B2] resolved; R2 HIGH-finding closed; R3 pre-flight update).** Vina is unreliable above 4-residue peptides (Hassan 2017, J. Comput. Chem. 38:1278, median pose accuracy drops below 2 Å for ≥ 5-mers). Strategy 3 therefore uses, in priority order:

- **HPEPDOCK** (Zhou *et al.* 2018, Nucleic Acids Res. 46:W443) as the *intended* primary docker for peptides ≥ 6 residues. Web service at `http://huanglab.phys.hust.edu.cn/hpepdock/`. **R3 pre-flight 2026-05-18: HPEPDOCK web server unreachable (curl exit code 7, connection failed).** Per the envelope below, this triggers immediate fallback to CABS-dock.
- **CABS-dock** (Kurcinski *et al.* 2015, Nucleic Acids Res. 43:W419) at `http://biocomp.chem.uw.edu.pl/CABSdock` — R3 pre-flight 2026-05-18: returns HTTP 302 (server alive). **Promoted to actual primary peptide docker for v2 execution.**
- **Vina** retained for ≤ 5-residue mimetics only.
- **Fragment-decomposition fallback** (R3 contingency): if both web services are down at execution time, the LR octapeptide is decomposed into overlapping 4-mers (positions 1–4, 2–5, 3–6, 4–7, 5–8) and each is docked into the interface box with Vina. Cross-strategy ranking remains broken — the result is qualitative ("does the canonical fragment landscape rank below the scrambled-control landscape?").

**HPEPDOCK web operations envelope (R2-closed).** Concrete failure handling:

| Setting | Value |
| --- | --- |
| Per-submission HTTP timeout (POST + initial response) | 60 s |
| Per-job result-polling interval | 5 min |
| Per-job total wall-time before TIMEOUT | **4 h** |
| Retries on 5xx / connection error | 3, with 30-min exponential backoff (30, 60, 120 min) |
| Behaviour on TIMEOUT | mark row `hpepdock_status = TIMEOUT`; do NOT silently drop |
| Behaviour if HPEPDOCK unreachable for > 15 min cumulative | switch to **CABS-dock** for remaining jobs in this strategy; record `engine = CABS-dock` in row |
| Behaviour if CABS-dock also unreachable | abort Strategy 3 peptide track entirely; document under **extended Stop Condition S3** (web service unavailable on date X) |
| Triggered abort threshold | if > 50 % of submitted peptide jobs hit TIMEOUT, abort the peptide track and document |

**HPEPDOCK A0 re-dock control (R2-closed).** The Strategy-3 equivalent of the small-molecule A0 gate: re-dock the LR octapeptide into 1HVY. Pose-RMSD criterion loosened to ≤ 4 Å (peptides are flexible; published HPEPDOCK benchmarks accept ≤ 4 Å). Cardinale 2011 supplementary figure provides the reference pose; if no Cardinale-validated pose is available for human-TYMS, the gate degrades to the scrambled-control check below.

**HPEPDOCK scrambled-sequence control (R2-closed concrete spec).** For each peptide submitted:
1. Generate a shuffled-residue control using `numpy.random.default_rng(seed=42).permutation(list(peptide_sequence))`.
2. Submit the shuffled control alongside as a separate HPEPDOCK job.
3. Compare top-1 HPEPDOCK score: expected behaviour is `canonical_top1 ≤ shuffled_top1 − 2 kcal/mol` (HPEPDOCK kcal/mol-equivalent units).
4. If not, flag `peptide_specificity_unreliable = True` in row; rank with caution.

**Why this tradeoff.** HPEPDOCK and Vina report different energy scales; cross-strategy ranking on raw kcal/mol is therefore broken for Strategy 3 peptides. We accept this: Strategy 3 is internally ranked (canonical vs scrambled control) but is not numerically compared to Strategies 1/2/4. The headline plot Section 2G figure 2 separates Strategy 3 onto its own subpanel.

### E — Pose analysis

For every compound that beats its strategy reference by Δ ≥ 0.85 kcal/mol:

1. **Re-dock RMSD vs crystal pose** (where a crystal exists): chain-A copy of the bound ligand in 1HVY (dUMP, raltitrexed), 1JU6 (5-FdUMP), 1HW3 (raltitrexed, if available), etc. RMSD computed via RDKit's symmetry-corrected `GetBestRMS` after Kabsch heavy-atom alignment. **Note:** this is *in addition to* the A0 gate, which only checks the native reference for the strategy.
2. **Interaction profile — both PROLIF and PLIP** (R1 §2 LOW; revised from PROLIF-only).
   - `prolif` (RDKit-based) computes per-residue interaction fingerprint (HBond donor/acceptor, hydrophobic, π-stacking, salt-bridge, halogen). Output: per-compound `prolif_fingerprint.json` + a heatmap of interaction-type × residue, with the strategy reference fingerprint as the comparison column.
   - **PLIP** (Salentin *et al.* 2015) is run on the same complexes for a human-readable per-complex report (`<compound>_<seed>_plip.txt` + `_plip.xml`). PLIP gives the field-standard inhibitor-paper interaction summary; PROLIF gives the cross-compound heatmap matrix. The two outputs are saved side-by-side and any disagreement (e.g. PROLIF reports an HBond that PLIP misses) is flagged in column `prolif_plip_agreement`.
   - For Strategy-1 active-site compounds, the **separate step E1b** below quantifies whether the top pose would have made a water-bridged H-bond in the crystal.
3. **Pose-cluster diversity** — DBSCAN on heavy-atom RMSD across the 5 seeds × top-3 modes. **Parameters: `eps = 2.0 Å`, `min_samples = 2`** (R1 §2 LOW). Report number of distinct clusters (>2 Å); single tight cluster across seeds = high-confidence pose, many loose clusters = Vina searching, not converging.

#### E1b — Crystal water-bridge check (Strategy 1 only) — R2 HIGH-finding fix

PROLIF analyses the *waterless* docked complex; it cannot, by construction, flag *missing* crystal-water bridges. Hence a separate MDAnalysis-based script (`scripts/v14/E1b_water_bridge_check.py`) runs only for Strategy-1 active-site compounds:

```python
# E1b_water_bridge_check.py
#
# IMPORTANT (R3-fix): Vina pose PDBQTs are LIGAND-ONLY. The pose is already in
# receptor-frame coordinates because Vina docked into the receptor's coord system.
# We do NOT load the pose-PDBQT as an mda.Universe — we load 1HVY directly,
# read the docked pose heavy-atom coordinates with RDKit (or by hand-parsing
# the PDBQT ATOM lines), then compute distances against 1HVY waters in their
# native frame. No alignment is needed because the receptor coord frame is
# preserved end-to-end (Phase 6c -> Phase 14).
import json, glob, numpy as np
import MDAnalysis as mda
from rdkit import Chem

xtal = mda.Universe("input/1HVY_with_waters.pdb")
xtal_waters_o = xtal.select_atoms("resname HOH and name O")
xtal_tyr258_oh = xtal.select_atoms("resid 258 and name OH and segid A")[0]

def parse_pose_pdbqt_heavy(path):
    """Return (N,3) ndarray of heavy-atom coords from a Vina pose PDBQT, MODEL 1 only."""
    coords = []
    in_model1 = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("MODEL 1"): in_model1 = True; continue
            if line.startswith("ENDMDL") and in_model1: break
            if in_model1 and (line.startswith("ATOM") or line.startswith("HETATM")):
                # AD-type is the last column; non-polar H is 'H' or 'HD', we skip H
                element = line[76:78].strip() if len(line) >= 78 else ""
                if element.startswith("H"): continue
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append((x, y, z))
    return np.array(coords)

for pose_pdbqt in glob.glob("01_active_site/docked/*_apo_seed42.pdbqt"):
    pose_heavy = parse_pose_pdbqt_heavy(pose_pdbqt)
    flagged_waters = []
    for water_o in xtal_waters_o:
        d_tyr = np.linalg.norm(water_o.position - xtal_tyr258_oh.position)
        if d_tyr > 3.5: continue
        d_pose = float(np.min(np.linalg.norm(pose_heavy - water_o.position, axis=1)))
        if d_pose <= 3.5:
            flagged_waters.append({"resid": int(water_o.resid), "d_tyr": float(d_tyr), "d_pose": d_pose})

    out = {
        "pose": pose_pdbqt,
        "water_bridge_lost": len(flagged_waters) > 0,
        "water_bridge_residues": flagged_waters,
    }
    json.dump(out, open(pose_pdbqt.replace(".pdbqt", "_waterbridge.json"), "w"), indent=2)
```

**Why no alignment is needed.** The Phase-6c receptor (`06f_receptor_fixed/protein_dimer_apo_fixed.pdbqt`) was built directly from 1HVY without any rigid-body transform. The Vina docking box (centre `(-0.137, 4.232, 15.159)`) is also in 1HVY frame. Vina poses are therefore returned in 1HVY frame and can be compared to 1HVY waters without alignment. We assert this once at pipeline start by checking that the receptor PDB header `ORIGX1/2/3` rows are identity matrices.

Per the v1 decision waters are still removed for docking; E1b *annotates* whether the removed waters would have mattered. Output: per-pose `_waterbridge.json` + a `water_bridge_lost` column in `results.csv`. The known Tyr258 ↔ dUMP O4 bridge is the most common one expected to fire.

### F — Strategy-specific structural metrics

| Strategy | Metric | Tool | Why |
| --- | --- | --- | --- |
| 1. Active-site | Pocket-occupancy ratio (ligand SASA buried / ligand SASA free) | `freesasa` (Python) | A genuine substrate-mimetic should bury ≥ 80 % of its surface, as dUMP does. |
| 2. Cofactor-site | ΔSASA on chain A residues 80, 130–134, 218–225 (mTHF binding face) | `freesasa` Δ between apo and complex | Antifolate occupancy should specifically de-expose these residues. |
| 3. Dimer interface | **PPI metrics**: BSA (buried surface area) on dimer interface residues from A3, plus per-residue contact count between chain A and chain B with the ligand present vs absent | `freesasa` + MDAnalysis (PRODIGY-LIG attempted as best-effort, see below) | A real PPI disruptor should *reduce* dimer interface BSA when bound. |
| 4. Allosteric | FPocket druggability score of the cavity the top pose lands in; distance from active-site Cα centroid (must be ≥ 8 Å) | `fpocket` (Homebrew) | An allosteric hit should land in a druggable pocket distinct from active/cofactor. |

**Tooling decisions (R1 [OPEN-F1], [OPEN-F2] resolved).**

- **FPocket**: install via `brew install fpocket` (arm64-darwin Homebrew bottle exists, R1 verdict). Pipeline gate: A0 for Strategy 4 fails if `which fpocket` returns nothing or if FPocket runs but produces zero cavities with druggability ≥ 0.5. Fallback if install fails: MDpocket (web service); install failure is logged with OS version.
- **PROLIF**: install via `PIP_BREAK_SYSTEM_PACKAGES=1 pip install prolif`. Verified runtime on Apple Silicon Py 3.14 (Phase 14 setup, 2026-05-18).
- **PLIP install gate (R2 MEDIUM-finding closed).** PLIP depends on the **OpenBabel Python bindings**, which on arm64-darwin Python 3.14 are *not* installable via pip alone — they require a Homebrew openbabel + a manual symlink. The install gate is:

```bash
# Step 1 — try pip
PIP_BREAK_SYSTEM_PACKAGES=1 pip3 install --user plip
# Step 2 — verify the openbabel Python binding loads
python3 -c "from openbabel import openbabel" 2>/dev/null
if [ $? -ne 0 ]; then
  # Step 3 — try the alternate import path some openbabel builds expose
  python3 -c "import openbabel" 2>/dev/null
  if [ $? -ne 0 ]; then
    echo "OpenBabel Python bindings missing on this host; PLIP cannot run."
    echo "skip_plip" > 14_inhibitor_design/00_roadmap/PLIP_STATUS
  fi
fi
# Step 4 — if step 2 or 3 succeeded, verify PLIP loads
python3 -c "from plip.structure.preparation import PDBComplex" 2>/dev/null || echo "skip_plip" > 14_inhibitor_design/00_roadmap/PLIP_STATUS
```

If `PLIP_STATUS = skip_plip` is present, all PLIP columns in `results.csv` are blanked with `null_reason = "plip_unavailable_arm64_py314_openbabel_binding"` per Stop Condition S3, and the pipeline relies on PROLIF alone. PROLIF is RDKit-based and has no OpenBabel dependency, so its operation is independent of the PLIP gate.
- **PRODIGY-LIG**: web service, ≤ 30 submissions/day. Used as best-effort for Strategy 3 only. Failure → column populated with `null_reason = "PRODIGY-LIG quota exceeded"` per Stop Condition S3.
- **FoldX 5**: Schymkowitz lab ships x86_64 binary only; on arm64-darwin this works under Rosetta 2 but is slow. **Decision: not used in v1.** The interface ΔΔG column is left blank with `null_reason = "FoldX 5 not available on arm64-darwin"`; BSA + interface residue-residue contact change (MDAnalysis) covers the same scientific question with a known smaller resolution. Documented limitation — accept null result per Stop Condition S1.
- **GNINA / autogrid4**: not available on Apple Silicon; not used. We document this in `TECHNICAL_NOTES.md` Phase-14 section. Per the original constraint, Vina is the only docking engine.

### G — Aggregate scoring + per-site Δ + ROC-AUC + apo-vs-holo

Per-strategy `results.csv` schema (R1 sign-off #10):

```
compound_name, pubchem_cid, inchikey, canonical_smiles,
mw, logp, hba, hbd, rotb, tpsa,
tier, anchor_for,
stereochem (S/R/racemate), enantiomer, tautomer_id, protomer_id,
n_states_docked, protonation_pH,
formal_charge_input, total_q_pdbqt, delta_q,
pains_flag, brenk_flag, nih_flag, lipinski_flag, veber_flag, ro3_flag,
state, seed, engine (Vina/HPEPDOCK),
top1_kcalmol, top3_mean, top5_mean,
mean_over_seeds_top1, sd_over_seeds_top1,
pose_cluster_count, pose_rmsd_vs_xtal,
sasa_buried, sasa_buried_pct,
ppi_bsa, ppi_delta_bsa,
fpocket_druggability, distance_from_active_site_ca,
prolif_overlap_with_reference, plip_summary, prolif_plip_agreement, water_bridge_lost,
delta_vs_reference_top1, delta_vs_reference_significant_p085,
apo_minus_holo_top1, cryptic_pocket_flag,
pose_distance_from_box_centre, off_site_minimum,
receptor_model_sensitivity,
notes, null_reason
```

(The `delta_vs_reference_*` columns are populated against the per-strategy reference per §0 — dUMP for site 1, raltitrexed for site 2, LR peptide for site 3, none for site 4.)

**Per-strategy enrichment metrics (R1 §"Additional findings" MEDIUM, sign-off #10).** For each strategy with Tier-1 actives + Tier-2 decoys, compute:

- `roc_auc_top1` — ROC area-under-curve on label = "Tier-1 active" vs score = `top1_kcalmol` (more negative = predicted active)
- `bedroc_alpha20` — Boltzmann-Enhanced Discrimination of ROC at α=20 (Truchon & Bayly 2007), captures early enrichment

Both go in a per-strategy summary file `01_active_site/enrichment.csv` (and equivalents for Strategies 2 and 3; Strategy 4 has no actives so no enrichment metric).

**Cryptic-pocket / induced-fit indicator (R1 §"Additional findings" MEDIUM).** For Strategies 1 and 2 (apo + holo states), `apo_minus_holo_top1` quantifies the score gap; rows with gap > 2 kcal/mol get `cryptic_pocket_flag = True`. This is the cheapest free signal for an induced-fit / cryptic-pocket effect and surfaces it without any extra docking.

Final `14_inhibitor_design/05_aggregate/master.csv` is the union of the four per-strategy CSVs plus a `strategy` column.

Headline plots (saved under `figures/`):

1. Per-strategy distribution of `top1_kcalmol` (violin × 4 strategies), with the per-strategy reference line — `fig_distributions.png`
2. **Per-strategy `delta_vs_reference_top1` ranked horizontal bar**, colour-coded by Tier and significance gate, with Strategies 1/2 on the same panel (both Vina-scored) and Strategy 3 on a separate subpanel (HPEPDOCK-scored, different scale), Strategy 4 by absolute `top1_kcalmol` instead — `fig_delta_ranking.png`
3. PROLIF + PLIP interaction-fingerprint heatmaps (one per strategy's top hit + the per-strategy reference) — `fig_interaction_heatmaps.png`
4. SASA / BSA / pocket-occupancy / FPocket panel (strategy-specific) — `fig_structural_metrics.png`
5. Per-strategy pose overlay PyMOL render (top hit vs reference vs Tier-1 anchor where applicable) — `fig_pose_overlays.png`
6. **ROC + BEDROC curves per strategy with actives + decoys** — `fig_enrichment.png`

---

## 3. Reviewer / corrector loop

> Exactly per user instruction: one biologist+bioinformatician **reviewer agent**, one **corrector agent**, unbounded loop until reviewer signs `PASS` on every section.

Each cycle:

1. **Reviewer agent** receives the current roadmap and writes `reviews/00_roadmap_R<n>_review.md` with verdicts `PASS / CONDITIONAL_PASS / FAIL` per section and per `[OPEN-*]` tag.
2. **Corrector agent** receives the review and writes the next roadmap version (`ROADMAP.md` overwritten; previous saved as `ROADMAP_v<n-1>.md`) plus `reviews/00_roadmap_R<n>_corrector_changelog.md`.
3. Repeat until the reviewer's overall verdict is `PASS` (no FAIL, no CONDITIONAL_PASS, every `[OPEN-*]` tag resolved or explicitly converted to a "documented limitation").

The same loop runs **again** after each strategy is executed (review the *results*, not the plan): cycle continues until reviewer signs the analysis off.

---

## 4. Compute budget and stop conditions (R1-revised)

| Stage | Estimated wall-time (single-Mac, M-series, exh=32) | Comments |
| --- | --- | --- |
| Compound assembly + prep (all 4 strategies, including stereochem + tautomer + protomer enumeration + PAINS/Brenk/NIH/Lipinski filtering) | ~30–45 min | PubChem REST + RDKit/Dimorphite-DL/meeko |
| **A0 + A1 + A2 + A3 gates** (re-dock dUMP, re-dock raltitrexed, contact-map, His tautomer assertion, CID verification, FPocket sanity for S4) | ~1–2 h | dUMP + raltitrexed re-dock = ~2 × Vina runs at exh=32; A3 contact map < 1 min; A2 PubChem fetch ~5 min |
| Strategy 1 docking (4 anchors + 5-FU sanity + 200 DUD-E decoys, each × N enumerated states × 2 receptor states × 5 seeds) | ~10–15 h | ~250 compound-states × 2 states × 5 seeds = ~2500 Vina runs (offset by fewer Tier-1 anchors but tautomer enumeration expands; DUD-E web takes 3–24 h turnaround — submit early) |
| Strategy 2 docking (6 anchors + 1 neg control + 20 PubChem analogs, each × N states × 2 receptor states × 5 seeds) | ~6–10 h | smaller decoy set than S1 (RDKit fallback) |
| Strategy 3 docking — HPEPDOCK for ≥ 6-mer peptides (3 peptide anchors + 10 mimics × 1 state × HPEPDOCK default ensemble) + Vina for ≤ 5-mer mimetics (× 5 seeds) | **~12–20 h** (revised up from v0's 3–5 h per R1 §4 MEDIUM; HPEPDOCK web latency dominates, not local CPU) | HPEPDOCK web submission queue + processing; budget assumes overnight quota use |
| Strategy 4 docking (200 ZINC15 fragments × 1 state × 1 seed; top 20 re-docked × 5 seeds) | ~6–10 h | small ligands fast |
| Pose analysis + PROLIF + PLIP + figures + enrichment | ~2–4 h | adds PLIP (per-complex) and ROC/BEDROC computation vs v0 |
| Flex-residue follow-up (≤ 5 compounds × 4 strategies = 20 runs at ~30 min each) | ~10 h max | gated per D1 resolution |
| Optional AlphaFold-receptor sanity dock for the top hit per strategy (4 runs × 5 seeds = 20 runs) | ~2 h | per R1 §"Additional findings" |
| **Total** | **~50–80 h sequential** (revised from v0's 20–35 h) | parallelisable across strategies if budget permits 4-way fork; HPEPDOCK web latency is unavoidable |

**Tradeoff note on the budget revision.** v0 underestimated by ~2×. The extra time comes from (i) tautomer + protomer enumeration multiplying the compound count, (ii) HPEPDOCK web turnaround for Strategy 3, (iii) DUD-E web turnaround for Strategy 1, (iv) ROC/BEDROC + PLIP + AlphaFold-sanity additions. None of these are negotiable — the reviewer flagged each as a required fix.

### Stop conditions (R1-revised; honoured *up-front*)

- **S1.** Rigid-receptor Vina cannot resolve binding-mode differences below Vina's ±0.85 kcal/mol noise floor. If a strategy yields no compound with `Δ vs strategy reference ≤ −0.85`, the reviewer must accept "null result, see Phase 14 Limitations" rather than loop indefinitely on the same compute.
- **S2.** Phase 14 inherits Phase 8b's documented flex-residue compute budget (30 min per compound). Flex-only-on-rigid-hits gate (§D above) is mandatory; capped at 5 compounds per strategy.
- **S3.** Tools that aren't installed or aren't compatible with arm64-darwin (FoldX 5, GNINA, autogrid4, PRODIGY-LIG over quota) are *documented as missing* and the corresponding column in the results CSV is left blank with a `null_reason` value — not silently skipped.
- **S4 (new, R1 §4 LOW).** **If a Tier-1 anchor (known active) fails to dock into its own site with `Δ vs reference ≤ −0.85`, the docking protocol is the suspect, not the chemistry** — flag the strategy as "protocol failure under investigation" and do not report any Tier-2 ranking from that strategy until the anchor docks correctly. The most likely culprits are (i) wrong tautomer/protomer selected at C-step, (ii) box size, (iii) receptor protonation regression. Pipeline diagnostics auto-runs on S4 firing: re-runs A2 (CID/InChIKey), prints the tautomer scan, re-runs A1 (His tautomer), re-runs A0 (re-dock RMSD).

---

## 5. Deliverables (what gets committed back)

- `14_inhibitor_design/00_roadmap/ROADMAP.md` — this file, final version after reviewer PASS.
- `14_inhibitor_design/00_roadmap/ROADMAP_v<n>.md` — every previous draft (audit trail).
- `14_inhibitor_design/00_roadmap/reviews/` — all review iterations + corrector changelogs verbatim.
- `14_inhibitor_design/{01..04}_<strategy>/`:
  - `A0_redock_gate/` — re-dock RMSD diagnostic, pose overlay PNG
  - `A1_protonation_check/` — His tautomer assertion log
  - `A2_cid_verification/` — InChIKey verification table + any mismatches
  - `A3_contact_map/` (Strategy 3 only) — interface residue lists
  - `ligands/` — all input SDF + PDBQT (per tautomer × protomer × enantiomer)
  - `docked/` — all output PDBQT + logs
  - `analysis/` — PROLIF + PLIP outputs, SASA tables, PPI metrics, pose-RMSD tables, FPocket cavities
  - `enrichment.csv` — ROC-AUC + BEDROC per strategy (S1–S3)
  - `results.csv` — strategy-level master with the schema in §2G
  - `README.md` — strategy-level educational write-up (limitations, Tier-1 chemistry, Δ-reference choice)
- `14_inhibitor_design/05_aggregate/master.csv` — cross-strategy union
- `14_inhibitor_design/figures/` — six headline plots (fig_distributions, fig_delta_ranking, fig_interaction_heatmaps, fig_structural_metrics, fig_pose_overlays, fig_enrichment)
- **Root `README.md` update** — educational summary of Phase 14 + commands + figures.
- **Root `TECHNICAL_NOTES.md` update** — Phase 14 caveats: HPEPDOCK-vs-Vina scale mismatch, FoldX 5 / GNINA / autogrid4 absence on arm64-darwin, crystal-water removal decision, 22³ off-site-minimum flag, AlphaFold-vs-1HVY sensitivity findings.
- **Root `CHANGELOG.md` update** — commit-level summary.

---

## 6. Resolution of all v0 `[OPEN-*]` tags (no open tags remain in v1)

| Tag | v0 question | v1 resolution |
| --- | --- | --- |
| `[OPEN-A1]` | Confirm or veto literature-prior dimer-interface residue list | **Resolved** — replaced with §A3 programmatic 4 Å contact-map computation in human-TYMS 1HVY numbering. Literature-prior list retained only as a visual-check expectation. |
| `[OPEN-B1]` | Add a CID-to-InChIKey verification gate | **Resolved** — §A2 is now a hard gate with literature-InChIKey reference table; pemetrexed CID corrected to 60843 (free acid). |
| `[OPEN-B2]` | Peptides in Vina | **Resolved** — Strategy 3 ≥ 6-mers use HPEPDOCK (Zhou 2018); ≤ 5-mer mimetics use Vina; cross-strategy ranking separated per §D tradeoff note. |
| `[OPEN-C1]` | Charge-delta surface in CSV | **Resolved** — §C now records `formal_charge_input`, `total_q_pdbqt`, `delta_q` and gate-fails on |Δq| > 0.05 e; rows with |Δq| > 0.10 e surfaced to a separate review-me file. |
| `[OPEN-D1]` | Flex-residue follow-up gating threshold | **Resolved** — rigid first for full panel; flex follow-up only on rigid `Δ ≤ −0.85`, capped at top 5 per strategy (≤ 20 flex runs total). |
| `[OPEN-F1]` | PRODIGY-LIG / FoldX availability + fallback | **Resolved** — PRODIGY-LIG used best-effort within ≤ 30 submissions/day; FoldX 5 not available on arm64-darwin → documented limitation, ΔΔG column blank with `null_reason`; MDAnalysis BSA + contact-count change covers the science. |
| `[OPEN-F2]` | `fpocket` install-or-skip | **Resolved** — `brew install fpocket` is the install step (arm64-darwin bottle exists); Strategy-4 A0 gate fails if FPocket is absent or returns zero druggable cavities. |

**Items the reviewer raised but were declined or scaled back:**

- **GNINA** — declined. Apple Silicon limit, project-wide engine constraint. Documented as a known limitation in TECHNICAL_NOTES.
- **AutoGrid4 maps** — declined. Same reason.
- **FoldX 5 ΔΔG** — declined as a hard pipeline step. Documented limitation; MDAnalysis BSA covers the same scientific question with smaller resolution.

---

## 7. Round-1 sign-off checklist (R1 §"Sign-off requirements") — status

For audit transparency, the 14 sign-off requirements from the R1 review map to v1 sections as follows:

| # | R1 sign-off requirement | v1 location |
| --- | --- | --- | 
| 1 | Swap Strategy 1 anchor list (drop 5-FU, move nolatrexed/ZD9331 to S2, add dUMP+BrdUMP) | §1.1 table |
| 2 | Verify every PubChem CID by InChIKey | §A2 + table in §A2 |
| 3 | Fix the dimer-interface peptide, switch primary docker to HPEPDOCK | §1.1 row 11 + §A3 + §D (Strategy 3 deviation) |
| 4 | Add step A0 = positive-control re-dock gate with ≤ 2 Å RMSD | §A0 (hard gate) |
| 5 | Add PAINS/Brenk filter step between B and C | §1.3 |
| 6 | Add tautomer + protomer enumeration at pH 7.4 ± 0.5 | §1.4 |
| 7 | Add enantiomer verification for stereoactive anchors | §1.4 |
| 8 | Document crystal-water handling explicitly | §0 "Crystal-water handling" |
| 9 | Split `Δ_vs_dUMP` into per-site Δ references | §0 "Per-site Δ reference convention" + §G plot 2 |
| 10 | Add roc_auc_top1 + bedroc_alpha20 + apo_minus_holo_top1 + pains_flag + lipinski_flag + tautomer_id + enantiomer columns | §G schema + §G "Per-strategy enrichment metrics" |
| 11 | Compute dimer-interface residue list from 4 Å contact map | §A3 |
| 12 | Pick one box size for Strategy 1 | §A "Footnote — box size choice" → 22³ chosen |
| 13 | Re-verify apo receptor's His tautomer assignment at pH 7.4 | §A1 |
| 14 | Re-estimate strategy 3 compute budget under HPEPDOCK | §4 budget table → 12–20 h |

End of v1 (post-R1 corrector pass). Awaiting R2 review.
