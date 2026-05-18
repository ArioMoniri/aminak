# Roadmap Review — Round 3 (biologist + bioinformatician)

Reviewer agent role: structural bioinformatician + medicinal-chemistry biologist. Severity: HIGH / MEDIUM / LOW. Round 3, against `ROADMAP.md` v2 (post-R2 corrector pass).

## Overall verdict

**CONDITIONAL_PASS** — one HIGH-severity hole remains (the E1b alignment script does not work on Vina pose PDBQT inputs as written), plus two MEDIUM and a couple of LOWs. Five of the six R2 sign-off items are genuinely closed. The roadmap is *almost* execution-ready: if the corrector fixes E1b's input-format assumption and pre-bakes the LR peptide sequence + the PubChem similarity REST endpoint, R4 should be PASS. We are NOT in cosmetic-churn territory yet but we are close.

## Status of R2 sign-off items (6 items)

**HIGH-1 (CID verification)** — **CLOSED with caveats**. The placeholder text is gone. `anchor_compounds_verified.json` exists and the v2 §A2 script compares against it as ground truth. The eight rejected wrong-CID audit trail is excellent and exactly what R2 demanded. *Caveats baked-in for the v2 audit below*: three rows (pemetrexed, plevitrexed, nolatrexed) still carry `inchikey: null` and rely on a RDKit-derived runtime InChIKey — that is acceptable iff the SMILES used for derivation is itself trusted (see "Verified-anchors JSON audit" below). One Tier-1 row (pemetrexed CID 135410875) is questionable — see NEW-issue [HIGH] below.

**HIGH-2 (PROLIF water-bridge)** — **CLOSED** in the conceptual sense (PROLIF is no longer asked to do impossible work; a separate E1b script was added) but the script itself has a NEW HIGH-severity bug. See "New v2 issues" — the script loads a `.pdbqt` file with `mda.Universe(pose_pdbqt)` and then `align.alignto(... select="protein and name CA and segid A")`. The Vina pose PDBQT contains *only the ligand*, no protein, no segid, no Cα. The alignment selection will return an empty AtomGroup and MDAnalysis will raise. The water-bridge gate cannot run as written.

**HIGH-3 (HPEPDOCK fallback)** — **CLOSED**. The §D envelope table (60 s POST timeout / 5-min polls / 4 h per-job TIMEOUT / 3× retry with exponential backoff / CABS-dock named fallback / > 50 % TIMEOUT triggers strategy abort / extended S3 for service-down) is a clean specification. The scrambled-sequence control (numpy seed=42 permutation, ≥ 2 kcal/mol separation) is now concrete. Two residual concerns flagged below: poll-count math, and CABS-dock-as-2026-fallback viability.

**MEDIUM-1 (cofactor box centre)** — **CLOSED**. §A row 2 now reads "computed once from the holo receptor `cofactor_A.pdbqt` and reused unchanged for both apo and holo dockings — the apo receptor has no D16 by definition." One sentence; exactly the fix R2 asked for.

**MEDIUM-2 (PLIP install gate)** — **CLOSED**. The four-step bash gate (pip → `from openbabel import openbabel` → fallback `import openbabel` → final `from plip.structure.preparation import PDBComplex`) is the right shape; `PLIP_STATUS = skip_plip` writes the sentinel; PROLIF is explicitly noted as RDKit-independent of this gate. Acceptable.

**MEDIUM-3 (HPEPDOCK scrambled control)** — **CLOSED**. Concrete spec in §D: `numpy.random.default_rng(seed=42).permutation(...)`; submitted as separate job; `canonical_top1 ≤ shuffled_top1 − 2 kcal/mol` else `peptide_specificity_unreliable = True`. One residual LOW: the R2 reviewer asked whether the Cys position should be preserved (Cys is the only nucleophile in LSCQLYQR — randomising it changes chemistry meaningfully). v2 does not address this; impact is bounded because it's an internal control, not a docking-mode question.

## Verified-anchors JSON audit

Per-anchor verdict against my own knowledge + the published literature:

- **dUMP CID 65063 / `JSRLJPSBLDHEIO-SHYZEUOFSA-N` / C9H13N2O8P** — **ACCEPT**. This is the canonical 2′-deoxyuridine-5′-monophosphate record. Formula and InChIKey match. v0/v1 CID 22848 was indeed wrong (steroid).
- **5-FdUMP CID 8642 / `HFEKDTCAMMOLQP-RRKCRQDMSA-N` / C9H12FN2O8P** — **ACCEPT**. Carreras & Santi 1995 active species; InChIKey matches the standard PubChem record. R2's note that v0/v1 CID 15718 was wrong is correct.
- **BrdUMP CID 93036 / `LHLHVDBXXZVYJT-RRKCRQDMSA-N` / C9H12BrN2O8P** — **ACCEPT**. Formula matches 5-bromo-dUMP (Santi & McHenry 1972 probe); stereochem layer `RRKCRQDMSA` matches the 2′-deoxyribose configuration of the other dU compounds.
- **Floxuridine CID 5790 / `ODKNJVUHOIMIIZ-RRKCRQDMSA-N` / C9H11FN2O5** — **ACCEPT**. Standard FdUR record.
- **5-FU CID 3385 / `GHASVSINZRGABV-UHFFFAOYSA-N` / C4H3FN2O2** — **ACCEPT**. Universally cited 5-FU InChIKey.
- **MTX CID 126941 / `FBOZXECLQNJBKD-ZDUSSCGKSA-N` / C20H22N8O5** — **ACCEPT**. (S)-methotrexate standard.
- **Raltitrexed CID 104758 / `IVTVGDXNLFLDRM-HNNXBMFYSA-N` / C21H22N4O6S** — **ACCEPT**. (S)-Tomudex / D16 — formula and InChIKey both verified. The 104758-vs-135400182 tautomer note is appropriate.
- **Pemetrexed (S) CID 135410875 / `inchikey: null` / C20H21N5O6** — **REVISE (HIGH)**. The widely cited PubChem CID for pemetrexed (Alimta, free acid, S configuration) is **CID 446556**, not 135410875. The v1 reviewer originally suggested CID 60843 (R-enantiomer per the v2 corrector). The v2 corrector now picks 135410875 but provides *no* InChIKey to verify it. This is the one Tier-1 anchor whose identity I am not confident in based on the JSON alone. Concrete fix in NEW-issue [HIGH] below.
- **Nolatrexed CID 135400184 / `inchikey: null` / C14H12N4OS** — **ACCEPT-WITH-CAVEAT**. Formula matches the published nolatrexed (Webber 1996, AG-337, 2-amino-6-methyl-5-(4-pyridylthio)-3H-quinazolin-4-one) and the IUPAC name in the JSON is right. v0/v1 CID 60198 → estrogen steroid is a real and damning catch. CAVEAT: the canonical PubChem CID for nolatrexed has historically been **CID 6604755** in some references; the corrector should at minimum note why 135400184 was chosen over 6604755 (it may be a tautomer choice). Not a blocker — formula + IUPAC are right.
- **Plevitrexed (ZD9331) CID 135430970 / `inchikey: null` / C26H25FN8O4** — **ACCEPT**. Formula C26H25FN8O4 matches Jackman 1997 BGC9331 with the fluorobenzoyl + propynyl + 2,7-dimethylquinazolinone framework. IUPAC name in the JSON is consistent with the literature structure.
- **Ibuprofen CID 3672 / `HEFNNWSXXWATRW-UHFFFAOYSA-N` / C13H18O2** — **ACCEPT**.
- **LR octapeptide** — **ACCEPT-WITH-LOW**. JSON entry notes "Build from sequence LSCQLYQR via RDKit `Chem.MolFromSequence`". Good that the sequence is now written down (which closes R2 sign-off #6 on the corrector side). But the roadmap §A3 still says the human-TYMS residue range is "computed at runtime via Biopython pairwise alignment to Cardinale 2011 *E. coli*"; the *E. coli* source range is not given anywhere. R2 sign-off #6 is technically half-closed.

## New v2 issues (in addition to R2 carryover)

- **[HIGH] E1b water-bridge script cannot align a Vina pose PDBQT.** The v2 script does `pose = mda.Universe(pose_pdbqt)` and then `align.alignto(pose, xtal, select="protein and name CA and segid A")`. A standard Vina output PDBQT contains the *ligand only*, no protein atoms, no chain identifier, no Cα atoms. The selection returns zero atoms; MDAnalysis raises `SelectionError`. Per TECHNICAL_NOTES Phase 7 §"Multi-replica Vina", Vina poses are ligand-only. Fix: load the receptor PDBQT separately, align *receptor → receptor* by Cα, then apply the transformation matrix to the ligand pose. Or simpler: load `xtal` (1HVY full structure), do not align at all (the docking box is in the same 1HVY frame as the crystal because the receptor was prepared from 1HVY without re-centring — verify against TECHNICAL_NOTES Phase 6c), and just check the ligand pose coordinates directly against the crystal water positions. Either approach works; the as-written script does not.

- **[HIGH] Pemetrexed CID 135410875 is unverified.** See JSON audit above. The v2 corrector chose this CID over the conventional CID 60843 or CID 446556 but provided no InChIKey to cross-check. If A2 runs, `Chem.MolToInchiKey(Chem.MolFromSmiles(pubchem.ConnectivitySMILES))` will return *some* InChIKey, but with `expected_key = None` the A2 gate (`if expected_key and actual_key != expected_key`) is a no-op for this row — same failure mode R2 flagged for v1. Fix: bake the literature InChIKey for (S)-pemetrexed (`QOFFJEBXNKRSPX-ZDUSSCGKSA-N` per the R1 reviewer's table) into the JSON `inchikey` field, and add a comment explaining which CID returned what when. If 135410875 returns a different InChIKey than the literature value, switch to the CID that does.

- **[MEDIUM] `ConnectivitySMILES` is not a standard PubChem PUG REST property name.** The standard property names are `CanonicalSMILES` and `IsomericSMILES`. `ConnectivitySMILES` is not in the PUG REST property table (PubChem documentation, accessed via the schema at `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/<cid>/property/<property>/JSON`). If the corrector intended the connectivity-only (no stereochem) SMILES, that is a *bad* choice for InChIKey derivation: an InChIKey computed from a stereo-stripped SMILES will not match the stereo-bearing literature InChIKey. The right field is **`IsomericSMILES`** (preserves @, /\\, E/Z, R/S). Fix in §A2: replace `ConnectivitySMILES` → `IsomericSMILES` in both the URL path and the dict access. This is a one-token fix but it matters: the three null-InChIKey rows (pemetrexed, plevitrexed, nolatrexed) all have defined stereocentres or geometric isomers.

- **[MEDIUM] HPEPDOCK 4 h timeout × 5-min polls × 13 jobs = 624 polls; CABS-dock 2026 viability not verified.** Math: 4 h / 5 min = 48 polls per job worst-case; 13 peptide jobs × 48 = 624 GET requests against the Zhou lab service in worst case (plus 13 more for scrambled-controls = 1248). At normal queue depth most jobs finish in 30–60 min, so the realistic poll count is ~5–10 per job (65–130 total), which is fine. The 624-poll worst-case is bounded and the abort-at-50%-TIMEOUT condition will fire before then. **However**: CABS-dock as the named fallback is itself a web service operated by the Kolinski lab; its uptime in 2026 is unverified by v2. If both HPEPDOCK *and* CABS-dock are down on execution date, the §D chain dead-ends to "abort under S3" which is the correct behaviour, but the roadmap should explicitly say "a fresh agent must `curl -I https://huanglab.phys.hust.edu.cn/hpepdock/` and `curl -I http://biocomp.chem.uw.edu.pl/CABSdock` *before* submitting the first job and abort early if both 5xx." One-paragraph fix in §D.

- **[MEDIUM] LR-peptide sequence is in the JSON but not in the roadmap §A3 prose.** R2 sign-off #6 was "write the Cardinale 2011 *E. coli* source peptide sequence into §A3 so the Biopython pairwise alignment is reproducible". The JSON now contains `LSCQLYQR` as the human-TYMS sequence to build, but §A3 still refers vaguely to "the Cardinale 2011 *E. coli* peptide" without writing it down. A fresh execution agent reading only the roadmap (not the JSON) cannot reproduce the alignment. Fix: paste the *E. coli* source sequence (Cardinale 2011 reports it explicitly in their supplementary) into §A3 prose.

- **[LOW] DUD-E web upload + HPEPDOCK web upload procedure not documented.** §1.2 says "DUD-E web service" and §D says "HPEPDOCK web service" but neither gives the endpoint URL, upload format (SMILES list? mol2? PDB? FASTA?), or result-retrieval mechanism (HTML scrape? REST poll? email-based?). A fresh agent will Google and may pick the wrong service version. Same for the **PubChem similarity REST endpoint** (§1.2): the actual endpoint is `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/<cid>/cids/JSON?Threshold=70` — write it down.

- **[LOW] DUD-E rate limit "~10 jobs/day, 3–24 h turnaround" is asserted without source.** Cite the DUD-E web FAQ (or a recent issue tracker entry) so a fresh agent can verify before submission.

## Execution-readiness verdict

**NO, not yet — but only barely.** A fresh agent armed with the roadmap + the repo + Vina/HPEPDOCK/RDKit/MDAnalysis/freesasa/FPocket/PROLIF/PLIP could execute Strategies 1, 2, and 4 essentially as written. The three blockers for full readiness are:

1. **E1b water-bridge script** crashes on first invocation because Vina pose PDBQT has no protein atoms to align by Cα (NEW-HIGH). Strategy 1 cannot complete its post-analysis as specified.
2. **Pemetrexed CID 135410875** identity is unverified; the A2 gate will silently pass for this row because `expected_key = null` (NEW-HIGH).
3. **`ConnectivitySMILES` field name** is wrong; A2 will fail with `KeyError` on the dict access for the three null-InChIKey rows (NEW-MEDIUM).

The LR-peptide sequence + DUD-E endpoint + HPEPDOCK endpoint + PubChem similarity endpoint omissions are LOW-cost rewordings that any competent agent can resolve by reading the JSON and one minute of documentation; they would not by themselves block execution but they would slow it.

## Sign-off requirements for R3 → R4

Corrector agent must, before R4 PASS:

1. **Fix the E1b script.** Either load the receptor PDB(QT) explicitly and do the Cα alignment on the receptor, then apply the transform to the pose; or accept that 1HVY-frame Vina runs need no alignment and check pose-vs-water distances directly. Verify the chosen path actually runs on a representative `*_apo_seed42.pdbqt` file.

2. **Resolve pemetrexed.** Bake the literature InChIKey for (S)-pemetrexed (free acid) into `anchor_compounds_verified.json`. If PubChem CID 135410875 returns the expected InChIKey, keep it; otherwise switch the CID. Add a one-line note explaining the choice between 135410875 / 446556 / 60843.

3. **Replace `ConnectivitySMILES` with `IsomericSMILES`** in §A2 (both the URL path and the dict access). This is a one-token edit.

4. **Bake the LR-peptide *E. coli* source sequence** into §A3 prose (not just the JSON).

5. **Write down the four web endpoints** explicitly: DUD-E URL + submission format, HPEPDOCK URL + submission format + result-retrieval, PubChem similarity REST endpoint with `Threshold=70`, CABS-dock URL. One paragraph in §1.2 + §D suffices.

6. **(LOW)** One sentence in §D pre-flighting HPEPDOCK + CABS-dock liveness with `curl -I` before the first submission, aborting early if both are down rather than burning 4 h to TIMEOUT.

If items 1–3 land cleanly, items 4–6 are minor and R4 should PASS. We are one corrector pass away from execution.

End of Round 3 review.
