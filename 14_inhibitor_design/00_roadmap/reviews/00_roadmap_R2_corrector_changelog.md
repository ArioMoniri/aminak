# R2 Corrector Changelog (v1 → v2)

> Maps every R2 reviewer finding to the v2 edits that close it. Format: **R2 finding** → *what changed* → *file/section*.

## R2 HIGH findings (3) — all CLOSED

### HIGH-1 — CID verification is a no-op
- **What R2 said.** The "to verify A0.2" placeholders and the literature InChIKey reference table in v1 §A2 were not actually verified by the corrector — they deferred verification to pipeline runtime, but with no hard-coded ground truth runtime cannot fail an assertion.
- **What v2 does.** Corrector performed all PubChem PUG REST queries *now* (2026-05-18), wrote the per-anchor record to `anchor_compounds_verified.json` with verified CID, InChIKey, formula, MW, IUPAC name, canonical SMILES, and an explicit `wrong_cids_rejected` audit list for each anchor.
  - **Eight of ten v0/v1 CIDs returned the wrong molecule.** Examples: dUMP `22848` → Solanum steroid; 5-FdUMP `15718` → tert-butylphenoxy acid; BrdUMP `135398598` → GTP; nolatrexed `60198` → estrogen steroid; plevitrexed `122478` → imidazole dimer.
  - The corrected ground-truth CIDs are: dUMP 65063, 5-FdUMP 8642, BrdUMP 93036, floxuridine 5790, 5-FU 3385, MTX 126941, raltitrexed 104758, pemetrexed-S 135410875, plevitrexed 135430970, nolatrexed 135400184, ibuprofen 3672.
- **Sections edited.** §1.1 (anchor table rewritten with verified CIDs + InChIKeys), §A2 (table replaced with verified InChIKeys + runtime script that compares against `anchor_compounds_verified.json` as ground truth).

### HIGH-2 — PROLIF cannot flag missing waters
- **What R2 said.** v1 claimed "PROLIF post-analysis explicitly flags whether the top pose would have made a water-bridged H-bond". PROLIF analyses the *waterless* docked complex and cannot, by construction, flag missing crystal waters.
- **What v2 does.** Added a *separate* step **E1b — Crystal water-bridge check** (Strategy 1 only). MDAnalysis-based: align docked pose to 1HVY chain A by Cα, iterate over `resname HOH and around 4.5 (resname UMP or resname DUR)`, apply dual 3.5 Å criterion (water within 3.5 Å of Tyr258:OH AND within 3.5 Å of any pose heavy atom); emit per-pose `_waterbridge.json` and a `water_bridge_lost` column. PROLIF is no longer responsible for the missing-water flag.
- **Sections edited.** §0 "Crystal-water handling" (one-line cross-reference), §E.2 (deferred to E1b), new §E1b with the full algorithm and example code.

### HIGH-3 — HPEPDOCK has no fallback, no timeout, no offline detection
- **What R2 said.** v1's 50–80 h compute budget assumes HPEPDOCK always responds; in practice the web service can queue, 503, or be unreachable for days, making the budget unbounded.
- **What v2 does.** Concrete failure-handling envelope in §D:
  - Per-submission HTTP timeout 60 s.
  - Per-job poll interval 5 min, total wall-time TIMEOUT 4 h.
  - Retries 3× on 5xx with 30-min exponential backoff.
  - TIMEOUT marks `hpepdock_status = TIMEOUT` (not silently dropped).
  - HPEPDOCK unreachable > 15 min cumulative → switch to **CABS-dock** (Kurcinski 2015) for remaining jobs, mark `engine = CABS-dock`.
  - CABS-dock also unreachable → abort peptide track under **extended Stop Condition S3**.
  - > 50 % TIMEOUT triggers full peptide-track abort.
  - HPEPDOCK A0 re-dock control: LR octapeptide into 1HVY, ≤ 4 Å pose RMSD pass criterion.
  - Scrambled-sequence control: numpy seed=42 permutation, submitted as separate job, `canonical_top1 ≤ shuffled_top1 − 2 kcal/mol` else flag `peptide_specificity_unreliable = True`.
- **Sections edited.** §D "Strategy 3 docking deviation" — fully rewritten with envelope table and CABS-dock fallback explicitly named.

## R2 MEDIUM findings (3) — all CLOSED

### MEDIUM-1 — Cofactor-site box centre ambiguity for apo docking
- **What R2 said.** The apo receptor has no D16; the v1 instruction "centroid of raltitrexed D16 (chain A, recomputed from `cofactor_A.pdbqt`)" needed to make explicit that the centroid comes from the *holo* receptor's D16 even when docking against the apo receptor.
- **What v2 does.** §A pocket-definition table row for Strategy 2 now reads "*computed once from the holo receptor `06f_receptor_fixed/cofactor_A.pdbqt` and reused unchanged for both apo and holo dockings — the apo receptor has no D16 by definition*".
- **Sections edited.** §A pocket-definition table.

### MEDIUM-2 — PLIP install plan for arm64-darwin Py 3.14
- **What R2 said.** PLIP depends on the OpenBabel Python bindings; on arm64-darwin Python 3.14 these are not installable via pip alone, requiring a Homebrew openbabel + manual symlink. v1's "pip install plip" line was inadequate.
- **What v2 does.** §F "Tooling decisions" PLIP entry replaced with a four-step install gate (pip install, verify `from openbabel import openbabel` or `import openbabel`, verify `from plip.structure.preparation import PDBComplex`); on failure writes `14_inhibitor_design/00_roadmap/PLIP_STATUS = skip_plip` and the pipeline blanks PLIP columns with `null_reason = "plip_unavailable_arm64_py314_openbabel_binding"` per Stop Condition S3. PROLIF is explicitly noted to be independent of this gate (RDKit-based).
- **Sections edited.** §F "Tooling decisions" PLIP bullet.

### MEDIUM-3 — HPEPDOCK scrambled-sequence control underdefined
- **What R2 said.** v1 mentioned "scrambled-sequence control" in §A0 but did not specify the shuffling algorithm or the pass/fail threshold.
- **What v2 does.** Concrete spec in §D: `numpy.random.default_rng(seed=42).permutation(list(peptide_sequence))`; submit as separate HPEPDOCK job; require `canonical_top1 ≤ shuffled_top1 − 2 kcal/mol`; else `peptide_specificity_unreliable = True`.
- **Sections edited.** §D Strategy 3 deviation block — scrambled-seq control subsection.

## Items the corrector declined or scaled back

None for R2. All six R2 findings (3 HIGH + 3 MEDIUM) were closed with concrete edits; no findings were declined.

## v2 file additions

- `00_roadmap/anchor_compounds_verified.json` — ground-truth CID + InChIKey + formula + MW + SMILES + rejected-CID audit list (11 anchors).
- `00_roadmap/ROADMAP_v1.md` — preserved v1 audit trail.
- `00_roadmap/reviews/00_roadmap_R2_corrector_changelog.md` — this file.

## v2 cross-cutting changes

- Header version bumped v1 → v2.
- Status paragraph rewritten to summarise the v2 closures.
- §1.1 "v0/v1 wrong CIDs rejected" subsection added (audit-trail surface of the eight wrong CIDs).
- Compute-budget table unchanged (R2 did not flag the 50–80 h estimate); the HPEPDOCK envelope adds at most ~1 h of polling overhead (5-min polls × 12 jobs × 1 retry).

## Awaiting

R3 reviewer audit. If R3 signs PASS, execution begins on Strategy 1.
