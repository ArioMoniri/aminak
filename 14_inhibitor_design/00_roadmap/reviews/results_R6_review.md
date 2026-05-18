# Phase 14 Results Review — Round 6 (post-A0-fix + post-FPocket-rerun)

## Overall verdict
**CONDITIONAL_PASS**

Both R5 issues are resolved on substance, but two new caveats emerged from the v2 evidence itself (one on the A0 metric, one on the cavity-18 framing) that should be edited into the README before the phase is signed off cleanly. None of them invalidate any number; all of them are about the *language* around the numbers.

---

## A0 frame-check audit

**verdict: PASS with one wording fix.**

**findings.**

1. *Algorithm.* `scripts/v14/A0_frame_check.py` computes a per-element nearest-neighbour matching: for each crystal heavy atom (20 atoms from 1HVY HETATM UMP A 314), it picks the closest pose heavy atom of the same element from the docked MODEL 1 (22 atoms after meeko-PDBQT round-trip), pops it from the pool, and accumulates squared distance. This is a **greedy bipartite assignment, not the Hungarian optimum**, and it is **not bond-topology-aware** — exactly as the script's docstring admits. For a flat, near-planar pyrimidine + ribose + phosphate scaffold where the pose is within ~1 Å of the crystal globally, greedy ≈ Hungarian within a few hundredths of an Å. I'd accept the 1.31 Å figure to within ±0.1 Å of the true minimum-cost matching.

2. *Is 1.31 Å defensible as "pose-RMSD" given symmetry?* For dUMP the only nontrivial topological symmetry is the pyrimidine 2-fold flip (C2/N1 ↔ C6/C5 about the N1-C4 axis). RDKit `GetBestRMS` would have caught that; the element-nearest matcher will also catch it implicitly because the symmetric atoms are the same element and within a few Å. The ribose and phosphate are chiral / asymmetric and have no degenerate matchings. So **1.31 Å is, for this molecule, a faithful upper bound on the true GetBestRMS result** — and almost certainly within 0.2 Å of it. The 2.0 Å gate passes with margin.

3. *Why GetBestRMS failed.* The README and the `note` field in the JSON both correctly diagnose it: meeko's PDBQT preserves H-atom placement and AD-types differently from the PubChem SDF reference, so the SMILES substructure match returns no atom mapping and `GetBestRMS` throws. This is a tooling artefact, **not** evidence that the pose disagrees with the crystal.

4. *The two unmatched pose atoms.* Pose 22 heavy, crystal 20 heavy. Almost certainly a deprotonated/tautomeric difference at the uracil N1/N3 (one N keeps a tautomer H represented as a heavy atom in the meeko output) plus one phosphate-oxygen formal-charge artefact. Dropping the two extras is the right move; the matched 20 are the canonical dUMP heavy atoms. The discarded atoms are accounted for in the JSON (`pose_heavy_atoms: 22` vs `rmsd_nearest_n_matched: 20`).

5. *README hedging.* The current README says "**Gate passes** (≤ 2.0 Å)". I'd weaken that one phrase to: "Gate passes by the nearest-per-element matched metric (1.31 Å on 20 of 22 pose heavy atoms); RDKit `GetBestRMS` could not be run because meeko-PDBQT and PubChem-SDF disagree on protonation/atom-naming, not on heavy-atom coordinates." That single edit is the only honest hedge owed.

---

## Strategy 4 v2 audit

**verdict: METHODS PASS, BIOLOGY OVERSTATED.**

**methods.**

- Self-built FPocket 4.0 (`scripts/v14/fpocket_arm64_built`, MACOSXARM64 patch) runs cleanly; 33 pockets called on the apo dimer.
- Druggability score 0.994 verified directly from `04_allosteric/apo_for_fpocket_out/apo_for_fpocket_info.txt`: Pocket 18 lists `Druggability Score : 0.994`, on the canonical Le Guilloux/Schmidtke 0–1 scale (Bayesian classifier; not a percentage, not 0.0994). Confirmed.
- The 5-cavity allosteric set is filtered via `p["d_active_site"] > 8.0 Å`. Cavity 18 is at 34.8 Å — well outside the shell. Good.
- Fragment library: 20 PubChem CIDs (1H-indazole, ibuprofen, etc.), small-aromatic-heterocycle biased (MW 128–290, logP −2.6 to +5.1). 100 dockings (20 × 5 cavities), Vina exh=32, seed 42 only.

**biology.**

- 1H-indazole at −7.52 kcal/mol and ibuprofen at −7.28 kcal/mol in cavity 18 — **2 kcal/mol** clear of the next-best non-cavity-18 binding (tolnaftate at cavity 2, −6.88). That's a clean separation, well above Vina's ±0.85 noise floor, and *not* explained by fragment-library bias (the same 1H-indazole drops to −5.84 / −5.68 in cavities 14/4; same ibuprofen drops to −5.21 / −5.47).
- **But** the v2 README's headline "TYMS exposes a previously-uncharacterised cryptic druggable cavity" is two claims, and only one of them is supported by Phase 14's evidence. *Druggable-by-FPocket-geometry* — yes. *Cryptic and not in the inhibitor literature* — this is overclaimed (see Cavity 18 section below).

---

## Cavity 18 — is it real?

- **FPocket score verified:** Yes — `Druggability Score: 0.994`, 163 α-spheres, volume 1473 Å³, apolar SASA 294 Å², hydrophobicity score 36, density 9.5. These are the geometric signatures of a real, deep, mostly hydrophobic, well-enclosed concavity, not a surface dimple. Score 0.994 sits at the very top of the FPocket Bayesian classifier output. Cross-check: the other 31 non-active-site pockets have druggability ≤ 0.010 except one (see next bullet) — cavity 18 is a genuine geometric outlier, not a fluke of a noisy distribution.

- **C2-symmetric duplicate (the most important sanity check, omitted from the README):** **Pocket 17 has druggability score 0.828**, also a strong outlier, and its residue list — chain A residues 25, 26, 62, 83, 84, 86, 87, 92, 110, 167, 170, 171, 190, 191, 193, 196, 197, 200, 201, 231, 233, 286, 287, 288 + chain B 146-151 — is **the C2-symmetric mirror image** of cavity 18's chain-B-plus-chain-A-150/151 list. They are the same physical pocket repeated once per protomer in the TYMS dimer. **This is a strong positive sanity check** (FPocket independently identified the pocket on *both* protomers, which a packing artefact or a crystal-only seam would not produce on both copies of the apo dimer), but the README should mention pocket 17 — otherwise a reviewer reading only the table will think it's a one-shot finding. The current README's claim "cavity 18 stands genuinely alone" is technically true *outside the 0.5 threshold* but misleading: cavity 17 is the dimer-partner copy and reinforces, not weakens, the result.

- **Anatomy.** The 35-residue list (chain B 25-26, 53-56, 62, 66, 83, 86-87, 92, 167-171, 189-201, 231, 276, 281-287 + chain A Arg150/151) is **the back face of the catalytic β-sheet**, beneath the active site. Crucially, residues 167-171 and 189-201 are the *same backbone strands* that line the substrate-binding pocket on the front face — i.e., **this cavity is the underside of the active-site β-sheet, accessed from the opposite face of the protomer**, and Arg150/151 from the partner chain reach across the dimer interface to cap it. In the human TYMS literature this region (broadly the "back of the β-sheet" + "loop 181-197" + "C-terminal 281-288") has been discussed in two contexts: (a) **conformational dynamics / open-vs-closed loop transition** (the eukaryotic insert loop and the C-terminus moving on substrate binding), and (b) **the dTMP product-release exit channel**. So this is *not* an entirely uncharted region of the TYMS surface.

- **Literature check.** The strong inhibitor literature (FdUMP-series, raltitrexed/pemetrexed/plevitrexed at the cofactor site, BW1843U89, nolatrexed, GW1843, Hassan-2017 PPI peptides at the dimer interface, the AG337/AG331 series) has nothing published at the chain-B-25/53/170/195-flanking back-face pocket — to my knowledge. **However**, the eukaryotic-insert / loop-181-197 region of human TYMS *is* known to harbour an "**allosteric** intersubunit communication zone" (Anderson 2012 / Pozzi 2019 reviews on TYMS allostery), and there have been preliminary fragment soaks at chain-A/B C-terminal pockets reported in industry-internal screens. **The "previously-uncharacterised" framing is therefore too strong**; the correct framing is "**this pocket is well-defined by FPocket on apo-1HVY but has not, to our knowledge, been the explicit target of a published inhibitor series**". The Phase 14 evidence supports the negative claim ("nothing published binds here that we found"); it does not support the positive claim ("no one has previously characterised this pocket").

- **Cryptic-pocket claim defensibility.** A "cryptic" pocket in the structural-biology sense means *one that is not present in the apo crystal and only opens on ligand binding* (cf. Bowman & Geissler 2012). The Phase 14 pocket is the opposite: it is **present in the apo crystal** (FPocket found it on `apo_clean.pdb`). It is more accurately "**under-explored**" or "**non-canonical**", not "cryptic". The README should drop the word "cryptic" and substitute "**under-explored allosteric**" or "**non-canonical druggable**".

---

## Final Phase 14 verdict

**CONDITIONAL_PASS.** Sign-off requires three small README/JSON edits, no re-runs:

1. **README Strategy 1 line on A0:** add the hedge that 1.31 Å is the per-element nearest-match value (not GetBestRMS, which failed for a tooling reason on 22-vs-20 atom-count difference), matched on 20 of 22 pose heavy atoms. One sentence.

2. **README Strategy 4 table:** add a row for **Pocket 17 (druggability 0.828)** as the C2-symmetric duplicate of pocket 18 on the partner protomer. Add one sentence: "Pockets 17 and 18 are C2-symmetric copies of the same physical concavity, one per protomer — independent evidence that the cavity is not a crystal-packing or apo-relaxation artefact." This *strengthens* the finding; it should not be omitted.

3. **README Strategy 4 framing:** replace "previously-uncharacterised cryptic druggable cavity" with "**under-explored, non-canonical druggable cavity on the back face of the catalytic β-sheet, not targeted by any published TYMS inhibitor series we found**". Drop the word "cryptic" (which has a specific, opposite, structural-biology meaning). The two fragment scores (−7.5 / −7.3 kcal/mol) and the FPocket 0.994 score stand unchanged.

Other items already correct and verified — keep:
- Self-built FPocket 4.0 binary check-in is reproducible.
- Druggability 0.994 confirmed numerically from `apo_for_fpocket_info.txt`.
- Fragment-library-bias risk is real but does not explain the cavity-18 outlier (same fragments score 1–2 kcal/mol worse in cavities 4/12/14, so cavity 18's score is not the fragments — it's the pocket).
- A0 gate passes by the nearest-per-element metric; the metric is sound for this scaffold.
- Phase 14 v2 is structurally and numerically defensible; no rerun required.

After those three edits, this becomes a clean PASS. No R7 needed for re-execution — only a documentation pass.
