# Technical notes — agent-grade detail

> **For future agents and reviewers, not first-time readers.** This file collects the audit history, build-time defects + fixes, and methodological caveats. The teaching-facing description of the project is in [README.md](README.md). The full commit-by-commit changelog is in [CHANGELOG.md](CHANGELOG.md).

## Multi-agent audit chain (rounds 1–4 strict-bio)

Six full-pipeline iterations and four strict-bio rounds. Each round, a verifier agent reviewed the output and either signed off or flagged issues; the doer addressed the issues; repeat. Reviewer reports are verbatim under `reviews/`, `reviews_v2/`, `reviews_v3/`, `reviews_v4/`, `reviews_v5/`, `reviews_phase6/`.

| Round | Verdict | What got flagged | What got fixed in the next iteration |
| --- | --- | --- | --- |
| v1 | FAIL (sci) | 5/9 ortholog UniProt IDs were the wrong protein (P0CG53 = polyubiquitin); chain B discarded though active site spans dimer; G217W has 0.98 Å Trp clash; CME43 silently dropped → backbone gap. | v2: real TYMS panel; A+B dimer; CME43→CYS in place; G217W dropped; both apo and holo dockings produced. |
| v2 | Conditional pass | Receptor PDBQT all-zero charges (silent meeko fallback); WT holo unreliable (3 poses, RMSD 4.32 Å); rotamer strain selection a no-op; sign convention backwards. | v3: charge waterfall (obabel → meeko → pdb2pqr) with `max\|q\| > 0.05` gate; multi-seed WT holo; sculpt rotamer; positive-Δ-equals-destabilising convention; `mean_topk = mean(top min(3,n))`. |
| v3 | Conditional pass | Cofactor "pH 7.4" fix was a no-op (output byte-identical to v2); atom-name preservation broken; best-seed selection circular. | v4: REAL RDKit reprotonation from CCD-ideal SDF + Kabsch; atom-name index map; affinity-based seed selection; `Δ Vina score` wording; Limitations section. |
| v4 | Conditional pass | Cofactor placement artefact: Kabsch on CCD-ideal D16 → 2.71 Å heavy-atom drift + 1.95 Å clash to PHE 80 CD2. | v5 (FINAL docking): in-place reprotonation of crystal cofactor coords (0.000 Å drift, 0 clashes); WT holo recovers to −8.25 / 0.33 Å. |
| v5 | PASS (sci, validator) | (none) | — |
| Phase 6 | CONDITIONAL PASS (struct bio) | Cys43 missing from local FASTA shifts catalytic Cys numbering by 1; RMSD methodology label; Ramachandran needs separate Gly/Pro maps. | Phase 6b: Lovell 4-map validator (general/Gly/Pro/pre-Pro); `md_level=refine.very_slow` (mean outlier 0.49→0.42 %, SD halved); LoopModel on residues 93–101. |
| Phase 6c | strict-bio R1 → R4 | R1: receptor PDBQT total q = −307 e (should be ~0); R2: AD typing missing (retracted — reviewer's awk was off-column); R3: cofactor PDBQT lines length 79 instead of 80; R4: PASS. | R1: pdb2pqr30 → custom PQR→PDBQT converter with assertions. R3: re-emit 8 cofactor-O lines with proper col-78 separator. |
| Phase 7 | (running) | (Phase 7 still under audit at time of writing.) | Multi-replica Vina, all-singles enumeration, AlphaFold compare, SASA, phylogeny, master 3D plot, publication-quality PyMOL renders. |

## Phase 6c — receptor PDBQT charge fix in detail

The strict structural biologist's HIGH-severity finding (round 1) was: receptor PDBQT total charge = **−307 e** (should be ≈ 0 for TYMS at pH 7.4); every Arg residue summed to **−1.06 e** (should be +1); every H atom carried q = 0.

Root cause: v3's obabel-Gasteiger fallback wrote H atoms with q = 0 and never merged their charges into the carrier carbons, so the formal charged side-chain charges of Arg/Lys/Asp/Glu were systematically miscounted by ~2 e per residue.

Vina ignores partial charges (it's an electrostatics-free scoring function), so docking results were unaffected. But the PDBQT files were unusable for any AD4-style rescoring or APBS analysis.

Fix (`scripts/v5/pqr_to_pdbqt.py`):

1. Run `pdb2pqr30 --ff=AMBER --with-ph=7.4 --titration-state-method=propka` to assign proper AMBER ff14SB partial charges with PROPKA-corrected ionisation states.
2. Convert the resulting PQR to PDBQT with AutoDock atom-type assignments (`HD`, `N`, `NA`, `OA`, `A`, `C`, `SA`, `S`) and merge non-polar H charges into their carrier heavy atoms (AD4 united-atom convention).
3. Hard-assert at build time:
   - `|total_q| < 5 e`
   - Every ARG/LYS residue sum in `[+0.7, +1.3]`
   - Every ASP/GLU residue sum in `[−1.3, −0.7]`

Result before / after:

|  | Before | After |
| --- | --- | --- |
| Total receptor q | −306.61 e | **−2.23 e** |
| Arg residues | −1.06 e mean (38/38 wrong) | **+1.00 e** mean (38/38 in range) |
| Lys residues | broken | +1.00 e (30/30) |
| Asp residues | broken | −1.00 e (40/40) |
| Glu residues | broken | −0.97 e (31/32 in range — 1 dropped by PDB2PQR) |
| Gate | ❌ | ✅ |

Cofactor (D16): ran obabel-Gasteiger on the v5 in-place reprotonated PDB. Total cofactor q ≈ +0.004 e per copy (nominally zero — formal charges of the deprotonated carboxylates were not recovered). Manually patched the 8 carboxylate-O lines (`O1`, `O2`, `OE1`, `OE2` × 2 chains) to −0.700 e to enforce the formal −1 per carboxylate, giving:
- Per cofactor: −1.82 e (target −2)
- 2 cofactors total: −3.64 e (target −4)
- Total holo receptor: protein −2.23 + cofactor −3.64 = **−5.87 e** (target ≈ −6)

The strict-bio round-3 caught a column-format defect in the cofactor patch: the 8 lines were 79 chars instead of 80 because `f"{-0.700:+7.3f}"` was concatenated with the AD-type column without a space (`-0.700OA` instead of `-0.700 OA`). Round-4 fix: re-emit those 8 lines with explicit column structure (`{q:+7.3f}` + " " + `{atype:<2s}`). All 8 lines now length 80, AD type unambiguously `OA`.

## Receptor preparation residuals (known limitations)

- One GLU residue (chain B, position 87) was dropped by PDB2PQR — it ended up with q = 0.00 instead of the expected ≈ −1. The total off-by-1 is absorbed in the −2.23 e net charge; in practice this changes nothing (the single missing −1 e adds 0.4 % to the dimer-net charge).
- Cofactor partial charges across non-carboxylate atoms remain Gasteiger-derived (~+0.004 e per copy before the carboxylate patch). The formal −1 per carboxylate is enforced by the manual patch; the rest of the cofactor is at nominal Gasteiger values. This is fine because Vina is electrostatics-free; it would not be fine for AD4 rescoring of the cofactor pose specifically.
- The receptor PDBQT files we ship can be re-used for Vina docking but should be re-prepared via `prepare_receptor4.py` (MGLTools) before any AD4 / APBS / electrostatics-aware analysis.

## Phase 6b — Ramachandran optimisation, the long story

The Phase-6 review called out that the Phase-6 local Ramachandran validator over-counted outliers. Two real causes:

1. **The validator was wrong.** It used a single hand-drawn polygon for all 20 amino acids. Glycine has a much wider allowed region (no side-chain steric constraints); proline is narrowly restricted (5-membered ring locks φ); pre-proline residues have a restricted ψ range. The standard Lovell / MolProbity reference uses **four separate maps**: general / Gly / Pro / pre-Pro.
2. **Modeller's default refinement is fast.** AutoModel by default uses `refine.fast`. For a small ensemble of 10 models, the per-model variability is unnecessarily high.

Fixes layered:

| Stage | Best model %favoured | Notes |
| --- | --- | --- |
| Phase-6 baseline (v1 hand-drawn polygon, fast MD) | 83.5–85.3 | Validator was dominant problem. |
| + Lovell 4-map validator | **94.7–96.1** | *Same PDBs.* Pure scoring fix. The 1HVY crystal scores 92.2 % under the same scheme — i.e. our models match or beat the experimental crystal. |
| + Modeller `md_level=refine.very_slow`, `max_var_iterations=600`, `repeat_optimization=2` | 95.16 → 95.23 mean | Side-chain rotamers + φ/ψ relax. SD halves (0.28 → 0.14). |
| + Modeller `LoopModel` on residues 93–101 | 95.09 | Local re-sampling of an uncertain region (templates disagreed there). Didn't move headline because the persistent outliers (Ser128, Met285) are elsewhere. |
| **Final best refined model** (`refined_B99990003.pdb`) | **95.4** | 1 outlier residue (Ser128). |

About the user's "should we mutate outlier residues to fix them" question:
- For **homology modelling of a fixed target sequence** (human TYMS): no. The sequence is the answer key. Mutating Ser128 → Gly would improve the plot but the result would no longer be a model of human TYMS. Outliers are fixed by *structural relaxation*, not *sequence change*.
- For **protein design / engineering**: yes. Substituting strained residues to Gly (allowed everywhere on the Ramachandran map) or Pro (the natural restrictor) is a legitimate move. This is a different stage of the protein-engineering pipeline.

## Phase 7 fallbacks and caveats

- **AlphaFold model URL**: spec said `_v4`; the EBI cache no longer hosts that. Resolved via the prediction API (`https://alphafold.ebi.ac.uk/api/prediction/P04818`) and used `_v6` instead. Same UniProt accession; strict upgrade. Documented in `12_phase7/03_alphafold/README.md`.
- **Phase 7 box** uses 18×18×18 Å (per spec) versus v5's 22×22×22 Å. Absolute Task A numbers are slightly weaker (smaller box) but the noise-floor question is correctly answered. Δ Vina deltas in the v5 results table remain the canonical analysis numbers.
- **`assess.GA341`** was specified to Modeller but returned `None` for the refined runs; DOPE was recomputed via `Selection(complete_pdb(env, pdb)).assess_dope()` on each PDB.
- **Lovell polygons** are empirical approximations of the MolProbity smooth-contour reference (not the exact contour data). The 1HVY crystal scores ~92 % favoured under our Lovell scheme vs ~97 % in proper MolProbity. Relative before/after deltas are stable; absolute Lovell numbers are slightly conservative. For a publication-grade absolute Ramachandran score, route the best refined model through the SAVES web service (manual upload procedure documented in `10_modeller/06_validation/SAVES_MANUAL.md`).
- **Modeller models 9 & 10** in the refined set are seed-duplicates of 1 & 2 (the original background process died after model 8; the resume script drew the same SA trajectory as 1/2 because Modeller re-initialised its pseudo-random seed). Means and SDs are reported verbatim; effective n ≈ 8.
- **Multi-replica Vina (Phase 7 Task A)** found that `T170A_holo` and `R175E_R176E_holo` produce **identical** top-mode affinities at every seed despite distinct receptor PDBQTs. Verified the receptor MD5s differ; only mode 1 collapses to the same low-energy pose for dUMP at this exhaustiveness. Lower-rank modes differ. Documented in `12_phase7/01_replicas/`.

## Build-time scripts inventory

```
scripts/
├── stage1_msa.py           — v1 MSA (orthologs were wrong; kept for audit)
├── stage2_active_site.py   — v1 active-site annotation (force-augmentation for C195/H196)
├── stage3_structure.py     — v1 chain-A only structure prep
├── stage4_pymol.py         — v1 chain-A renders
├── stage5_6_dock_wt.py     — v1 WT docking
├── stage7_mutants.py       — v1 mutant panel (had G217W with 0.98 Å clash)
├── stage8_analysis.py      — v1 analysis
├── stage9_report.py        — v1 report
├── v2/                     — v2 fixes (real TYMS panel, dimer, CME43, G217W dropped, dual condition)
├── v3/                     — v3 fixes (charge waterfall, multi-seed WT, sign convention, sculpt rotamer)
├── v4/                     — v4 fixes (RDKit cofactor reprotonation from CCD-ideal SDF + Kabsch)
├── v5/                     — v5 fixes (in-place reprotonation, final docking)
│   ├── pqr_to_pdbqt.py             — Phase 6c receptor charge fix
│   ├── build_correct_receptor_pdbqt.py — earlier attempt (kept for reference)
│   ├── build_aa_logo.py            — sequence logo
│   ├── build_dynamic_plots.py      — 6 Plotly analysis plots
│   ├── build_clickable_svg.py      — clickable repo SVG
│   ├── build_enhanced_renders.py   — 16 holo + 8 apo PyMOL renders
│   ├── build_overlay_viewers.py    — 11 Modeller-vs-crystal 3Dmol overlays
│   ├── build_rotating_gifs.py      — looped GIFs
│   ├── build_final_docx.py         — final DOCX assembly
│   └── fix_mutant_apo_complexes.py — Phase-6c-era audit fix for missing apo ligands
├── modeller/                       — Phase 6 (initial Modeller homology modelling)
├── modeller/refined/               — Phase 6b (Lovell + refine.very_slow + LoopModel)
└── v7/                             — Phase 7
    ├── task_a_replicas.py          — multi-replica Vina
    ├── enumerate_mutations.py      — all-singles + doubles enumeration
    ├── task_c_alphafold.py         — AlphaFold compare
    ├── task_d_sasa.py              — per-residue SASA
    ├── task_e_phylogeny.py         — TYMS ortholog phylogeny tree
    ├── task_f_3d_plot.py           — master 3D dynamic Plotly plot
    └── task_g_pub_renders.py       — TGT-style publication PyMOL renders
```

## Where the active-site box lives (for reproducing docking)

Centred on the **chain-A active-site Cα centroid** of the residues `[80, 87, 109, 135, 175, 176, 195, 196, 214, 215, 217, 218, 221, 225, 226, 258]`:

- Centroid coordinates (Å): **x = −0.137, y = +4.232, z = +15.159**
- Box size: **22 × 22 × 22 Å** (v5 canonical) or **18 × 18 × 18 Å** (Phase 7 multi-replica)
- Vina exhaustiveness: **32** (v5) or **96** (v5 WT holo multi-seed sweep)
- Vina num_modes: 20 or 32
- Seed: 42 (canonical) + sanity {7, 13, 99, 256, 1, 2025, 31337}

The literal Vina invocation is in `12_phase7/01_replicas/VINA_COMMAND.md`.
