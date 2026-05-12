# Pipeline Roadmap

This is the roadmap the agent followed to take a single small/medium enzyme from raw UniProt accession all the way to a docking-and-mutagenesis report. It's broken into 9 stages plus environment bootstrap. Each stage produces files in its own numbered subfolder so a downstream stage can pick up where the previous one left off.

## Phase 0 — Environment bootstrap
- macOS Darwin / arm64, Python 3.11.9 (pyenv).
- Native binaries via Homebrew: `mafft`, `open-babel`, `pymol` 3.1.0, `brewsci/bio/autodock-vina` 1.2.7, plus their deps (`boost@1.85`, `glew`, `libxml2`).
- Python venv at `.venv/`: biopython, rdkit, prody, mdanalysis, meeko, gemmi, pyfamsa, matplotlib, seaborn, pandas, requests, jinja2, weasyprint, python-docx.
- The `vina` Python package failed to build; the project uses the Vina CLI via `subprocess` instead.
- The `pymol-open-source` pip wheel was installed but is broken (hardcoded paths from another user's mamba env); the project uses the Homebrew PyMOL binary headlessly via `subprocess` instead.

Full version manifest in [`00_setup/installed_libraries.md`](00_setup/installed_libraries.md).

## Phase 1 — Target choice
- Enzyme: **human Thymidylate Synthase (TYMS)**, UniProt **P04818**, 313 aa.
- Indication relevance: TYMS is the molecular target of **5-fluorouracil (5-FU)**, the backbone of colorectal-cancer chemotherapy.
- Reference structure: **PDB 1HVY** (1.9 Å, dUMP + raltitrexed). Natural substrate dUMP is the docking ligand.
- Rationale logged to [`pipeline.log`](pipeline.log).

## Phase 2 — Stages 1–9
| Stage | Folder | Inputs | Outputs |
| --- | --- | --- | --- |
| 1. MSA & per-residue Jensen-Shannon conservation | `01_msa/` | UniProt FASTA × 8 orthologs | `aligned.fa`, `conservation_scores.csv`, `conservation_plot.png` |
| 2. Active-site annotation (UniProt features + PDBe binding-site graph API) | `02_active_site/` | UniProt JSON, PDBe REST | `active_site_residues.csv`, `selected_residues.csv`, `overlap_figure.png` |
| 3. Structure prep (chain isolation, ligand split, protonation) | `03_structure/` | RCSB PDB | `protein_h.pdb`, `ligand.pdb`, `cofactor.pdb` |
| 4. Headless PyMOL ray-traced renders | `04_pymol/` | `protein_h.pdb`, `ligand.pdb`, conservation scores | 4× 1600×1200 PNGs (overview, closeup, conservation cartoon, cavity) |
| 5. Ligand prep (Gasteiger charges, PDBQT) | `05_ligand/` | `ligand.pdb` | `ligand.pdbqt` |
| 6. Wild-type Vina docking (exhaustiveness 16, seed 42, box 22³ Å) | `06_docking_wt/` | receptor + ligand PDBQT | `wt_poses.pdbqt`, `wt_result.json`, `wt_topdock.png` |
| 7. Mutagenesis panel + redocking (21 mutants) | `07_mut_docking/` | WT box + receptor preparation pipeline | per-mutant PDBQT + dock + screenshot, master `results_full.csv` |
| 8. Aggregate analysis & figures | `08_analysis/` | `results_full.csv`, conservation scores | `analysis.md` + 3 plots |
| 9. Report generation | `09_report/` | everything above | `report.html` (Jinja2, embedded base64 PNG), `report.pdf` (WeasyPrint), `report.docx` (python-docx) |

## Phase 3 — Multi-agent review
After the pipeline completed, four reviewers were spawned in parallel (read-only) and audited the output:
- [`reviews/01_validator.md`](reviews/01_validator.md) — file integrity & number reproduction.
- [`reviews/02_code_review.md`](reviews/02_code_review.md) — Python correctness & robustness.
- [`reviews/03_scientific_officer.md`](reviews/03_scientific_officer.md) — peer-review-grade scientific defensibility.
- [`reviews/04_structural_bioinformatician.md`](reviews/04_structural_bioinformatician.md) — deep technical methods audit.

The reviews surfaced two **critical** issues that invalidate the conservation half of the analysis (most "ortholog" UniProt accessions point to the wrong proteins; chain B was discarded even though catalytic residues come from the dimer partner) plus several lower-severity issues. These are documented openly in this repo and in the pipeline report; the docking / mutagenesis half (rigid receptor, single chain, Gasteiger charges) remains internally self-consistent.

## Phase 4 — Final deliverables
- `09_report/report.pdf` (2.5 MB) — primary report.
- `09_report/report.html` (3.5 MB) — self-contained, embedded PNGs.
- `09_report/report.docx` — Word version with embedded images.
- `workflow_diagram.png` — pipeline overview (this file embedded in the README below).
- All raw data, scripts, logs, and per-stage artefacts under their numbered folders.
