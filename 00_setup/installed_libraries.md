# Installed Libraries & Tools — Conserved Site Pipeline

Generated 2026-05-12. Project: ~/conserved_site_project/

## System

| Component | Version | Source |
| --- | --- | --- |
| macOS | Darwin 25.4.0 (arm64) | system |
| Python | 3.11.9 | pyenv |
| Homebrew | latest | system |

## Native binaries (Homebrew)

| Tool | Version | Formula |
| --- | --- | --- |
| AutoDock Vina | AutoDock Vina v1.2.7 | brewsci/bio/autodock-vina |
| PyMOL | 3.1.0_3 (headless render OK) | homebrew/core/pymol |
| MAFFT | v7.526 (2024/Apr/26) | homebrew/core/mafft |
| OpenBabel (CLI) | Open Babel 3.1.0 -- May  8 2020 -- 15:48:11 | homebrew/core/open-babel |
| Boost C++ libs | 1.85.0_3 (vina dep) | brewsci/bio/boost@1.85 |
| GLEW | 2.3.1 (pymol GL dep) | homebrew/core/glew |
| libxml2 | system + homebrew/core/libxml2 | homebrew/core/libxml2 |

## Python venv (.venv, Python 3.11.9)

Selected scientific packages (full freeze in `pip_freeze.txt`).

| Package | Version | Used for |
| --- | --- | --- |
| biopython | 1.87 | MSA parsing, PDB I/O, FASTA |
| rdkit | 2026.3.1 | ligand handling, charge fallback |
| prody | 2.6.1 | structure analysis (loaded as sanity check) |
| MDAnalysis | 2.10.0 | trajectory libs (kept for future MD extension) |
| meeko | 0.7.1 | AutoDock Vina prep helpers (mk_prepare_*) |
| gemmi | 0.7.5 | structural file I/O (meeko dep) |
| pyfamsa | 0.7.0 | MSA fallback (unused; MAFFT used) |
| matplotlib | 3.10.9 | figures |
| matplotlib-venn | 1.1.2 | Stage 2 overlap diagram |
| seaborn | 0.13.2 | heatmap, scatter |
| pandas | 3.0.3 | results tables |
| numpy | 2.4.4 | JSD + numerics |
| scipy | 1.17.1 | stats helpers |
| requests | 2.34.0 | UniProt / PDBe REST |
| jinja2 | 3.1.6 | report templating |
| weasyprint | 68.1 | HTML→PDF |
| python-docx | 1.2.0 | HTML→DOCX (final report) |
| openbabel-wheel | 3.1.1.23 | Python OpenBabel bindings (CLI used) |

## Tools attempted but NOT successfully installed (with workarounds)

| Attempted | Failure mode | Workaround used |
| --- | --- | --- |
| `vina` Python pkg (pip) | Build failed: missing Boost headers in pip wheel | Used Vina CLI v1.2.7 via subprocess |
| `pymol-open-source` (pip) | Wheel hardcoded paths from another user's mamba env (libGLEW.2.1, libxml2 from /Users/Martin/) | Used Homebrew PyMOL 3.1.0 binary headless via subprocess |
| Modeller / Rosetta | Not installed (heavy, license) | PyMOL Mutagenesis Wizard rotamer pick (acknowledged caveat) |
