#!/usr/bin/env python3
"""Package the PPTX with all deep-link target files into a single downloadable ZIP.

User unzips the bundle → all relative-path hyperlinks in the PPTX resolve
because the linked files (PDB, HTML viewers, CSV, MD configs) live in the
right relative locations alongside the PPTX.
"""
from pathlib import Path
import shutil, zipfile

REPO = Path(__file__).resolve().parents[2]
PRES_DIR = REPO / "14_inhibitor_design" / "presentation"
BUNDLE_DIR = PRES_DIR / "_bundle"
BUNDLE_ZIP = REPO / "14_inhibitor_design" / "presentation" / "aminak_phase14_bundle.zip"

# Files/folders to include in the bundle
# Layout inside the ZIP mirrors the repo so the relative PPTX hyperlinks work
INCLUDE = [
    # PPTX itself
    "14_inhibitor_design/presentation/aminak_phase14_summary.pptx",
    # Receptor PDBQTs (linked from Phase 3+4 slide)
    "06f_receptor_fixed/protein_dimer_apo_fixed.pdbqt",
    "06f_receptor_fixed/protein_dimer_holo_fixed.pdbqt",
    "06f_receptor_fixed/cofactor_A.pdbqt",
    "06f_receptor_fixed/cofactor_B.pdbqt",
    # Phase 14 master CSV + per-strategy summaries
    "14_inhibitor_design/05_aggregate/master.csv",
    "14_inhibitor_design/01_active_site/results_summary.csv",
    "14_inhibitor_design/02_cofactor_site/results_summary.csv",
    "14_inhibitor_design/03_dimer_interface/results_summary.csv",
    "14_inhibitor_design/04_allosteric/results_summary.csv",
    # Roadmap + reviews
    "14_inhibitor_design/00_roadmap/ROADMAP.md",
    "14_inhibitor_design/00_roadmap/anchor_compounds_verified.json",
    # Cavity 18 evidence — viewers + downloads
    "14_inhibitor_design/04_allosteric/cavity18_evidence/viewers/cavity18_apo.html",
    "14_inhibitor_design/04_allosteric/cavity18_evidence/viewers/cavity18_indazole.html",
    "14_inhibitor_design/04_allosteric/cavity18_evidence/viewers/cavity18_ibuprofen.html",
    "14_inhibitor_design/04_allosteric/cavity18_evidence/downloads/cavity18_apo.pdb",
    "14_inhibitor_design/04_allosteric/cavity18_evidence/downloads/cavity18_apo_pocket.pdb",
    "14_inhibitor_design/04_allosteric/cavity18_evidence/downloads/cavity18_indazole_complex.pdb",
    "14_inhibitor_design/04_allosteric/cavity18_evidence/downloads/cavity18_ibuprofen_complex.pdb",
    "14_inhibitor_design/04_allosteric/cavity18_evidence/downloads/cavity18_residues.csv",
    "14_inhibitor_design/04_allosteric/cavity18_evidence/downloads/cavity18_mutations_per_taxon.json",
    "14_inhibitor_design/04_allosteric/cavity18_residues.txt",
    # Smina rescore artefacts
    "14_inhibitor_design/06_smina_rescore/custom_scoring_q.txt",
    "14_inhibitor_design/06_smina_rescore/custom_scoring_qamp.txt",
    "14_inhibitor_design/06_smina_rescore/rescore_results.csv",
    "14_inhibitor_design/06_smina_rescore/rescore_summary.csv",
    "14_inhibitor_design/06_smina_rescore/rescore_plot.png",
    # Advanced methods — MM-GBSA + HADDOCK3 configs
    "14_inhibitor_design/07_advanced_methods/mmgbsa/README.md",
    "14_inhibitor_design/07_advanced_methods/mmgbsa/tleap.in",
    "14_inhibitor_design/07_advanced_methods/mmgbsa/mmpbsa.in",
    "14_inhibitor_design/07_advanced_methods/haddock3/README.md",
    "14_inhibitor_design/07_advanced_methods/haddock3/haddock3_config.cfg",
    "14_inhibitor_design/07_advanced_methods/haddock3/active_residues_chainA.txt",
    "14_inhibitor_design/07_advanced_methods/ppi_per_residue_bsa.csv",
    "14_inhibitor_design/07_advanced_methods/ppi_summary.json",
    "14_inhibitor_design/07_advanced_methods/smina_full_panel.csv",
    "14_inhibitor_design/07_advanced_methods/modeller_vs_alphafold.csv",
    # README files for navigation
    "README.md",
    "14_inhibitor_design/README.md",
    "TECHNICAL_NOTES.md",
]

# Write a bundle README so the user knows what's inside
BUNDLE_README = """# Aminak Phase 14 — PPTX bundle

Unzip this folder anywhere. Open `14_inhibitor_design/presentation/aminak_phase14_summary.pptx`.
Click the **deep-linked text** on slides to open the underlying files (PDBs in PyMOL/ChimeraX,
HTMLs in your browser, CSVs in Excel/Numbers).

## Quick links
- `14_inhibitor_design/presentation/aminak_phase14_summary.pptx`  — the slide deck
- `README.md`  — full project README
- `14_inhibitor_design/README.md`  — Phase 14 educational README
- `14_inhibitor_design/00_roadmap/ROADMAP.md`  — roadmap + reviewer chain
- `14_inhibitor_design/05_aggregate/master.csv`  — 86-row cross-strategy result table

## What the deep links open
- The dimer-aware Phase 3+4 slide links to the prepared receptor PDBQT.
- The Cavity-18 pose slides link to the 3Dmol.js HTML viewers — open them in any modern browser.
- The Smina rescoring slide links to the custom scoring files and CSV tables.
- The HADDOCK3 / MM-GBSA roadmap slides link to ready-to-run configs.

## What's NOT in the bundle
- The raw docked PDBQTs (1500+ files, ~100 MB) — they remain in the GitHub repo under
  `14_inhibitor_design/*/docked/`.
- The Phase 7 / 8 mutant pose libraries — same; in the repo under `13_phase8/` and `12_phase7/`.

Browse the live version at https://github.com/ArioMoniri/aminak.
"""

def main():
    BUNDLE_ZIP.parent.mkdir(parents=True, exist_ok=True)
    if BUNDLE_DIR.exists(): shutil.rmtree(BUNDLE_DIR)
    BUNDLE_DIR.mkdir(parents=True)
    # Copy files into the bundle dir, preserving repo structure
    n_files = 0; total_size = 0
    for rel in INCLUDE:
        src = REPO / rel
        dst = BUNDLE_DIR / rel
        if not src.exists():
            print(f"  ! missing: {src}"); continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        n_files += 1; total_size += dst.stat().st_size
    # Write bundle README
    (BUNDLE_DIR / "README_BUNDLE.md").write_text(BUNDLE_README)
    n_files += 1
    print(f"  → {n_files} files staged in {BUNDLE_DIR} ({total_size/1e6:.1f} MB)")
    # Zip it
    if BUNDLE_ZIP.exists(): BUNDLE_ZIP.unlink()
    with zipfile.ZipFile(BUNDLE_ZIP, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in BUNDLE_DIR.rglob("*"):
            if f.is_file():
                zf.write(f, f.relative_to(BUNDLE_DIR.parent))
    print(f"  → {BUNDLE_ZIP}  ({BUNDLE_ZIP.stat().st_size/1e6:.1f} MB)")

if __name__ == "__main__":
    main()
