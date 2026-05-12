#!/usr/bin/env python3
"""
Audit fix: every mutant apo complex PDB has 0 ligands.
Root cause: the v3 apo top-pose files use residue name UNL (Vina default).
When the v5 complex-builder concatenated protein + top_pose for the
apo condition, it relied on residue-name=UMP filter, so nothing went into
the complex as a HETATM.

Fix:
1. Take the WT crystal dUMP atom-name map (20 heavy atoms, names like
   N1, C2, O2, ..., O5', P, OP1, OP2, OP3) from 05b_ligand_v2/dump.pdb.
2. For every mutant apo top_pose.pdb (which has UNL residue) — rewrite each
   ATOM/HETATM line by:
     - residue name UMP, chain X, resi 414
     - atom names mapped from the WT reference by index (the PDBQT atom
       order is preserved across mutants for the same ligand).
3. Build a new <mut>_apo_complex.pdb by concatenating the corresponding
   mutant protein (chains A+B) + the rewritten top pose.

Run from project root after source 00_setup/env.sh.
"""
import os, pathlib, re

ROOT = pathlib.Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
VIEWER_DIR = ROOT / "07e_mut_docking_v5" / "viewer_files"
WT_REF = ROOT / "05b_ligand_v2" / "dump.pdb"
WT_HOLO_TOP = ROOT / "06e_docking_wt_v5" / "wt_apo_top_pose.pdb"

# Build the canonical atom-name list from the WT reference
ref_atom_names = []
for line in WT_REF.read_text().splitlines():
    if line.startswith(("ATOM","HETATM")):
        name = line[12:16].strip()
        elem = line[76:78].strip()
        if elem == "H" or name.startswith("H"): continue
        ref_atom_names.append(name)
print(f"reference UMP heavy-atom names ({len(ref_atom_names)}): {ref_atom_names}")

def rewrite_top_pose(text: str) -> str:
    """Take a top_pose.pdb with UNL residue, rewrite as proper UMP."""
    out = []
    i = 0
    for line in text.splitlines():
        if not line.startswith(("HETATM","ATOM")):
            # Drop PDBQT REMARK / MODEL / ENDMDL / ROOT / BRANCH lines
            continue
        # Skip hydrogens
        elem_field = line[76:78].strip() if len(line) >= 78 else ""
        name_field = line[12:16].strip()
        if elem_field == "H" or name_field.startswith("H"):
            continue
        if i >= len(ref_atom_names):
            break
        new_name = ref_atom_names[i]
        # Standard PDB atom-name padding: 4-char field, right-justified for elements
        if len(new_name) <= 3:
            atom_name = f" {new_name:<3s}"
        else:
            atom_name = f"{new_name:<4s}"
        # Rebuild the line with proper UMP / chain X / resi 414
        new_line = (
            "HETATM"                              # 1-6  record
            + f"{i+1:>5d}"                        # 7-11 serial
            + " "                                  # 12   blank
            + atom_name                            # 13-16 atom name
            + " "                                  # 17   altloc
            + "UMP"                                # 18-20 resname
            + " "                                  # 21   blank
            + "X"                                  # 22   chain
            + f"{414:>4d}"                         # 23-26 resi
            + " "                                  # 27   iCode
            + "   "                                # 28-30 blank
            + line[30:54]                          # 31-54 coords
            + "  1.00"                             # 55-60 occ
            + "  0.00"                             # 61-66 bf
            + "          "                         # 67-76 blank
            + f"{new_name[0]:>2s}"                 # 77-78 element
        )
        out.append(new_line)
        i += 1
    return "\n".join(out) + "\n"

def fix_one_mutant_apo(mut_dir: pathlib.Path):
    """Given a viewer-files-style basename like 'C195A_apo', rebuild
    the top_pose + complex."""
    pass

# Find every <mut>_apo_top_pose.pdb in viewer_files
top_poses = sorted(VIEWER_DIR.glob("*_apo_top_pose.pdb"))
print(f"\nfound {len(top_poses)} apo top_pose files")

# Determine the mutant directory layout to find the protein PDB
mut_root = ROOT / "07e_mut_docking_v5"

fixed = 0
for tp in top_poses:
    mut = tp.name.replace("_apo_top_pose.pdb","")
    # Find the protonated mutant protein. Skip if it doesn't exist.
    cand = [
        mut_root / mut / f"{mut}_mut_h.pdb",
        mut_root / mut / f"{mut}_mut.pdb",
        # v3 apo source (if reused)
        ROOT / "07c_mut_docking_v3" / mut / f"{mut}_mut_h.pdb",
        ROOT / "07c_mut_docking_v3" / mut / f"{mut}_mut.pdb",
    ]
    prot = next((c for c in cand if c.exists()), None)
    if prot is None:
        print(f"  [skip] {mut}: no protein PDB found in any of {[str(c.relative_to(ROOT)) for c in cand]}")
        continue
    # Rewrite top pose
    new_pose = rewrite_top_pose(tp.read_text())
    if not new_pose.strip():
        print(f"  [warn] {mut}: rewritten pose empty")
        continue
    (tp.parent / tp.name).write_text(new_pose)   # overwrite
    # Build complex: protein lines + TER + dUMP + END
    prot_text = prot.read_text()
    prot_lines = [l for l in prot_text.splitlines()
                  if l.startswith(("ATOM","TER","HETATM")) and not l.startswith("END")]
    complex_path = tp.parent / f"{mut}_apo_complex.pdb"
    out = "\n".join(prot_lines) + "\nTER\n" + new_pose + "TER\nEND\n"
    complex_path.write_text(out)
    fixed += 1
    print(f"  [ok]   {mut}: top_pose rewritten (UMP X 414), complex rebuilt ({complex_path.stat().st_size//1024} KB)")

print(f"\nfixed {fixed}/{len(top_poses)} mutant apo complexes")
