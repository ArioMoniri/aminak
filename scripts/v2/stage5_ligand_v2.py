#!/usr/bin/env python3
"""Stage 5 v2: Multi-format ligand prep (PDB, MOL2, SDF, PDBQT)."""
import os, sys, subprocess
from datetime import datetime

PROJECT = os.path.expanduser("~/conserved_site_project")
LIG_DIR = os.path.join(PROJECT, "05b_ligand_v2")
STR_DIR = os.path.join(PROJECT, "03b_structure_v2")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v2_05_ligand.log")
OBABEL = os.environ.get("OBABEL", "/Users/ario/conserved_site_project/.venv/bin/obabel")
VENV_PY = os.path.join(PROJECT, ".venv/bin/python")
MK_PREP_LIG = os.path.join(PROJECT, ".venv/bin/mk_prepare_ligand.py")


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V2] STAGE5: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def main():
    os.makedirs(LIG_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 5 v2 starting")

    src_pdb = os.path.join(STR_DIR, "ligand.pdb")          # raw UMP from crystal (no H)
    src_h_pdb = os.path.join(STR_DIR, "ligand_h.pdb")      # with H

    dst_pdb = os.path.join(LIG_DIR, "dump.pdb")
    dst_mol2 = os.path.join(LIG_DIR, "dump.mol2")
    dst_sdf = os.path.join(LIG_DIR, "dump.sdf")
    dst_pdbqt = os.path.join(LIG_DIR, "dump.pdbqt")

    # 1) Copy PDB (no-H crystal)
    subprocess.run(["cp", src_pdb, dst_pdb], check=True)
    log(f"wrote {dst_pdb}")

    # 2) MOL2 from H-added PDB
    proc = subprocess.run([OBABEL, src_h_pdb, "-O", dst_mol2], capture_output=True, text=True)
    log(f"obabel mol2 rc={proc.returncode} size={os.path.getsize(dst_mol2) if os.path.exists(dst_mol2) else 0}")

    # 3) SDF with --gen3d
    proc = subprocess.run([OBABEL, src_h_pdb, "-O", dst_sdf, "--gen3d"], capture_output=True, text=True)
    log(f"obabel sdf --gen3d rc={proc.returncode} size={os.path.getsize(dst_sdf) if os.path.exists(dst_sdf) else 0}")
    if not os.path.exists(dst_sdf) or os.path.getsize(dst_sdf) < 100:
        # fallback without --gen3d
        log("  sdf --gen3d failed, trying without --gen3d")
        proc = subprocess.run([OBABEL, src_h_pdb, "-O", dst_sdf], capture_output=True, text=True)
        log(f"  obabel sdf no-gen3d rc={proc.returncode} size={os.path.getsize(dst_sdf) if os.path.exists(dst_sdf) else 0}")

    # 4) PDBQT (Vina-ready, partial charges)
    proc = subprocess.run([OBABEL, src_h_pdb, "-O", dst_pdbqt,
                          "--partialcharge", "gasteiger", "-p", "7.4"],
                         capture_output=True, text=True)
    log(f"obabel pdbqt rc={proc.returncode} size={os.path.getsize(dst_pdbqt) if os.path.exists(dst_pdbqt) else 0}")
    if proc.returncode != 0 or os.path.getsize(dst_pdbqt) < 100:
        log("  obabel pdbqt failed, trying meeko fallback")
        if os.path.exists(MK_PREP_LIG):
            proc = subprocess.run([VENV_PY, MK_PREP_LIG, "-i", src_h_pdb, "-o", dst_pdbqt],
                                 capture_output=True, text=True)
            log(f"  meeko rc={proc.returncode}")

    # Validate (no UNK atom types)
    if os.path.exists(dst_pdbqt):
        with open(dst_pdbqt) as f:
            content = f.read()
        log(f"pdbqt has_UNK={' UNK ' in content}, n_atoms={content.count('ATOM ') + content.count('HETATM')}")

    log("Stage 5 v2 DONE")


if __name__ == "__main__":
    main()
