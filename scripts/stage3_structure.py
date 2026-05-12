#!/usr/bin/env python3
"""Stage 3: Download 1HVY, clean to chain A, split ligand, add hydrogens."""
import os, sys, subprocess, json
from datetime import datetime
import requests
from Bio.PDB import PDBParser, PDBIO, Select

PROJECT = os.path.expanduser("~/conserved_site_project")
STR_DIR = os.path.join(PROJECT, "03_structure")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
OBABEL = os.environ.get("OBABEL", "/Users/ario/conserved_site_project/.venv/bin/obabel")

def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] STAGE3: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

def main():
    os.makedirs(STR_DIR, exist_ok=True)
    raw = os.path.join(STR_DIR, "1hvy.pdb")
    log("Downloading 1HVY")
    if not os.path.exists(raw) or os.path.getsize(raw) < 1000:
        try:
            r = requests.get("https://files.rcsb.org/download/1HVY.pdb", timeout=60)
            r.raise_for_status()
            with open(raw, "w") as f:
                f.write(r.text)
            log("Got from RCSB")
        except Exception as e:
            log(f"RCSB failed: {e}, trying PDBe")
            r = requests.get("https://www.ebi.ac.uk/pdbe/entry-files/download/pdb1hvy.ent", timeout=60)
            r.raise_for_status()
            with open(raw, "w") as f:
                f.write(r.text)
    log(f"raw 1hvy size: {os.path.getsize(raw)} bytes")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("1HVY", raw)

    # Identify chain A
    model = structure[0]
    chain_a = model["A"]

    # Splits
    class ProteinChainA(Select):
        def accept_chain(self, chain): return chain.id == "A"
        def accept_residue(self, residue):
            # standard amino acids only (hetflag == ' ')
            return residue.id[0] == " "
        def accept_model(self, model): return model.id == 0

    class LigandUMP(Select):
        def accept_chain(self, chain): return chain.id == "A"
        def accept_residue(self, residue):
            return residue.get_resname() == "UMP"
        def accept_model(self, model): return model.id == 0

    class CofactorD16(Select):
        def accept_chain(self, chain): return chain.id == "A"
        def accept_residue(self, residue):
            return residue.get_resname() == "D16"
        def accept_model(self, model): return model.id == 0

    io = PDBIO()
    io.set_structure(structure)

    prot_chainA = os.path.join(STR_DIR, "protein_chainA.pdb")
    lig_pdb = os.path.join(STR_DIR, "ligand.pdb")
    cof_pdb = os.path.join(STR_DIR, "cofactor.pdb")
    io.save(prot_chainA, ProteinChainA()); log(f"wrote {prot_chainA}")
    io.save(lig_pdb, LigandUMP()); log(f"wrote {lig_pdb}")
    io.save(cof_pdb, CofactorD16()); log(f"wrote {cof_pdb}")

    # Add hydrogens
    log("Adding hydrogens to protein with obabel")
    prot_h = os.path.join(STR_DIR, "protein_h.pdb")
    proc = subprocess.run([OBABEL, prot_chainA, "-O", prot_h, "-h"],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        log(f"obabel protein H failed: {proc.stderr[:500]}")
    else:
        log(f"protein_h.pdb size {os.path.getsize(prot_h)}")

    log("Adding hydrogens to ligand with obabel")
    lig_h = os.path.join(STR_DIR, "ligand_h.pdb")
    proc = subprocess.run([OBABEL, lig_pdb, "-O", lig_h, "-h"],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        log(f"obabel ligand H failed: {proc.stderr[:500]}")
    else:
        log(f"ligand_h.pdb size {os.path.getsize(lig_h)}")

    # Verify selected residues exist in chain A
    sel = json.load(open(os.path.join(PROJECT, "02_active_site/selected_meta.json")))["selected"]
    ph_struct = parser.get_structure("ph", prot_h)
    chainA_resids = {r.id[1] for r in ph_struct[0]["A"].get_residues() if r.id[0] == " "}
    missing = [p for p in sel if p not in chainA_resids]
    log(f"Selected residues check: {len(sel)} requested, {len(missing)} missing: {missing}")
    if missing:
        log("WARNING: some active site residues missing from cleaned structure")

    log("Stage 3 DONE")

if __name__ == "__main__":
    main()
