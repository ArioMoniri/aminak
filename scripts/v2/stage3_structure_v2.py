#!/usr/bin/env python3
"""Stage 3 v2: Keep dimer (chains A & B), restore CME43 -> CYS, split UMP & D16 per chain."""
import os, sys, subprocess, json
from datetime import datetime
import requests
from Bio.PDB import PDBParser, PDBIO, Select

PROJECT = os.path.expanduser("~/conserved_site_project")
STR_DIR_V1 = os.path.join(PROJECT, "03_structure")
STR_DIR = os.path.join(PROJECT, "03b_structure_v2")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v2_03_structure.log")
OBABEL = os.environ.get("OBABEL", "/Users/ario/conserved_site_project/.venv/bin/obabel")

# Standard CYS heavy-atom names
CYS_HEAVY = {"N", "CA", "C", "O", "CB", "SG"}


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V2] STAGE3: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def main():
    os.makedirs(STR_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 3 v2 starting")

    # Reuse the v1 1HVY raw file
    raw = os.path.join(STR_DIR_V1, "1hvy.pdb")
    if not os.path.exists(raw):
        log("v1 1hvy.pdb missing — refetching")
        r = requests.get("https://files.rcsb.org/download/1HVY.pdb", timeout=60)
        r.raise_for_status()
        with open(raw, "w") as f:
            f.write(r.text)
    log(f"using raw 1hvy: {raw}")

    # Pre-process raw PDB lines: rewrite CME records to CYS records (drop hydroxyethyl atoms)
    cleaned_raw = os.path.join(STR_DIR, "1hvy_cme_to_cys.pdb")
    n_cme_lines = 0
    n_kept = 0
    n_dropped = 0
    with open(raw) as fin, open(cleaned_raw, "w") as fout:
        for line in fin:
            rec = line[:6]
            if rec in ("ATOM  ", "HETATM"):
                resname = line[17:20]
                atom_name = line[12:16].strip()
                chain_id = line[21]
                # CME on any chain: rename to CYS, drop atoms not in CYS_HEAVY
                if resname == "CME":
                    n_cme_lines += 1
                    if atom_name in CYS_HEAVY:
                        # Convert HETATM -> ATOM, rename CME -> CYS
                        new_line = "ATOM  " + line[6:17] + "CYS" + line[20:]
                        fout.write(new_line)
                        n_kept += 1
                    else:
                        n_dropped += 1
                    continue
            fout.write(line)
    log(f"CME->CYS conversion: {n_cme_lines} CME lines, {n_kept} kept (heavy), {n_dropped} hydroxyethyl atoms dropped")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("1HVY", cleaned_raw)

    # Selectors
    class ProteinDimer(Select):
        def accept_chain(self, chain): return chain.id in ("A", "B")
        def accept_residue(self, residue): return residue.id[0] == " "
        def accept_model(self, model): return model.id == 0

    class LigandUMP_chainA(Select):
        def accept_chain(self, chain): return chain.id == "A"
        def accept_residue(self, residue): return residue.get_resname() == "UMP"
        def accept_model(self, model): return model.id == 0

    class CofactorD16_chainA(Select):
        def accept_chain(self, chain): return chain.id == "A"
        def accept_residue(self, residue): return residue.get_resname() == "D16"
        def accept_model(self, model): return model.id == 0

    class CofactorD16_chainB(Select):
        def accept_chain(self, chain): return chain.id == "B"
        def accept_residue(self, residue): return residue.get_resname() == "D16"
        def accept_model(self, model): return model.id == 0

    io = PDBIO()
    io.set_structure(structure)

    prot_dimer = os.path.join(STR_DIR, "protein_dimer.pdb")
    lig_pdb = os.path.join(STR_DIR, "ligand.pdb")
    cof_a_pdb = os.path.join(STR_DIR, "cofactor_chainA.pdb")
    cof_b_pdb = os.path.join(STR_DIR, "cofactor_chainB.pdb")

    io.save(prot_dimer, ProteinDimer()); log(f"wrote {prot_dimer}")
    io.save(lig_pdb, LigandUMP_chainA()); log(f"wrote {lig_pdb}")
    io.save(cof_a_pdb, CofactorD16_chainA()); log(f"wrote {cof_a_pdb}")
    io.save(cof_b_pdb, CofactorD16_chainB()); log(f"wrote {cof_b_pdb}")

    # Add hydrogens to dimer
    prot_h = os.path.join(STR_DIR, "protein_dimer_h.pdb")
    proc = subprocess.run([OBABEL, prot_dimer, "-O", prot_h, "-h"],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        log(f"obabel protein_h failed: {proc.stderr[:500]}")
        sys.exit(1)
    log(f"protein_dimer_h.pdb size {os.path.getsize(prot_h)}")

    # Add hydrogens to ligand
    lig_h = os.path.join(STR_DIR, "ligand_h.pdb")
    proc = subprocess.run([OBABEL, lig_pdb, "-O", lig_h, "-h"],
                          capture_output=True, text=True)
    log(f"obabel ligand_h rc={proc.returncode}, size {os.path.getsize(lig_h) if os.path.exists(lig_h) else 0}")

    # Add hydrogens to cofactors
    for src, dst in [(cof_a_pdb, os.path.join(STR_DIR, "cofactor_chainA_h.pdb")),
                     (cof_b_pdb, os.path.join(STR_DIR, "cofactor_chainB_h.pdb"))]:
        proc = subprocess.run([OBABEL, src, "-O", dst, "-h"], capture_output=True, text=True)
        log(f"  {os.path.basename(src)} -> {os.path.basename(dst)} rc={proc.returncode}")

    # Verify both chains and target residues present
    s = parser.get_structure("ph", prot_h)
    chains = list(s[0].get_chains())
    log(f"chains in dimer_h: {[c.id for c in chains]}")
    for c in chains:
        n_res = sum(1 for r in c.get_residues() if r.id[0] == " ")
        log(f"  chain {c.id}: {n_res} residues")
        # Check Cys43 (was CME) is present
        for r in c.get_residues():
            if r.id[1] == 43 and r.get_resname() == "CYS":
                atoms = [a.get_name() for a in r if a.element != "H"]
                log(f"  chain {c.id} Cys43 heavy atoms: {sorted(atoms)}")
                break

    sel = json.load(open(os.path.join(PROJECT, "02b_active_site_v2/selected_meta.json")))["selected"]
    chainA_resids = {r.id[1] for r in s[0]["A"].get_residues() if r.id[0] == " "}
    chainB_resids = {r.id[1] for r in s[0]["B"].get_residues() if r.id[0] == " "}
    missing_A = [p for p in sel if p not in chainA_resids]
    missing_B = [p for p in sel if p not in chainB_resids]
    log(f"selected residues missing in A: {missing_A}; in B: {missing_B}")

    log("Stage 3 v2 DONE")


if __name__ == "__main__":
    main()
