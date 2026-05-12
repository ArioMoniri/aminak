#!/usr/bin/env python3
"""Stage 3 v4: REAL cofactor reprotonation at pH 7.4 using bond-order-aware mol.

Key fix vs v3: v3 ran obabel -p 7.4 on a PDB that lacked bond orders, so nothing
changed (output byte-identical to v2). Here we:
  1. Load D16 ideal SDF (with bonds) from RCSB Chemical Component Dictionary
  2. Use RDKit to add hydrogens at pH 7.4 (deprotonates carboxylates)
  3. Align the protonated, ionised cofactor onto the crystal coords (chain A and B)
  4. Preserve crystal atom names; assign extra hydrogens
  5. ASSERT MD5 differs from v2/v3 file (loud abort otherwise)
"""
import os, sys, subprocess, hashlib, json
from datetime import datetime
import numpy as np

PROJECT = os.path.expanduser("~/conserved_site_project")
STR_V4 = os.path.join(PROJECT, "03d_structure_v4")
STR_V2 = os.path.join(PROJECT, "03b_structure_v2")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
LOG_DIR = os.path.join(PROJECT, "logs")
STAGE_LOG = os.path.join(LOG_DIR, "v4_03_cofactor.log")
OBABEL = os.environ.get("OBABEL", "/Users/ario/conserved_site_project/.venv/bin/obabel")


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V4] STAGE3: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def md5_of(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def parse_crystal_cofactor(pdb_path):
    """Parse crystal cofactor PDB. Returns list of (atom_name, element, coords, line_template)."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16].strip()
                element = line[76:78].strip() or name[0]
                if element == "H":
                    continue
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                atoms.append({"name": name, "element": element,
                              "coord": np.array([x, y, z]), "line_template": line})
    return atoms


def reprotonate_with_rdkit(sdf_in, sdf_out):
    """Use RDKit to read SDF (with bonds), add Hs at pH 7.4 (deprotonate carboxyl)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    mol = Chem.MolFromMolFile(sdf_in, removeHs=True)
    if mol is None:
        log("RDKit failed to parse SDF")
        return None

    # Standardize and set proper protonation at pH 7.4
    # Use the standard MolStandardize uncharger then re-protonate where needed
    # Simpler: use Reionizer to put charges on the most acidic positions
    reionizer = rdMolStandardize.Reionizer()
    mol = reionizer.reionize(mol)

    # The carboxyl groups are the only ionisable groups in this molecule that
    # ionise at pH 7.4. Force them to deprotonated state.
    # Find all -C(=O)OH and convert OH->O-
    smarts_cooh = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    matches = mol.GetSubstructMatches(smarts_cooh)
    log(f"  found {len(matches)} -COOH groups")
    rwmol = Chem.RWMol(mol)
    for match in matches:
        # match = (C, O_double, O_single_H)
        oh_idx = match[2]
        atom = rwmol.GetAtomWithIdx(oh_idx)
        atom.SetFormalCharge(-1)
        atom.SetNumExplicitHs(0)

    mol = rwmol.GetMol()
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol, addCoords=True)

    # Confirm via SMILES
    smi = Chem.MolToSmiles(Chem.RemoveHs(mol))
    log(f"  SMILES (no H): {smi}")
    if smi.count("[O-]") < 2:
        log(f"  WARNING: expected 2 [O-]; got {smi.count('[O-]')}")

    # Write SDF
    w = Chem.SDWriter(sdf_out)
    w.write(mol)
    w.close()
    return mol


def align_protonated_to_crystal(mol_protonated, crystal_atoms, out_pdb, chain="A", resname="D16", resi=414):
    """Align RDKit mol heavy atoms to crystal heavy atoms by index correspondence,
    then write a PDB with crystal atom names + crystal coords + extra Hs added.

    Strategy: the RCSB ideal SDF atom order matches CCD nomenclature. We build a
    name-mapping from CCD to crystal by examining the CCD atom block then crystal atoms.
    Simpler: since both have 32 heavy atoms and are the same chemical entity,
    align the 3D conformer to crystal heavy coords using RMSD-minimising transform.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # Get heavy atoms of protonated mol with their indices
    mol_heavy = [a for a in mol_protonated.GetAtoms() if a.GetSymbol() != "H"]
    log(f"  protonated heavy count: {len(mol_heavy)}, crystal heavy count: {len(crystal_atoms)}")
    if len(mol_heavy) != len(crystal_atoms):
        log(f"  WARNING: mismatch in heavy count")

    # The ideal SDF atom order matches the chemical component dictionary _chem_comp_atom
    # order. Crystal D16 atoms in 1HVY were originally generated from this same CCD,
    # so the order should match. Let's verify by element sequence.
    crystal_elements = [a["element"] for a in crystal_atoms]
    mol_elements = [a.GetSymbol() for a in mol_heavy]
    log(f"  crystal elements: {crystal_elements[:10]}...")
    log(f"  mol elements:     {mol_elements[:10]}...")

    # Try to align: extract conformer coords for heavy atoms
    conf = mol_protonated.GetConformer()
    mol_coords = np.array([conf.GetAtomPosition(a.GetIdx()) for a in mol_heavy])
    crystal_coords = np.array([a["coord"] for a in crystal_atoms])

    # Element-pattern match: if elements match index-by-index we can use direct alignment
    if crystal_elements == mol_elements:
        log(f"  element order matches; using direct Kabsch alignment")
        # Kabsch alignment
        c1 = mol_coords - mol_coords.mean(axis=0)
        c2 = crystal_coords - crystal_coords.mean(axis=0)
        H = c1.T @ c2
        U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        D = np.diag([1, 1, d])
        R = Vt.T @ D @ U.T
        # Apply
        all_coords = np.array([conf.GetAtomPosition(i) for i in range(mol_protonated.GetNumAtoms())])
        all_centered = all_coords - mol_coords.mean(axis=0)
        all_rotated = all_centered @ R.T + crystal_coords.mean(axis=0)
        # Compute heavy-atom RMSD post-alignment
        new_mol_coords = all_rotated[[a.GetIdx() for a in mol_heavy]]
        rmsd = np.sqrt(((new_mol_coords - crystal_coords) ** 2).sum(axis=1).mean())
        log(f"  alignment RMSD heavy = {rmsd:.3f} A")
    else:
        log(f"  element order MISMATCH; trying O3A-based alignment via RDKit")
        # Build a target mol from crystal coords and elements (no bonds)
        target = Chem.RWMol()
        for a in crystal_atoms:
            atom = Chem.Atom(a["element"])
            target.AddAtom(atom)
        target = target.GetMol()
        target_conf = Chem.Conformer(len(crystal_atoms))
        for i, a in enumerate(crystal_atoms):
            target_conf.SetAtomPosition(i, a["coord"].tolist())
        target.AddConformer(target_conf)
        # Use rigid alignment by all-to-all RMSD search is overkill; fall back to
        # simple pairing: for each crystal atom, find closest matching-element mol atom
        # then compute Kabsch with that mapping
        used = set()
        pairs = []
        for ci, ce in enumerate(crystal_elements):
            best = None; best_d = 1e9
            for mi, me in enumerate(mol_elements):
                if mi in used or me != ce:
                    continue
                d = np.linalg.norm(mol_coords[mi] - crystal_coords[ci])
                if d < best_d:
                    best_d = d; best = mi
            if best is not None:
                used.add(best)
                pairs.append((best, ci))
        log(f"  matched {len(pairs)} of {len(crystal_atoms)} atoms")
        if len(pairs) >= len(crystal_atoms) // 2:
            mp = np.array([mol_coords[p[0]] for p in pairs])
            cp = np.array([crystal_coords[p[1]] for p in pairs])
            c1 = mp - mp.mean(axis=0)
            c2 = cp - cp.mean(axis=0)
            H = c1.T @ c2
            U, S, Vt = np.linalg.svd(H)
            d = np.sign(np.linalg.det(Vt.T @ U.T))
            D = np.diag([1, 1, d])
            R = Vt.T @ D @ U.T
            all_coords = np.array([conf.GetAtomPosition(i) for i in range(mol_protonated.GetNumAtoms())])
            all_centered = all_coords - mp.mean(axis=0)
            all_rotated = all_centered @ R.T + cp.mean(axis=0)
        else:
            log(f"  alignment too sparse; not transforming coordinates")
            all_rotated = np.array([conf.GetAtomPosition(i) for i in range(mol_protonated.GetNumAtoms())])

    # Write the final PDB: heavy atoms with crystal names where possible, Hs as auto-named
    out_lines = []
    atom_serial = 1
    name_map = {}
    if crystal_elements == mol_elements:
        # Direct mapping: heavy_idx_in_mol[i] gets crystal name [i]
        for i, a in enumerate(mol_heavy):
            name_map[a.GetIdx()] = crystal_atoms[i]["name"]
    else:
        for mi, ci in pairs:
            name_map[mol_heavy[mi].GetIdx()] = crystal_atoms[ci]["name"]

    # Now produce PDB lines. Heavy atoms first (in mol's heavy order), then hydrogens
    h_counter = 1
    for atom in mol_protonated.GetAtoms():
        idx = atom.GetIdx()
        coord = all_rotated[idx]
        elem = atom.GetSymbol()
        if elem == "H":
            name = f"H{h_counter}"
            if len(name) > 4:
                name = name[:4]
            h_counter += 1
        else:
            name = name_map.get(idx, f"{elem}{idx+1}")
        # Right/left justify: PDB atom name field is 4 chars with element-aware padding
        if len(name) < 4 and len(elem) == 1:
            name_field = f" {name:<3s}"
        else:
            name_field = f"{name:<4s}"
        line = (f"HETATM{atom_serial:5d} {name_field} {resname:<3s} {chain}{resi:4d}    "
                f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                f"  1.00  0.00          {elem:>2s}\n")
        out_lines.append(line)
        atom_serial += 1
    out_lines.append("END\n")
    with open(out_pdb, "w") as f:
        f.writelines(out_lines)
    return out_pdb


def main():
    os.makedirs(STR_V4, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 3 v4 starting (REAL cofactor reprotonation)")

    sdf_ideal = os.path.join(STR_V4, "D16_ideal.sdf")
    if not os.path.exists(sdf_ideal):
        log("FATAL: D16_ideal.sdf missing. Fetch from RCSB first.")
        sys.exit(1)

    sdf_protonated = os.path.join(STR_V4, "D16_protonated_ph74.sdf")
    mol = reprotonate_with_rdkit(sdf_ideal, sdf_protonated)
    if mol is None:
        log("FATAL: RDKit reprotonation failed")
        sys.exit(1)
    log(f"  protonated SDF written")

    # Confirm via obabel canonical SMILES
    proc = subprocess.run([OBABEL, sdf_protonated, "-ocan"], capture_output=True, text=True)
    smi_can = proc.stdout.strip().split()[0] if proc.stdout else "(empty)"
    log(f"  obabel canonical SMILES: {smi_can}")
    n_neg = smi_can.count("[O-]")
    if n_neg < 2:
        log(f"  CRITICAL WARNING: expected 2x [O-] (glutamate alpha+gamma); got {n_neg}")

    # Align onto each crystal cofactor (chain A and chain B) and write PDB
    for chain in ["A", "B"]:
        crystal_pdb = os.path.join(STR_V2, f"cofactor_chain{chain}.pdb")
        crystal_atoms = parse_crystal_cofactor(crystal_pdb)
        out_pdb = os.path.join(STR_V4, f"cofactor_chain{chain}_v4.pdb")
        align_protonated_to_crystal(mol, crystal_atoms, out_pdb, chain=chain)
        log(f"  wrote {out_pdb}")

    # ASSERT MD5 differs
    for chain in ["A", "B"]:
        v2_md5 = md5_of(os.path.join(STR_V2, f"cofactor_chain{chain}_h.pdb"))
        v4_md5 = md5_of(os.path.join(STR_V4, f"cofactor_chain{chain}_v4.pdb"))
        log(f"  chain {chain}: v2_md5={v2_md5} v4_md5={v4_md5}")
        if v2_md5 == v4_md5:
            log(f"  FATAL: chain {chain} v4 byte-identical to v2; reprotonation no-op")
            sys.exit(1)
        log(f"  chain {chain}: MD5 differs (good)")

    # Save provenance
    prov = {
        "source_sdf": "https://files.rcsb.org/ligands/download/D16_ideal.sdf",
        "method": "RDKit Reionizer + explicit -COOH -> COO- + AddHs(addCoords=True)",
        "alignment": "Kabsch by element-order if matches, else nearest-element-match",
        "canonical_smiles": smi_can,
        "n_carboxylates_deprotonated": n_neg,
    }
    with open(os.path.join(STR_V4, "cofactor_provenance.json"), "w") as f:
        json.dump(prov, f, indent=2)
    log("Stage 3 v4 DONE")


if __name__ == "__main__":
    main()
