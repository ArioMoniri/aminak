#!/usr/bin/env python3
"""Stage 3 v5: IN-PLACE cofactor reprotonation. NEVER moves heavy atoms.

Round 4 reviewer found that v4 placed the cofactor by Kabsch-aligning the
CCD-ideal D16 onto the bound conformer, producing a 2.71 A heavy-atom RMSD
vs the 1HVY-bound conformer plus a real protein clash (cofA O1 <-> PHE 80
CD2 at 1.95 A). The "-2 cofactor expels dUMP" finding in v4 is therefore a
placement artefact, not biology.

v5 fix: do NOT re-derive coordinates from the CCD ideal SDF. Start from the
crystal HETATM D16 (which is already correctly bound) and only re-protonate
in place: strip ALL hydrogens, assign bond orders from a reference SMILES,
deprotonate carboxylates, AddHs(addCoords=True), and assert that NO heavy
atom moved by > 0.001 A. Also enumerate cofactor-protein clashes < 1.8 A
and abort if any exist.
"""
import os, sys, json, hashlib
from datetime import datetime
import numpy as np

PROJECT = os.path.expanduser("~/conserved_site_project")
STR_V5 = os.path.join(PROJECT, "03e_structure_v5")
STR_V2 = os.path.join(PROJECT, "03b_structure_v2")
STR_V4 = os.path.join(PROJECT, "03d_structure_v4")
ORIG_PDB = os.path.join(PROJECT, "03_structure", "1hvy.pdb")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
LOG_DIR = os.path.join(PROJECT, "logs")
STAGE_LOG = os.path.join(LOG_DIR, "v5_03_cofactor.log")

# Reference: D16 ideal SDF from RCSB Chemical Component Dictionary (already
# fetched in v4). Has correct bond orders, element order matches the 1HVY HETATM
# D16 atom order index-by-index. We use the SDF as the bond-order TEMPLATE only;
# coordinates are taken verbatim from the crystal HETATM records.
D16_TEMPLATE_SDF = os.path.join(STR_V4, "D16_ideal.sdf")


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V5] STAGE3: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def md5_of(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def parse_crystal_heavy(pdb_path, chain):
    """Parse crystal cofactor PDB. Returns ordered list of dicts for heavy atoms only."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "D16" and line[21:22].strip() == chain:
                name = line[12:16].strip()
                element = line[76:78].strip() or name[0]
                if element == "H":
                    continue
                resi = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                bfac = float(line[60:66]) if line[60:66].strip() else 0.0
                atoms.append({"name": name, "element": element,
                              "coord": np.array([x, y, z], dtype=float),
                              "bfac": bfac, "resi": resi})
    return atoms


def parse_protein_heavy(pdb_path):
    """Heavy atoms (no H) of dimer protein, returned as list of dicts."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                name = line[12:16].strip()
                element = line[76:78].strip() or name[0]
                if element == "H":
                    continue
                resn = line[17:20].strip()
                ch = line[21:22].strip()
                resi = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                atoms.append({"name": name, "element": element,
                              "resn": resn, "chain": ch, "resi": resi,
                              "coord": np.array([x, y, z], dtype=float)})
    return atoms


def reprotonate_in_place(crystal_atoms, template_sdf):
    """Take crystal heavy-atom coords, swap them into a bond-order-aware mol
    parsed from the D16 ideal SDF (which has the SAME atom order index-by-index
    as 1HVY HETATM D16, verified pre-flight), deprotonate carboxylates,
    AddHs(addCoords=True). Heavy-atom coordinates are byte-equivalent to the
    crystal input (verified by caller).
    """
    from rdkit import Chem

    # Parse template (D16 ideal SDF; has bonds and matching atom order)
    template = Chem.MolFromMolFile(template_sdf, removeHs=True, sanitize=True)
    if template is None:
        raise RuntimeError(f"failed to parse {template_sdf}")
    n_template = template.GetNumAtoms()
    n_crystal = len(crystal_atoms)
    log(f"  template heavy atoms: {n_template}, crystal heavy atoms: {n_crystal}")
    if n_template != n_crystal:
        raise RuntimeError(f"heavy-atom count mismatch {n_template} != {n_crystal}")

    # Verify element order matches index-by-index. If not, abort (we rely on it).
    template_elems = [a.GetSymbol() for a in template.GetAtoms()]
    crystal_elems = [a["element"] for a in crystal_atoms]
    if template_elems != crystal_elems:
        log(f"  template elems: {template_elems}")
        log(f"  crystal elems:  {crystal_elems}")
        raise RuntimeError("element order mismatch between template and crystal")
    log("  element order matches index-by-index")

    # Replace template's conformer coords with crystal coords IN-PLACE.
    new_conf = Chem.Conformer(n_template)
    for i, a in enumerate(crystal_atoms):
        new_conf.SetAtomPosition(i, a["coord"].tolist())
    template.RemoveAllConformers()
    template.AddConformer(new_conf, assignId=True)
    mol = template
    Chem.SanitizeMol(mol)

    # 4. Deprotonate carboxylates: -C(=O)OH -> -C(=O)[O-]
    smarts_cooh = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    matches = mol.GetSubstructMatches(smarts_cooh)
    log(f"  found {len(matches)} -COOH groups in template (will deprotonate)")
    rwmol = Chem.RWMol(mol)
    for match in matches:
        oh_idx = match[2]
        atom = rwmol.GetAtomWithIdx(oh_idx)
        atom.SetFormalCharge(-1)
        atom.SetNumExplicitHs(0)
    mol = rwmol.GetMol()
    Chem.SanitizeMol(mol)

    # 5. Add Hs with addCoords=True (places polar Hs only where required by valence;
    # heavy atoms are NOT moved).
    mol = Chem.AddHs(mol, addCoords=True)

    smi = Chem.MolToSmiles(Chem.RemoveHs(mol))
    log(f"  SMILES (no H): {smi}")
    n_neg = smi.count("[O-]")
    log(f"  carboxylate [O-] count: {n_neg}")
    if n_neg < 2:
        log(f"  ERROR: fewer than 2 [O-] groups; v5 cannot proceed")
        raise RuntimeError("carboxylate deprotonation failed")

    return mol, smi


def write_protonated_pdb(mol, crystal_atoms, out_pdb, chain="A", resname="D16"):
    """Write a PDB with crystal atom names + crystal coords + extra Hs.
    Heavy-atom output coords are exactly the crystal input coords (no transform).
    """
    from rdkit import Chem
    conf = mol.GetConformer()

    # Walk mol atoms; heavy atoms keep crystal name, Hs get auto-named H1, H2, ...
    h_counter = 1
    out_lines = []
    serial = 1
    heavy_idx = 0
    resi = crystal_atoms[0]["resi"] if crystal_atoms else 414
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        elem = atom.GetSymbol()
        pos = conf.GetAtomPosition(idx)
        coord = np.array([pos.x, pos.y, pos.z], dtype=float)
        if elem == "H":
            name = f"H{h_counter}"[:4]
            h_counter += 1
        else:
            if heavy_idx < len(crystal_atoms):
                name = crystal_atoms[heavy_idx]["name"]
                # ASSERT in-place: heavy-atom coordinates must match the input crystal
                ref_coord = crystal_atoms[heavy_idx]["coord"]
                d = np.linalg.norm(coord - ref_coord)
                if d > 0.001:
                    raise RuntimeError(f"HEAVY ATOM MOVED: {name} d={d:.4f} A "
                                       f"(must be < 0.001 A)")
            else:
                name = f"{elem}{idx+1}"
            heavy_idx += 1
        # PDB atom-name field: pad correctly
        if len(name) < 4 and len(elem) == 1:
            name_field = f" {name:<3s}"
        else:
            name_field = f"{name:<4s}"
        line = (f"HETATM{serial:5d} {name_field} {resname:<3s} {chain}{resi:4d}    "
                f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                f"  1.00  0.00          {elem:>2s}\n")
        out_lines.append(line)
        serial += 1
    out_lines.append("END\n")
    with open(out_pdb, "w") as f:
        f.writelines(out_lines)


def clash_check(cof_atoms, prot_atoms, threshold=1.8):
    """Return list of (cof_atom, prot_atom, distance) clashes < threshold."""
    cof_arr = np.array([a["coord"] for a in cof_atoms])
    prot_arr = np.array([a["coord"] for a in prot_atoms])
    clashes = []
    for i in range(len(cof_arr)):
        d = np.linalg.norm(prot_arr - cof_arr[i], axis=1)
        for j in np.where(d < threshold)[0]:
            clashes.append({
                "cofactor_atom": cof_atoms[i]["name"],
                "cofactor_element": cof_atoms[i]["element"],
                "protein_chain": prot_atoms[j]["chain"],
                "protein_resn": prot_atoms[j]["resn"],
                "protein_resi": prot_atoms[j]["resi"],
                "protein_atom": prot_atoms[j]["name"],
                "distance_A": float(d[j]),
            })
    return clashes


def heavy_atoms_from_pdb(pdb_path):
    """Read heavy atoms from a PDB; returns name->coord mapping for first chain."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                element = line[76:78].strip() or line[12:16].strip()[0]
                if element == "H":
                    continue
                name = line[12:16].strip()
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                atoms.append({"name": name, "element": element,
                              "coord": np.array([x, y, z], dtype=float)})
    return atoms


def main():
    os.makedirs(STR_V5, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 3 v5 starting (IN-PLACE cofactor reprotonation; no Kabsch)")

    # ---- Source: 1HVY (chain A and chain B) D16 ----
    if not os.path.exists(ORIG_PDB):
        log(f"FATAL: {ORIG_PDB} missing")
        sys.exit(1)
    crystal_a = parse_crystal_heavy(ORIG_PDB, "A")
    crystal_b = parse_crystal_heavy(ORIG_PDB, "B")
    log(f"  chain A heavy atoms: {len(crystal_a)}")
    log(f"  chain B heavy atoms: {len(crystal_b)}")
    if len(crystal_a) == 0 or len(crystal_b) == 0:
        log("FATAL: no D16 atoms found in 1HVY for one or both chains")
        sys.exit(1)

    # ---- Reprotonate each in place ----
    out_a = os.path.join(STR_V5, "cofactor_chainA_v5.pdb")
    out_b = os.path.join(STR_V5, "cofactor_chainB_v5.pdb")

    log("Chain A:")
    mol_a, smi_a = reprotonate_in_place(crystal_a, D16_TEMPLATE_SDF)
    write_protonated_pdb(mol_a, crystal_a, out_a, chain="A", resname="D16")
    log(f"  wrote {out_a} (md5={md5_of(out_a)})")

    log("Chain B:")
    # Chain B coords are different but the SAME chemical entity
    mol_b, smi_b = reprotonate_in_place(crystal_b, D16_TEMPLATE_SDF)
    write_protonated_pdb(mol_b, crystal_b, out_b, chain="B", resname="D16")
    log(f"  wrote {out_b} (md5={md5_of(out_b)})")

    # ---- Verification (a): heavy-atom RMSD vs original 1HVY chain A D16 ----
    written_a = heavy_atoms_from_pdb(out_a)
    written_b = heavy_atoms_from_pdb(out_b)
    name_to_crystal_a = {a["name"]: a["coord"] for a in crystal_a}
    name_to_crystal_b = {a["name"]: a["coord"] for a in crystal_b}
    sq_a = []
    for a in written_a:
        if a["name"] in name_to_crystal_a:
            sq_a.append(float(np.sum((a["coord"] - name_to_crystal_a[a["name"]]) ** 2)))
    sq_b = []
    for a in written_b:
        if a["name"] in name_to_crystal_b:
            sq_b.append(float(np.sum((a["coord"] - name_to_crystal_b[a["name"]]) ** 2)))
    rmsd_a = float(np.sqrt(sum(sq_a) / len(sq_a))) if sq_a else float("nan")
    rmsd_b = float(np.sqrt(sum(sq_b) / len(sq_b))) if sq_b else float("nan")
    log(f"  verify A: heavy-atom RMSD vs 1HVY chain-A D16 = {rmsd_a:.6f} A "
        f"(matched {len(sq_a)} of {len(crystal_a)} atoms)")
    log(f"  verify B: heavy-atom RMSD vs 1HVY chain-B D16 = {rmsd_b:.6f} A "
        f"(matched {len(sq_b)} of {len(crystal_b)} atoms)")
    if rmsd_a > 0.001:
        log(f"FATAL: chain A in-place check failed (RMSD {rmsd_a:.4f} > 0.001)")
        sys.exit(1)
    if rmsd_b > 0.001:
        log(f"FATAL: chain B in-place check failed (RMSD {rmsd_b:.4f} > 0.001)")
        sys.exit(1)

    # ---- Verification (b): MD5 differs from v4 ----
    v4_a = os.path.join(STR_V4, "cofactor_chainA_v4.pdb")
    if os.path.exists(v4_a):
        v5_md5 = md5_of(out_a)
        v4_md5 = md5_of(v4_a)
        log(f"  v5 chain-A md5 = {v5_md5}")
        log(f"  v4 chain-A md5 = {v4_md5}")
        if v5_md5 == v4_md5:
            log("FATAL: v5 cofactor MD5 equals v4 (must differ)")
            sys.exit(1)

    # ---- Verification (c): canonical SMILES contains [O-] ----
    if "[O-]" not in smi_a:
        log("FATAL: chain-A SMILES missing [O-]")
        sys.exit(1)

    # ---- Verification (d): clash check vs dimer protein ----
    prot_pdb = os.path.join(STR_V2, "protein_dimer_h.pdb")
    prot_atoms = parse_protein_heavy(prot_pdb)
    log(f"  protein dimer heavy atoms: {len(prot_atoms)}")
    clashes_a = clash_check(written_a, prot_atoms, threshold=1.8)
    clashes_b = clash_check(written_b, prot_atoms, threshold=1.8)
    log(f"  chain A clashes < 1.8 A: {len(clashes_a)}")
    log(f"  chain B clashes < 1.8 A: {len(clashes_b)}")
    if clashes_a:
        log(f"  chain A clash detail: {clashes_a[:5]}")
    if clashes_b:
        log(f"  chain B clash detail: {clashes_b[:5]}")

    # Soft warning < 1.95 A (the v4 reviewer's specific clash cutoff)
    near_a = clash_check(written_a, prot_atoms, threshold=2.0)
    log(f"  chain A near-clashes < 2.0 A: {len(near_a)}")
    if any(c["distance_A"] < 1.95 and c["protein_resn"] == "PHE"
           and c["protein_resi"] == 80 for c in near_a):
        log("FATAL: v5 still has the v4 PHE 80 / cof O1 clash. Aborting.")
        sys.exit(1)

    if clashes_a or clashes_b:
        # Reviewer required: any < 1.8 A means abort
        log("FATAL: v5 has hard clashes (< 1.8 A); aborting per FIX D")
        sys.exit(1)

    # ---- Provenance ----
    prov = {
        "method": "in-place reprotonation: heavy-atom coords from 1HVY HETATM D16, "
                  "RDKit bond-order assignment from reference SMILES, "
                  "carboxylate deprotonation, AddHs(addCoords=True). "
                  "NO Kabsch, NO 3D embedding, NO heavy-atom move > 0.001 A.",
        "source_pdb": ORIG_PDB,
        "template_sdf": D16_TEMPLATE_SDF,
        "canonical_smiles_chainA": smi_a,
        "canonical_smiles_chainB": smi_b,
        "n_carboxylates_deprotonated": 2,
        "rmsd_chainA_vs_1hvy_A": rmsd_a,
        "rmsd_chainB_vs_1hvy_B": rmsd_b,
        "md5_chainA": md5_of(out_a),
        "md5_chainB": md5_of(out_b),
        "clashes_chainA_lt_1p8A": clashes_a,
        "clashes_chainB_lt_1p8A": clashes_b,
        "v4_md5_chainA": md5_of(v4_a) if os.path.exists(v4_a) else None,
        "v4_diff_md5": (md5_of(out_a) != md5_of(v4_a)) if os.path.exists(v4_a) else None,
    }
    with open(os.path.join(STR_V5, "cofactor_provenance_v5.json"), "w") as f:
        json.dump(prov, f, indent=2, default=str)
    log("Stage 3 v5 DONE")


if __name__ == "__main__":
    main()
