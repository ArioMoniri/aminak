#!/usr/bin/env python3
"""Custom flex-residue PDBQT generator for AutoDock Vina.

Why a custom tool?  Meeko's mk_prepare_receptor.py and its prody/RDKit
backends repeatedly fail on the OpenBabel-derived PDB files in this
project ("Updated 1 H positions but deleted 7", "Explicit valence ...
greater than permitted").  AutoDockTools' prepare_flexreceptor4.py is
not installed (no MGLTools binaries on Apple Silicon).

The Vina flex-residue PDBQT format is well-defined: per residue, a
torsion tree of ROOT (CA) -> BRANCH (chi1, CA-CB) -> ATOM CB ->
BRANCH (chi2, CB-CG) -> ATOM CG -> ... -> ENDBRANCH ... -> ENDROOT.
Backbone N, CA, C, O stay rigid (CA appears in both rigid and flex but
the rigid CA is what holds the side chain in place).

This module provides:
- Hardcoded chi-rotation chain templates for the 20 standard amino acids
  (covers the 14-residue active-site panel of dTMP synthase).
- split_clean_pdbqt(rec_pdbqt, flex_residues) -> (rigid_str, flex_str)
  which parses an OpenBabel + Gasteiger PDBQT (chains A/B intact) and
  emits the two PDBQT strings.

Each amino-acid template is a list of "branch arms".  Each arm is a
sequence of (parent_atom, child_atom) pairs that defines the torsion
chain from CA outward.  The first pair (CA -> CB) is always present
(except GLY/ALA).  Subsequent pairs add another rotatable bond.
"""
from __future__ import annotations

from pathlib import Path

# --- amino-acid side-chain topology --------------------------------------
# For each residue: an ordered list of *bonds* representing the rotatable
# side-chain spine from CA outward.  Atoms not in the spine but attached
# to a spine atom belong to the same BRANCH segment (no further rotation).
#
# Format: list of dicts {"bond": (parent, child), "atoms": [child, ...]}
# where "atoms" is the heavy-atom set rigidly attached to `child`
# (e.g., for ARG chi1 BRANCH: bond=(CA,CB), atoms=[CB]; chi2: bond=(CB,CG),
# atoms=[CG]; chi3: bond=(CG,CD), atoms=[CD]; chi4: bond=(CD,NE),
# atoms=[NE,CZ,NH1,NH2]).
#
# H atoms attached to a spine atom go with that atom's BRANCH.

SIDECHAIN_TEMPLATES: dict[str, list[dict]] = {
    "ALA": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
    ],
    "ARG": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"), "atoms": ["CG"]},
        {"bond": ("CG", "CD"), "atoms": ["CD"]},
        {"bond": ("CD", "NE"), "atoms": ["NE", "CZ", "NH1", "NH2"]},
    ],
    "ASN": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"), "atoms": ["CG", "OD1", "ND2"]},
    ],
    "ASP": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"), "atoms": ["CG", "OD1", "OD2"]},
    ],
    "CYS": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "SG"), "atoms": ["SG"]},
    ],
    "GLN": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"), "atoms": ["CG"]},
        {"bond": ("CG", "CD"), "atoms": ["CD", "OE1", "NE2"]},
    ],
    "GLU": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"), "atoms": ["CG"]},
        {"bond": ("CG", "CD"), "atoms": ["CD", "OE1", "OE2"]},
    ],
    "GLY": [],  # no side chain, cannot be flex
    "HIS": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"), "atoms": ["CG", "ND1", "CD2", "CE1", "NE2"]},
    ],
    "ILE": [
        {"bond": ("CA", "CB"), "atoms": ["CB", "CG2"]},
        {"bond": ("CB", "CG1"), "atoms": ["CG1"]},
        {"bond": ("CG1", "CD1"), "atoms": ["CD1"]},
    ],
    "LEU": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"), "atoms": ["CG", "CD1", "CD2"]},
    ],
    "LYS": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"), "atoms": ["CG"]},
        {"bond": ("CG", "CD"), "atoms": ["CD"]},
        {"bond": ("CD", "CE"), "atoms": ["CE"]},
        {"bond": ("CE", "NZ"), "atoms": ["NZ"]},
    ],
    "MET": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"), "atoms": ["CG"]},
        {"bond": ("CG", "SD"), "atoms": ["SD"]},
        {"bond": ("SD", "CE"), "atoms": ["CE"]},
    ],
    "PHE": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"), "atoms": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]},
    ],
    "PRO": [],  # ring constraint; not rotated
    "SER": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "OG"), "atoms": ["OG"]},
    ],
    "THR": [
        {"bond": ("CA", "CB"), "atoms": ["CB", "CG2"]},
        {"bond": ("CB", "OG1"), "atoms": ["OG1"]},
    ],
    "TRP": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"),
         "atoms": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]},
    ],
    "TYR": [
        {"bond": ("CA", "CB"), "atoms": ["CB"]},
        {"bond": ("CB", "CG"),
         "atoms": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]},
    ],
    "VAL": [
        {"bond": ("CA", "CB"), "atoms": ["CB", "CG1", "CG2"]},
    ],
}

BACKBONE = {"N", "CA", "C", "O", "OXT"}


def _hydrogens_for(atom_name: str) -> set[str]:
    """Return set of H atom names that are typically bonded to atom_name."""
    h: set[str] = set()
    if atom_name == "CA":
        h.update({"HA", "HA1", "HA2", "HA3"})
    elif atom_name == "CB":
        h.update({"HB", "HB1", "HB2", "HB3", "1HB", "2HB", "3HB"})
    elif atom_name == "CG":
        h.update({"HG", "HG1", "HG2", "HG3", "1HG", "2HG", "3HG"})
    elif atom_name == "CG1":
        h.update({"HG11", "HG12", "HG13", "1HG1", "2HG1", "3HG1"})
    elif atom_name == "CG2":
        h.update({"HG21", "HG22", "HG23", "1HG2", "2HG2", "3HG2"})
    elif atom_name == "CD":
        h.update({"HD", "HD1", "HD2", "HD3", "1HD", "2HD", "3HD"})
    elif atom_name == "CD1":
        h.update({"HD11", "HD12", "HD13", "1HD1", "2HD1", "3HD1", "HD1"})
    elif atom_name == "CD2":
        h.update({"HD21", "HD22", "HD23", "1HD2", "2HD2", "3HD2", "HD2"})
    elif atom_name == "CE":
        h.update({"HE", "HE1", "HE2", "HE3", "1HE", "2HE", "3HE"})
    elif atom_name == "CE1":
        h.update({"HE1"})
    elif atom_name == "CE2":
        h.update({"HE2"})
    elif atom_name == "CE3":
        h.update({"HE3"})
    elif atom_name == "CH2":
        h.update({"HH2"})
    elif atom_name == "CZ":
        h.update({"HZ"})
    elif atom_name == "CZ2":
        h.update({"HZ2"})
    elif atom_name == "CZ3":
        h.update({"HZ3"})
    elif atom_name == "NE":
        h.update({"HE", "HNE"})
    elif atom_name == "NE1":
        h.update({"HE1", "HNE1"})
    elif atom_name == "NE2":
        h.update({"HE2", "HNE2", "1HE2", "2HE2"})
    elif atom_name == "ND1":
        h.update({"HD1", "HND1"})
    elif atom_name == "ND2":
        h.update({"HD2", "1HD2", "2HD2"})
    elif atom_name == "NH1":
        h.update({"HH11", "HH12", "1HH1", "2HH1"})
    elif atom_name == "NH2":
        h.update({"HH21", "HH22", "1HH2", "2HH2"})
    elif atom_name == "NZ":
        h.update({"HZ1", "HZ2", "HZ3", "1HZ", "2HZ", "3HZ"})
    elif atom_name == "OH":
        h.update({"HH"})
    elif atom_name == "OG":
        h.update({"HG"})
    elif atom_name == "OG1":
        h.update({"HG1"})
    elif atom_name == "SG":
        h.update({"HG"})
    elif atom_name == "SD":
        pass  # no H on Met SD
    return h


def _parse_pdbqt_atoms(pdbqt_text: str) -> list[dict]:
    """Parse ATOM lines of a PDBQT string into structured records."""
    atoms = []
    for ln in pdbqt_text.splitlines():
        if not ln.startswith("ATOM"):
            continue
        try:
            serial = int(ln[6:11])
            aname = ln[12:16].strip()
            resname = ln[17:20].strip()
            chain = ln[21]
            resnum = int(ln[22:26])
            x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            q = float(ln[66:76])
            atype = ln[77:79].strip() if len(ln) >= 79 else ln[76:].strip()
        except (ValueError, IndexError):
            continue
        atoms.append({
            "serial": serial, "name": aname, "resname": resname,
            "chain": chain, "resnum": resnum,
            "x": x, "y": y, "z": z, "q": q, "type": atype,
            "raw": ln,
        })
    return atoms


def _atom_line(rec: dict, new_serial: int) -> str:
    return (
        f"ATOM  {new_serial:>5d} {rec['name']:>4s} {rec['resname']:>3s} "
        f"{rec['chain']}{rec['resnum']:>4d}    "
        f"{rec['x']:8.3f}{rec['y']:8.3f}{rec['z']:8.3f}"
        f"  0.00  0.00    {rec['q']:+7.3f} {rec['type']:>2s}"
    )


def split_clean_pdbqt(
    pdbqt_text: str,
    flex_residues: list[tuple[str, int]],
) -> tuple[str, str, list[str]]:
    """Split a clean (chains A/B preserved) PDBQT into (rigid_str, flex_str).

    flex_residues: list of (chain, resnum) tuples to make flexible.

    Returns (rigid_str, flex_str, warnings).
    """
    atoms = _parse_pdbqt_atoms(pdbqt_text)
    warnings: list[str] = []

    flex_set = set(flex_residues)
    # Group atoms by (chain, resnum) for flex residues
    by_res: dict[tuple[str, int], list[dict]] = {}
    rigid_atoms: list[dict] = []
    for a in atoms:
        key = (a["chain"], a["resnum"])
        if key in flex_set:
            by_res.setdefault(key, []).append(a)
        else:
            rigid_atoms.append(a)

    # For each flex residue, also keep the backbone (N, CA, C, O, H, HA) in rigid
    flex_blocks: list[str] = []
    for (chain, resnum) in flex_residues:
        residue_atoms = by_res.get((chain, resnum), [])
        if not residue_atoms:
            warnings.append(f"flex residue {chain}:{resnum} not found in receptor")
            continue
        resname = residue_atoms[0]["resname"]
        template = SIDECHAIN_TEMPLATES.get(resname)
        if template is None:
            warnings.append(f"no template for {resname} ({chain}:{resnum}); skipping")
            # keep all atoms in rigid
            rigid_atoms.extend(residue_atoms)
            continue
        if not template:
            warnings.append(f"residue {resname} {chain}:{resnum} has no rotatable side chain; keeping rigid")
            rigid_atoms.extend(residue_atoms)
            continue

        # Build name -> atom map
        by_name = {a["name"]: a for a in residue_atoms}

        # Backbone atoms stay rigid (N, C, O, H1/H2/H3 if any, plus CA also
        # appears here so the residue stays anchored in the protein).
        for an in BACKBONE | {"H", "H1", "H2", "H3", "HN", "HA"}:
            if an in by_name:
                rigid_atoms.append(by_name[an])

        # CA atom is needed in the flex block too (as ROOT).
        ca = by_name.get("CA")
        if ca is None:
            warnings.append(f"no CA atom for {chain}:{resnum}; skipping flex")
            rigid_atoms.extend([a for a in residue_atoms if a["name"] not in BACKBONE])
            continue

        # Walk the template and build BRANCH tree.
        # We need atom serial numbers consistent within the flex block.
        # We'll renumber serials at write time.
        flex_block_atoms: list[tuple[str, dict, int]] = []
        # tag each placement with depth so we can write open/close BRANCHes
        # Simple linear walk: for each chi step, write BRANCH parent_serial child_serial
        # then ATOM child + any auxiliary atoms attached (no further rotation).
        # When done with all, write reverse ENDBRANCHes.

        # First emit ROOT (with CA only).
        lines: list[str] = []
        lines.append("ROOT")
        ca_idx = 1  # CA gets local serial 1 in flex residue
        flex_block_atoms.append(("ATOM", ca, ca_idx))
        lines.append(("ATOM", ca, ca_idx))  # placeholder; we'll render later
        lines.append("ENDROOT")

        next_serial = ca_idx + 1
        # parent name -> serial map
        name_to_serial: dict[str, int] = {"CA": ca_idx}
        branch_stack: list[tuple[int, int]] = []  # (parent_serial, child_serial)

        skipped_arm = False
        for arm in template:
            parent_name, child_name = arm["bond"]
            atom_names = arm["atoms"]
            parent_atom = by_name.get(parent_name)
            child_atom = by_name.get(child_name)
            if parent_atom is None or child_atom is None:
                warnings.append(f"{chain}:{resnum} {resname}: missing {parent_name} or {child_name}; truncating arm")
                skipped_arm = True
                break
            parent_serial = name_to_serial.get(parent_name)
            if parent_serial is None:
                warnings.append(f"{chain}:{resnum} {resname}: parent {parent_name} not yet placed; truncating")
                skipped_arm = True
                break
            child_serial = next_serial
            next_serial += 1
            name_to_serial[child_name] = child_serial
            lines.append(("BRANCH", parent_serial, child_serial))
            branch_stack.append((parent_serial, child_serial))
            lines.append(("ATOM", child_atom, child_serial))
            # Attach H to child
            for hname in _hydrogens_for(child_name):
                hatom = by_name.get(hname)
                if hatom is not None and hatom["name"] not in name_to_serial:
                    s = next_serial; next_serial += 1
                    name_to_serial[hname] = s
                    lines.append(("ATOM", hatom, s))
            # Attach other atoms in this arm (heavy atoms in same rigid block as child)
            for an in atom_names:
                if an == child_name:
                    continue
                aatom = by_name.get(an)
                if aatom is None:
                    warnings.append(f"{chain}:{resnum} {resname}: missing {an} in arm")
                    continue
                if an in name_to_serial:
                    continue
                s = next_serial; next_serial += 1
                name_to_serial[an] = s
                lines.append(("ATOM", aatom, s))
                # H on this atom
                for hname in _hydrogens_for(an):
                    hatom = by_name.get(hname)
                    if hatom is not None and hatom["name"] not in name_to_serial:
                        sh = next_serial; next_serial += 1
                        name_to_serial[hname] = sh
                        lines.append(("ATOM", hatom, sh))

        # Close all branches in reverse
        for parent_serial, child_serial in reversed(branch_stack):
            lines.append(("ENDBRANCH", parent_serial, child_serial))

        # Render the flex block
        torsion_count = len(branch_stack)
        block_lines = [
            f"BEGIN_RES {resname} {chain}{resnum:>4d}",
            f"REMARK   {torsion_count} active torsions:",
        ]
        for ln in lines:
            if isinstance(ln, str):
                block_lines.append(ln)
            elif ln[0] == "BRANCH":
                _, ps, cs = ln
                block_lines.append(f"BRANCH {ps:>3d} {cs:>3d}")
            elif ln[0] == "ENDBRANCH":
                _, ps, cs = ln
                block_lines.append(f"ENDBRANCH {ps:>3d} {cs:>3d}")
            elif ln[0] == "ATOM":
                _, atom_rec, s = ln
                block_lines.append(_atom_line(atom_rec, s))
        block_lines.append(f"END_RES {resname} {chain}{resnum:>4d}")
        flex_blocks.append("\n".join(block_lines))

    # Build rigid output: original PDBQT minus flex side-chain atoms
    # Keep all REMARK/HEADER lines from input, then renumber rigid atoms.
    out_rigid_lines: list[str] = []
    for ln in pdbqt_text.splitlines():
        if ln.startswith("ATOM"):
            continue
        out_rigid_lines.append(ln)
    # Append rigid atoms in original order
    rigid_atoms.sort(key=lambda a: a["serial"])
    serial = 1
    for a in rigid_atoms:
        out_rigid_lines.append(_atom_line(a, serial))
        serial += 1
    out_rigid_lines.append("TER")

    rigid_str = "\n".join(out_rigid_lines) + "\n"
    flex_str = "\n".join(flex_blocks) + "\n" if flex_blocks else ""

    return rigid_str, flex_str, warnings


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 4:
        print("Usage: flex_split.py <input.pdbqt> <out_rigid.pdbqt> <out_flex.pdbqt> chain:resnum [chain:resnum...]")
        sys.exit(1)
    inp = Path(sys.argv[1])
    out_rigid = Path(sys.argv[2])
    out_flex = Path(sys.argv[3])
    flex_res = []
    for tok in sys.argv[4:]:
        c, n = tok.split(":")
        flex_res.append((c, int(n)))
    rigid_str, flex_str, warnings = split_clean_pdbqt(inp.read_text(), flex_res)
    out_rigid.write_text(rigid_str)
    out_flex.write_text(flex_str)
    for w in warnings:
        print(f"WARN: {w}")
    print(f"wrote rigid={out_rigid} flex={out_flex}")
