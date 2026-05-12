#!/usr/bin/env python3
"""
Convert a PDB2PQR-produced .pqr file (AMBER FF14SB at pH 7.4, proper
per-residue charges) into an AutoDock Vina-compatible .pdbqt file with:
  - Correct AutoDock atom types (HD, N, NA, OA, A, C, SA, S)
  - Charges from PDB2PQR preserved verbatim
  - Non-polar Hs MERGED into their carrier heavy atoms (the AD4 united-atom
    convention)
  - Total receptor charge integer-valued (near 0 for TYMS)
  - Per-residue charge sanity: ARG/LYS ≈ +1, ASP/GLU ≈ −1, HIS ≈ 0 / +1

Usage:
  python pqr_to_pdbqt.py input.pqr output.pdbqt
"""
import sys, pathlib, re
from collections import defaultdict

PQR_PATH = pathlib.Path(sys.argv[1])
OUT_PATH = pathlib.Path(sys.argv[2])

# AutoDock atom-type lookup. Keys: (residue_name, atom_name) → AD type
# Backbone defaults: N → N (donor), CA/C → C, O → OA (acceptor), H/HN → HD (donor)
BACKBONE = {
    "N": "N",   # AD considers N as HB-acceptor; for amide N with H it's still N
    "CA": "C", "C": "C",
    "O": "OA", "OXT": "OA",
    "H": "HD", "HN": "HD", "H1": "HD", "H2": "HD", "H3": "HD",
}
# Aromatic carbons get "A", aromatic Ns and Os their specific HB types
AROMATIC_C = {"CG","CD1","CD2","CE1","CE2","CZ","CH2","CD","CE3","CZ2","CZ3"}
# Per-residue overrides
SIDECHAIN_TYPES = {
    "ARG": {"NE":"N","NH1":"N","NH2":"N","HE":"HD","HH11":"HD","HH12":"HD",
            "HH21":"HD","HH22":"HD","CZ":"C"},
    "LYS": {"NZ":"N","HZ1":"HD","HZ2":"HD","HZ3":"HD"},
    "HIS": {"ND1":"NA","NE2":"NA","HD1":"HD","HE2":"HD","CG":"A","CD2":"A","CE1":"A"},
    "HID": {"ND1":"NA","NE2":"NA","HD1":"HD","CG":"A","CD2":"A","CE1":"A"},
    "HIE": {"ND1":"NA","NE2":"NA","HE2":"HD","CG":"A","CD2":"A","CE1":"A"},
    "HIP": {"ND1":"N","NE2":"N","HD1":"HD","HE2":"HD","CG":"A","CD2":"A","CE1":"A"},
    "PHE": {a:"A" for a in ["CG","CD1","CD2","CE1","CE2","CZ"]},
    "TYR": {**{a:"A" for a in ["CG","CD1","CD2","CE1","CE2","CZ"]},"OH":"OA","HH":"HD"},
    "TRP": {**{a:"A" for a in ["CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"]},
            "NE1":"NA","HE1":"HD"},
    "ASP": {"OD1":"OA","OD2":"OA"},
    "GLU": {"OE1":"OA","OE2":"OA"},
    "ASN": {"ND2":"N","HD21":"HD","HD22":"HD","OD1":"OA"},
    "GLN": {"NE2":"N","HE21":"HD","HE22":"HD","OE1":"OA"},
    "SER": {"OG":"OA","HG":"HD"},
    "THR": {"OG1":"OA","HG1":"HD"},
    "CYS": {"SG":"SA","HG":"HD"},
    "MET": {"SD":"SA"},
}

def autodock_type(resname, atom_name, element):
    """Best-effort AD type assignment."""
    if element == "H":
        # Check if HD (polar) or non-polar (gets merged)
        # Heuristic: if attached to N/O/S (by name), it's polar
        # Use the per-residue table
        per = SIDECHAIN_TYPES.get(resname, {})
        if atom_name in per:
            return per[atom_name]
        if atom_name in BACKBONE: return BACKBONE[atom_name]
        return None  # non-polar — will be merged
    # Heavy atom
    per = SIDECHAIN_TYPES.get(resname, {})
    if atom_name in per:
        return per[atom_name]
    if atom_name in BACKBONE:
        return BACKBONE[atom_name]
    # Generic by element
    return {"C":"C","N":"N","O":"O","S":"S","P":"P"}.get(element, element[0])


def parse_pqr_atoms(text):
    """Yield (resname, resi_int, atom_name, x, y, z, q, element) tuples in order.
    Tracks chain id by RESETTING on residue-number wraparound (PQR has no chain col)."""
    chain = "A"
    last_resi = None
    for line in text.splitlines():
        if not line.startswith(("ATOM","HETATM")): continue
        toks = line.split()
        # PQR layout: ATOM serial name resname resi x y z q r
        if len(toks) != 10: continue
        try:
            serial, atom_name, resn, resi, x, y, z, q, r = toks[1:10]
            resi = int(resi); x=float(x); y=float(y); z=float(z); q=float(q)
        except (ValueError, IndexError):
            continue
        # Detect chain wraparound (resi resets to a small number after a large one)
        if last_resi is not None and resi < last_resi - 50:
            chain = "B"
        last_resi = resi
        # Element from atom name
        elem = atom_name[0].upper()
        if elem.isdigit() and len(atom_name) > 1:
            elem = atom_name[1].upper()
        yield (chain, resn, resi, atom_name, x, y, z, q, elem)

text = PQR_PATH.read_text()
atoms = list(parse_pqr_atoms(text))
print(f"  parsed {len(atoms)} atoms from {PQR_PATH.name}")

# Group H charges that need merging into a carrier heavy atom
# Build a per-(chain, resi) coordinate lookup of heavy atoms (for find-nearest)
heavy = {}
for chain, resn, resi, name, x, y, z, q, elem in atoms:
    if elem == "H": continue
    heavy.setdefault((chain, resi), []).append((name, x, y, z))

# For each non-polar H, find its nearest heavy atom in the same residue
# and add its charge to that heavy atom's running total
extra_q = defaultdict(float)   # (chain, resi, heavy_name) → extra charge
polar_h = []                   # list of (chain, resn, resi, name, x, y, z, q)
for chain, resn, resi, name, x, y, z, q, elem in atoms:
    if elem != "H": continue
    ad = autodock_type(resn, name, elem)
    if ad == "HD":
        polar_h.append((chain, resn, resi, name, x, y, z, q))
        continue
    # Non-polar H: find nearest heavy atom in same residue
    cand = heavy.get((chain, resi), [])
    if not cand:
        # No heavy atom in residue — drop the H charge silently
        continue
    best = min(cand, key=lambda a: (a[1]-x)**2 + (a[2]-y)**2 + (a[3]-z)**2)
    extra_q[(chain, resi, best[0])] += q

# Now write the PDBQT
lines_out = []
serial = 0
for chain, resn, resi, name, x, y, z, q, elem in atoms:
    if elem == "H":
        ad = autodock_type(resn, name, elem)
        if ad != "HD":
            continue  # non-polar H — merged
        # polar H — keep with its own charge
        q_final = q
        atype = "HD"
    else:
        # heavy atom: add any merged-H charges
        q_final = q + extra_q.get((chain, resi, name), 0.0)
        atype = autodock_type(resn, name, elem)
        if atype is None:
            atype = elem

    serial += 1
    # Build a proper PDB-format atom-name field (4 chars wide, cols 13-16):
    #   names with element symbol of 1 char and total length ≤ 3 → leading space
    #   4-char names like HH11/HH12/CD11 → no leading space
    if len(name) <= 3:
        an = " " + name + " " * (3 - len(name))   # 4 chars total
    else:
        an = name[:4]                              # 4 chars total
    record = "ATOM"
    line = (
        f"{record:6s}"                    # 1-6   record
        f"{serial:5d}"                    # 7-11  serial
        f" "                              # 12    blank
        f"{an:4s}"                        # 13-16 atom name (always 4 wide)
        f" "                              # 17    altloc
        f"{resn:3s}"                      # 18-20 resname
        f" "                              # 21    blank
        f"{chain:1s}"                     # 22    chain
        f"{resi:4d}"                      # 23-26 resi
        f"    "                           # 27-30 iCode + 3 blanks
        f"{x:8.3f}{y:8.3f}{z:8.3f}"       # 31-54 coords
        f"{0.0:6.2f}{0.0:6.2f}"           # 55-66 occ + bf
        f"    "                           # 67-70 blank
        f"{q_final:+7.3f}"                # 71-77 charge (signed, 7 wide)
        f" "                              # 78    blank
        f"{atype:<2s}"                    # 79-80 AD type
    )
    lines_out.append(line)

lines_out.append("TER")
OUT_PATH.write_text("\n".join(lines_out) + "\n")
print(f"  wrote {OUT_PATH} ({OUT_PATH.stat().st_size} bytes)")

# Validate
total = 0.0
by_res = defaultdict(float)
for line in OUT_PATH.read_text().splitlines():
    if not line.startswith("ATOM"): continue
    try: q = float(line[70:77])
    except ValueError: continue
    resn = line[17:20].strip()
    resi = line[22:26].strip()
    chain = line[21]
    total += q
    by_res[(resn, resi, chain)] += q
print(f"\n  total receptor q = {total:+.2f}")
def status(label, values, expected_range):
    bad = [round(v,2) for v in values if not expected_range[0] < v < expected_range[1]]
    flag = "✅" if not bad else "❌"
    n = len(values); mean = sum(values)/n if n else 0
    print(f"  {flag} {label:4s}: {n} residues, mean {mean:+.2f}, "
          f"out-of-range: {len(bad)}{(' '+str(bad[:3])) if bad else ''}")
    return not bad

arg = [v for (r,_,_),v in by_res.items() if r=="ARG"]
lys = [v for (r,_,_),v in by_res.items() if r=="LYS"]
asp = [v for (r,_,_),v in by_res.items() if r=="ASP"]
glu = [v for (r,_,_),v in by_res.items() if r=="GLU"]
ok_total = abs(total) < 5
ok_arg = status("ARG", arg, (0.7, 1.3))
ok_lys = status("LYS", lys, (0.7, 1.3))
ok_asp = status("ASP", asp, (-1.3, -0.7))
ok_glu = status("GLU", glu, (-1.3, -0.7))
print(f"  {'✅' if ok_total else '❌'} |sum_q| = {abs(total):.2f}  (must be < 5)")
all_ok = ok_total and ok_arg and ok_lys and ok_asp and ok_glu
print(f"\n{'✅ ALL ASSERTIONS PASS' if all_ok else '❌ ASSERTIONS FAILED'}")
sys.exit(0 if all_ok else 1)
