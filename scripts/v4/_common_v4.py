"""Shared v4 helpers."""
import os, re, math, hashlib, subprocess
import numpy as np
from datetime import datetime

PROJECT = os.path.expanduser("~/conserved_site_project")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
OBABEL = os.environ.get("OBABEL", os.path.join(PROJECT, ".venv/bin/obabel"))
VINA = os.environ.get("VINA", "/opt/homebrew/bin/vina")
VENV_PY = os.path.join(PROJECT, ".venv/bin/python")
MK_PREP_REC = os.path.join(PROJECT, ".venv/bin/mk_prepare_receptor.py")


def md5_of(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def parse_vina_pdbqt(pdbqt):
    """Return list of all REMARK VINA RESULT affinities."""
    affs = []
    if not os.path.exists(pdbqt):
        return affs
    with open(pdbqt) as f:
        for line in f:
            m = re.match(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", line)
            if m:
                affs.append(float(m.group(1)))
    return affs


def parse_pdbqt_models(pdbqt):
    """Return list of models, each a list of (atom_name, coord_array) for heavy atoms."""
    models = []
    cur = []
    with open(pdbqt) as f:
        for line in f:
            if line.startswith("MODEL"):
                cur = []
            elif line.startswith("ENDMDL"):
                models.append(cur)
                cur = []
            elif line.startswith(("ATOM", "HETATM")):
                try:
                    nm = line[12:16].strip()
                    elem = line[77:79].strip() if len(line) >= 79 else ""
                    if elem in ("HD", "H"):
                        continue
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    cur.append((nm, np.array([x, y, z])))
                except Exception:
                    pass
    return models


def native_heavy(pdb):
    """Return list of (atom_name, coord) for heavy atoms in a PDB ligand."""
    atoms = []
    with open(pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                name = line[12:16].strip()
                elem = line[76:78].strip() or name[0]
                if elem == "H":
                    continue
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                atoms.append((name, np.array([x, y, z])))
    return atoms


def rmsd_by_name(pose_atoms, native_atoms):
    """Compute RMSD by matching atom names. Returns NaN if too few matches."""
    nat_map = {n: c for n, c in native_atoms}
    matched = [(c, nat_map[n]) for n, c in pose_atoms if n in nat_map]
    if len(matched) < max(3, len(native_atoms) // 2):
        return float("nan")
    sq = [np.sum((a - b) ** 2) for a, b in matched]
    return float(math.sqrt(sum(sq) / len(sq)))


def rmsd_nearest(pose_atoms, native_atoms):
    """Fallback RMSD by nearest-neighbor matching."""
    pose_arr = np.array([c for _, c in pose_atoms])
    nat_arr = np.array([c for _, c in native_atoms])
    if len(pose_arr) == 0 or len(nat_arr) == 0:
        return float("nan")
    used = set()
    pairs = []
    for i in range(len(pose_arr)):
        best = None; best_d = 1e9
        for j in range(len(nat_arr)):
            if j in used:
                continue
            d = np.linalg.norm(pose_arr[i] - nat_arr[j])
            if d < best_d:
                best_d = d; best = j
        if best is not None:
            used.add(best)
            pairs.append((pose_arr[i], nat_arr[best]))
    if not pairs:
        return float("nan")
    sq = [np.sum((a - b) ** 2) for a, b in pairs]
    return float(math.sqrt(sum(sq) / len(pairs)))


def rmsd_top(pose_atoms, native_atoms):
    """RMSD by name first; fall back to nearest if names don't match."""
    nat_map = {n: c for n, c in native_atoms}
    n_match = sum(1 for n, _ in pose_atoms if n in nat_map)
    if n_match >= max(3, len(native_atoms) // 2):
        return rmsd_by_name(pose_atoms, native_atoms)
    return rmsd_nearest(pose_atoms, native_atoms)


def split_top(out_pdbqt, top_pdbqt):
    """Extract the first MODEL from a multi-model PDBQT."""
    with open(out_pdbqt) as f, open(top_pdbqt, "w") as g:
        in_first = False
        for line in f:
            if line.startswith("MODEL 1"):
                in_first = True
            if in_first:
                g.write(line)
                if line.startswith("ENDMDL"):
                    break


def _heavy_records_from_pdbqt(pdbqt_path):
    """Return list of (x, y, z, elem) for heavy (non H/HD) atoms in PDBQT, in file order."""
    out = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                elem = line[77:79].strip() if len(line) >= 79 else ""
                if elem in ("HD", "H"):
                    continue
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                out.append((x, y, z, elem))
    return out


def _heavy_names_from_pdb(pdb_path):
    """Return list of (atom_name_field, resname, element, x, y, z) for heavy atoms."""
    out = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                elem = line[76:78].strip()
                if not elem:
                    nm = line[12:16].strip()
                    elem = nm[0]
                if elem == "H":
                    continue
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                out.append((line[12:16], line[17:20].strip(), elem, x, y, z))
    return out


def restore_atom_names(top_pdbqt, ref_ligand_pdb, out_pdb,
                       reference_dock_pdbqt=None):
    """Take the top PDBQT pose and restore atom names by mapping via the
    reference (pre-docking) ligand PDBQT and the original named PDB.

    Strategy:
    - Build a mapping from each pre-docking PDBQT heavy-atom INDEX to the matching
      atom in the original named PDB by minimum Euclidean distance (the pre-dock
      PDBQT shares the same coordinates as the input PDB; obabel only re-orders).
    - Then walk the docked PDBQT and use that index->name mapping to write a PDB
      with proper atom names.
    """
    if reference_dock_pdbqt is None:
        # Default: assume the pre-dock PDBQT lives next to the docked one in the
        # standard pipeline layout (caller should pass it explicitly to be safe)
        reference_dock_pdbqt = ref_ligand_pdb.replace(".pdb", ".pdbqt")
        if not os.path.exists(reference_dock_pdbqt):
            # Fallback: try the canonical 05b_ligand_v2/dump.pdbqt
            reference_dock_pdbqt = os.path.join(
                PROJECT, "05b_ligand_v2", "dump.pdbqt")

    # Build index -> (name, element) mapping using the pre-dock PDBQT
    pre_heavy = _heavy_records_from_pdbqt(reference_dock_pdbqt)
    ref_atoms = _heavy_names_from_pdb(ref_ligand_pdb)

    # For each pre-dock heavy atom, find nearest by coords in the named PDB
    idx_to_name = []
    used = set()
    for (px, py, pz, pe) in pre_heavy:
        best = None; best_d = 1e9
        for j, (name, resname, elem, rx, ry, rz) in enumerate(ref_atoms):
            if j in used:
                continue
            d = (rx - px) ** 2 + (ry - py) ** 2 + (rz - pz) ** 2
            if d < best_d:
                best_d = d; best = j
        if best is not None:
            used.add(best)
            idx_to_name.append((ref_atoms[best][0], ref_atoms[best][1], ref_atoms[best][2]))
        else:
            idx_to_name.append((f" {pe:<3s}", "UNK", pe))

    # Now write out the docked pose with these names
    docked_heavy = _heavy_records_from_pdbqt(top_pdbqt)
    n = min(len(docked_heavy), len(idx_to_name))
    out_lines = []
    for i in range(n):
        x, y, z, _elem = docked_heavy[i]
        name_field, resname, e = idx_to_name[i]
        out_line = (f"HETATM{i+1:5d} {name_field} {resname:<3s} X 414    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {e:>2s}\n")
        out_lines.append(out_line)
    out_lines.append("END\n")
    with open(out_pdb, "w") as f:
        f.writelines(out_lines)
    return out_pdb


def receptor_max_abs_charge(pdbqt_path):
    if not os.path.exists(pdbqt_path):
        return 0.0
    m = 0.0
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    c = float(line[66:76].strip())
                    if abs(c) > m:
                        m = abs(c)
                except (ValueError, IndexError):
                    pass
    return m


def crystal_dump_centroid(ligand_pdb):
    coords = []
    with open(ligand_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                elem = line[76:78].strip() or line[12:16].strip()[0]
                if elem == "H":
                    continue
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                coords.append([x, y, z])
    arr = np.array(coords)
    return arr.mean(axis=0)


def make_holo_dimer(prot_h, cof_a, cof_b, out_pdb):
    with open(out_pdb, "w") as out:
        for src in [prot_h, cof_a, cof_b]:
            with open(src) as f:
                for line in f:
                    if line.startswith(("ATOM  ", "HETATM", "TER")):
                        out.write(line)
            out.write("TER\n")
        out.write("END\n")
    return out_pdb


def prepare_receptor_with_charges(prot_in, prot_out, label, log_fn=print):
    """Try Gasteiger via obabel first, then meeko."""
    proc = subprocess.run(
        [OBABEL, prot_in, "-O", prot_out, "-xr", "-p", "7.4", "--partialcharge", "gasteiger"],
        capture_output=True, text=True
    )
    if proc.returncode == 0 and os.path.exists(prot_out):
        m = receptor_max_abs_charge(prot_out)
        log_fn(f"  [{label}] obabel-gasteiger rc={proc.returncode} max|q|={m:.3f}")
        if m > 0.05:
            return True, "obabel_gasteiger"

    if os.path.exists(MK_PREP_REC):
        base = prot_out.rsplit(".", 1)[0]
        proc = subprocess.run([VENV_PY, MK_PREP_REC, "-i", prot_in, "-o", base],
                              capture_output=True, text=True)
        log_fn(f"  [{label}] meeko rc={proc.returncode}")
        if proc.returncode == 0 and os.path.exists(prot_out):
            m = receptor_max_abs_charge(prot_out)
            log_fn(f"  [{label}] meeko max|q|={m:.3f}")
            if m > 0.05:
                return True, "meeko"

    return False, None
