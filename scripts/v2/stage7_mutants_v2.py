#!/usr/bin/env python3
"""Stage 7 v2: Mutant panel x (apo, holo).
- Ala-scan singles, opposite-class singles, double mutants, Arg-clamp probes, control T170A.
- PyMOL Mutagenesis Wizard with rotamer-strain selection.
- G217W clash check (drop if heavy-atom clash <1.8 A).
- Re-dock with Vina (18^3 box, exhaustiveness 32, parsing from PDBQT REMARK lines).
- Build viewer files: <mut>_<cond>_complex.pdb and <mut>_<cond>_top_pose.pdb.
"""
import os, sys, subprocess, json, re, math, csv
from datetime import datetime
import numpy as np
from Bio.PDB import PDBParser

# Reuse helpers
sys.path.insert(0, os.path.join(os.path.expanduser("~/conserved_site_project"), "scripts/v2"))
from stage6_dock_wt_v2 import (parse_vina_result_pdbqt, parse_pdbqt_models,
                                native_heavy_coords, rmsd_to_native, prepare_receptor,
                                compute_centroid_chainA, split_top_pose_pdbqt,
                                make_holo_dimer_with_cofactor)

PROJECT = os.path.expanduser("~/conserved_site_project")
MUT_DIR = os.path.join(PROJECT, "07b_mut_docking_v2")
VIEWER_DIR = os.path.join(MUT_DIR, "viewer_files")
WT_DIR = os.path.join(PROJECT, "06b_docking_wt_v2")
LIG_DIR = os.path.join(PROJECT, "05b_ligand_v2")
STR_DIR = os.path.join(PROJECT, "03b_structure_v2")
AS_DIR = os.path.join(PROJECT, "02b_active_site_v2")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v2_07_mutants.log")
PYMOL = os.environ.get("PYMOL", "/opt/homebrew/bin/pymol")
VINA = os.environ.get("VINA", "/opt/homebrew/bin/vina")
OBABEL = os.environ.get("OBABEL", "/Users/ario/conserved_site_project/.venv/bin/obabel")

AA1to3 = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
          "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
          "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V2] STAGE7: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def build_mutation_panel(selected, ref_seq):
    """Build the v2 panel per spec.
    selected = list of selected residue ints (chain A numbering = ref).
    ref_seq = ungapped human TYMS sequence (1-indexed in code via -1).
    Returns list of dicts with id, mutations (list of (resi, new_aa3)), category."""
    panel = []

    # 1) Alanine scan of selected residues (skip those already A or G)
    ala_targets = [p for p in selected if ref_seq[p-1] not in ("A", "G")]
    # Cap at 8
    for p in ala_targets[:8]:
        wt = ref_seq[p-1]
        if wt == "A": continue
        mid = f"{wt}{p}A"
        panel.append({"id": mid, "mutations": [(p, "ALA")], "category": "ala_scan"})

    # 2) Chemically-opposite singles (7) — chemistry-meaningful
    opp_map = {"R": "E", "C": "S", "H": "A", "N": "D", "Y": "F",
               "K": "E", "D": "K", "E": "K"}
    seen = set([m["id"] for m in panel])
    opp_count = 0
    for p in selected:
        if opp_count >= 7: break
        wt = ref_seq[p-1]
        if wt not in opp_map: continue
        new = opp_map[wt]
        mid = f"{wt}{p}{new}"
        if mid in seen: continue
        panel.append({"id": mid, "mutations": [(p, AA1to3[new])], "category": "opposite"})
        seen.add(mid); opp_count += 1

    # 3) Doubles (5)
    doubles = [
        ("C195A_H196A", [(195, "ALA"), (196, "ALA")], "double_dyad"),
        ("R175E_R176E", [(175, "GLU"), (176, "GLU")], "double_phosclamp"),
        ("C195S_H196N", [(195, "SER"), (196, "ASN")], "double_polar_neutral"),
        ("R215A_N226A", [(215, "ALA"), (226, "ALA")], "double_substrate_orient"),
        ("Y258F_F225Y", [(258, "PHE"), (225, "TYR")], "double_aromatic_swap"),
    ]
    for mid, muts, cat in doubles:
        panel.append({"id": mid, "mutations": muts, "category": cat})

    # 4) Arg-clamp probes (4) — per scientific officer
    for p in [50, 175, 176, 215]:
        wt = ref_seq[p-1]
        mid = f"{wt}{p}A"
        if mid not in seen:
            panel.append({"id": mid, "mutations": [(p, "ALA")], "category": "arg_clamp"})
            seen.add(mid)

    # 5) Control: T170A (or whatever residue 170 actually is)
    wt170 = ref_seq[170-1]
    mid = f"{wt170}170A"
    panel.append({"id": mid, "mutations": [(170, "ALA")], "category": "control_surface"})

    return panel


def pymol_mutate_with_rotamer_pick(in_pdb, out_pdb, mutations):
    """Apply mutations using PyMOL Mutagenesis Wizard with strain-based rotamer pick.
    For each mutation:
      - try set_rotamer i for i in 0..n_rotamers
      - keep the rotamer with lowest get_strain()
    Falls back to default rotamer if API not available.
    """
    cmds = [
        f"load {in_pdb}, prot",
        "wizard mutagenesis",
        "refresh_wizard",
    ]
    for resi, new_aa3 in mutations:
        cmds += [
            f'cmd.get_wizard().set_mode("{new_aa3}")',
            f'cmd.get_wizard().do_select("/prot//A/{resi}/")',
            # Pick lowest-strain rotamer
            "best_strain = 1e30",
            "best_idx = 0",
            "n_rot = 0",
            "try:",
            "    while True:",
            "        cmd.get_wizard().set_rotamer(n_rot)",
            "        try:",
            "            s = cmd.get_wizard().get_strain()",
            "        except Exception:",
            "            s = 1e30",
            "        if s < best_strain:",
            "            best_strain = s",
            "            best_idx = n_rot",
            "        n_rot += 1",
            "        if n_rot > 30:",
            "            break",
            "except Exception:",
            "    pass",
            "try:",
            "    cmd.get_wizard().set_rotamer(best_idx)",
            "except Exception:",
            "    pass",
            f'print("ROTAMER_PICK resi={resi} aa={new_aa3} best_idx=" + str(best_idx) + " strain=" + str(best_strain))',
            "cmd.get_wizard().apply()",
        ]
    cmds += ["set_wizard", f"save {out_pdb}, prot"]

    pml_path = out_pdb.replace(".pdb", ".pml")
    with open(pml_path, "w") as f:
        f.write("\n".join(cmds))
    proc = subprocess.run([PYMOL, "-cq", pml_path], capture_output=True, text=True)
    return proc


def check_clashes(pdb_path, resi_list, threshold=1.8):
    """For each mutated residue, check side-chain heavy-atom clashes < threshold A
    against all OTHER residues' heavy atoms. Returns list of clashing residue ids."""
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("p", pdb_path)
    target_atoms = []   # (resi, atom_obj)
    other_atoms = []    # atom_obj
    for chain in s[0]:
        if chain.id != "A": continue
        for res in chain:
            if res.id[0] != " ": continue
            ri = res.id[1]
            for atom in res:
                if atom.element == "H": continue
                if ri in resi_list and atom.get_name() not in ("N", "CA", "C", "O"):
                    target_atoms.append((ri, atom))
                elif ri not in resi_list:
                    other_atoms.append(atom)
    clashes = []
    for ri, ta in target_atoms:
        ta_coord = ta.get_coord()
        for oa in other_atoms:
            d = np.linalg.norm(ta_coord - oa.get_coord())
            if d < threshold:
                clashes.append((ri, ta.get_name(), oa.get_parent().id[1], oa.get_name(), d))
                break  # one clash is enough to flag this resi
    return clashes


def build_complex_pdb(receptor_pdb, top_pose_pdb, out_complex):
    """Concatenate receptor + ligand pose with TER separation and ligand on chain X."""
    with open(out_complex, "w") as out:
        with open(receptor_pdb) as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM", "TER")):
                    out.write(line)
        out.write("TER\n")
        with open(top_pose_pdb) as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM")):
                    new_line = line[:21] + "X" + line[22:]
                    out.write(new_line)
        out.write("TER\nEND\n")


def dock_one(rec_pdbqt, lig_pdbqt, out_pdbqt, log_file, centroid):
    cx, cy, cz = centroid
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", "18", "--size_y", "18", "--size_z", "18",
           "--exhaustiveness", "32", "--num_modes", "20", "--seed", "42",
           "--out", out_pdbqt, "--cpu", "4"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_file, "w") as f:
        f.write("STDOUT\n" + proc.stdout + "\nSTDERR\n" + proc.stderr)
    return proc


def main():
    os.makedirs(MUT_DIR, exist_ok=True)
    os.makedirs(VIEWER_DIR, exist_ok=True)
    open(STAGE_LOG, "w").close()
    log("Stage 7 v2 starting")

    selected = json.load(open(os.path.join(AS_DIR, "selected_meta.json")))["selected"]

    # Get reference sequence
    import requests
    rseq = requests.get("https://rest.uniprot.org/uniprotkb/P04818.fasta", timeout=30).text
    ref_seq = "".join(rseq.split("\n")[1:])
    log(f"ref seq len={len(ref_seq)}")

    panel = build_mutation_panel(selected, ref_seq)
    log(f"Built panel of {len(panel)} mutants")

    # Add G217W exploratory variant (per spec — even though not in selected)
    panel.append({"id": "G217W", "mutations": [(217, "TRP")], "category": "explore_g217w"})
    log(f"Final panel size: {len(panel)}")

    # Get WT centroid (already computed)
    prot_h = os.path.join(STR_DIR, "protein_dimer_h.pdb")
    cof_a_h = os.path.join(STR_DIR, "cofactor_chainA_h.pdb")
    cof_b_h = os.path.join(STR_DIR, "cofactor_chainB_h.pdb")
    centroid, _ = compute_centroid_chainA(prot_h, set(selected))
    log(f"centroid: {centroid}")
    lig_pdbqt = os.path.join(LIG_DIR, "dump.pdbqt")
    native_lig = os.path.join(STR_DIR, "ligand.pdb")
    wt_apo_aff = json.load(open(os.path.join(WT_DIR, "wt_apo.json")))["top_affinity"]
    wt_holo_aff = json.load(open(os.path.join(WT_DIR, "wt_holo.json")))["top_affinity"]

    rows = []
    skipped = []
    for mi, m in enumerate(panel):
        mid = m["id"]
        muts = m["mutations"]
        cat = m["category"]
        log(f"[{mi+1}/{len(panel)}] {mid} ({cat})")

        mut_subdir = os.path.join(MUT_DIR, mid)
        os.makedirs(mut_subdir, exist_ok=True)

        # 1) PyMOL mutate with rotamer pick
        mutated_pdb = os.path.join(mut_subdir, f"{mid}_mut.pdb")
        proc = pymol_mutate_with_rotamer_pick(prot_h, mutated_pdb, muts)
        if proc.returncode != 0 or not os.path.exists(mutated_pdb):
            log(f"  PyMOL mutation failed rc={proc.returncode}: {proc.stderr[:200]}")
            skipped.append((mid, "pymol_mutate_failed"))
            continue
        # extract rotamer pick info
        rot_info = [l for l in proc.stdout.split("\n") if "ROTAMER_PICK" in l]
        for ri in rot_info:
            log(f"  {ri}")

        # 2) Clash check
        resi_list = [r for r, _ in muts]
        clashes = check_clashes(mutated_pdb, resi_list, threshold=1.8)
        if clashes:
            log(f"  clashes detected: {clashes[:3]}...")
            if mid == "G217W":
                log(f"  G217W has heavy-atom clashes <1.8A — DROPPING per spec")
                skipped.append((mid, f"g217w_clash:{len(clashes)}"))
                continue
            else:
                log(f"  warning: {len(clashes)} clashes, but proceeding")

        # 3) Add hydrogens
        mutated_h = os.path.join(mut_subdir, f"{mid}_mut_h.pdb")
        subprocess.run([OBABEL, mutated_pdb, "-O", mutated_h, "-h"],
                      capture_output=True, text=True)

        # 4) Prep apo & holo receptors
        rec_apo_pdbqt = os.path.join(mut_subdir, f"{mid}_apo.pdbqt")
        rec_holo_pdbqt = os.path.join(mut_subdir, f"{mid}_holo.pdbqt")
        ok_apo = prepare_receptor(mutated_h, rec_apo_pdbqt, f"{mid}_apo")
        # Build holo with cofactors
        holo_pdb = os.path.join(mut_subdir, f"{mid}_holo.pdb")
        make_holo_dimer_with_cofactor(mutated_h, cof_a_h, cof_b_h, holo_pdb)
        ok_holo = prepare_receptor(holo_pdb, rec_holo_pdbqt, f"{mid}_holo")

        if not ok_apo:
            log(f"  apo receptor prep failed — skipping mutant")
            skipped.append((mid, "apo_prep_failed"))
            continue
        if not ok_holo:
            log(f"  holo receptor prep failed — using apo for holo")
            rec_holo_pdbqt = rec_apo_pdbqt

        # 5) Dock both conditions
        for cond, rec in [("apo", rec_apo_pdbqt), ("holo", rec_holo_pdbqt)]:
            out_pdbqt = os.path.join(mut_subdir, f"{mid}_{cond}.pdbqt")
            log_file = os.path.join(mut_subdir, f"{mid}_{cond}.log")
            proc = dock_one(rec, lig_pdbqt, out_pdbqt, log_file, centroid)
            if proc.returncode != 0:
                log(f"  vina {mid} {cond} rc={proc.returncode}: {proc.stderr[:200]}")
                rows.append({"mutant": mid, "category": cat, "condition": cond,
                            "top_affinity": float("nan"), "mean_top3": float("nan"),
                            "rmsd_to_native": float("nan"),
                            "ddG_vs_wt": float("nan"), "n_clashes": len(clashes),
                            "error": proc.stderr[:200]})
                continue

            affs = parse_vina_result_pdbqt(out_pdbqt)
            top = affs[0] if affs else float("nan")
            mean3 = float(np.mean(affs[:3])) if len(affs) >= 3 else float("nan")

            poses = parse_pdbqt_models(out_pdbqt)
            native = native_heavy_coords(native_lig)
            rmsd = rmsd_to_native(poses[0], native) if poses else float("nan")

            wt_aff = wt_apo_aff if cond == "apo" else wt_holo_aff
            ddg = top - wt_aff if not math.isnan(top) else float("nan")

            # Build viewer files
            top_pdbqt = os.path.join(mut_subdir, f"{mid}_{cond}_top.pdbqt")
            split_top_pose_pdbqt(out_pdbqt, top_pdbqt)
            top_pose_pdb = os.path.join(VIEWER_DIR, f"{mid}_{cond}_top_pose.pdb")
            subprocess.run([OBABEL, top_pdbqt, "-O", top_pose_pdb], capture_output=True, text=True)

            complex_pdb = os.path.join(VIEWER_DIR, f"{mid}_{cond}_complex.pdb")
            # Use mutated_h for receptor (apo) or holo_pdb (holo) so user sees same env
            recv_pdb = mutated_h if cond == "apo" else holo_pdb
            build_complex_pdb(recv_pdb, top_pose_pdb, complex_pdb)

            rows.append({"mutant": mid, "category": cat, "condition": cond,
                        "top_affinity": top, "mean_top3": mean3,
                        "rmsd_to_native": rmsd, "ddG_vs_wt": ddg,
                        "n_clashes": len(clashes),
                        "complex_pdb": complex_pdb, "top_pose_pdb": top_pose_pdb})
            log(f"  {mid} {cond}: top={top:.2f}, ddG={ddg:+.2f}, rmsd={rmsd:.2f}")

    # Write summary CSV
    fieldnames = ["mutant", "category", "condition", "top_affinity", "mean_top3",
                 "rmsd_to_native", "ddG_vs_wt", "n_clashes", "complex_pdb",
                 "top_pose_pdb", "error"]
    csv_path = os.path.join(MUT_DIR, "mutant_results_v2.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    log(f"Wrote {csv_path} ({len(rows)} rows)")

    # Summary JSON
    with open(os.path.join(MUT_DIR, "summary_v2.json"), "w") as f:
        json.dump({"n_panel": len(panel), "n_completed": len(rows),
                  "skipped": skipped,
                  "wt_apo": wt_apo_aff, "wt_holo": wt_holo_aff},
                  f, indent=2)

    # Verify viewer files load in PyMOL (sample one)
    if rows:
        sample = next((r for r in rows if "complex_pdb" in r and r.get("complex_pdb")), None)
        if sample:
            test_pml = os.path.join(MUT_DIR, "_verify_viewer.pml")
            with open(test_pml, "w") as f:
                f.write(f"load {sample['complex_pdb']}, c\nprint('VIEWER_OK')\n")
            proc = subprocess.run([PYMOL, "-cq", test_pml], capture_output=True, text=True)
            ok = "VIEWER_OK" in proc.stdout
            log(f"viewer file PyMOL load test: {'OK' if ok else 'FAIL'} (sample={os.path.basename(sample['complex_pdb'])})")

    log(f"Stage 7 v2 DONE: {len(rows)} successful runs, {len(skipped)} skipped")


if __name__ == "__main__":
    main()
