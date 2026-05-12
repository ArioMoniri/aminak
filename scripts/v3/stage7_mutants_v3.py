#!/usr/bin/env python3
"""Stage 7 v3: Mutant panel x (apo, holo) — fixed.

Fixes:
- FIX 1: receptors built with Gasteiger charges
- FIX 3: delta_vina_vs_wt = top_aff_mut - top_aff_wt (positive = destabilising)
- FIX 4: sculpt-based rotamer minimisation (avoiding broken get_strain in headless)
- FIX 5: mean_topk = mean(affinities[:min(3, len(affinities))]); n_modes column
- FIX 8: mis_docked flag for RMSD > 3.0 A
- FIX 9: rename column ddG -> delta_vina
"""
import os, sys, subprocess, json, re, math, csv
from datetime import datetime
import numpy as np
from Bio.PDB import PDBParser

sys.path.insert(0, os.path.join(os.path.expanduser("~/conserved_site_project"), "scripts/v3"))
from stage6_dock_wt_v3 import (parse_vina_pdbqt, parse_pdbqt_models, native_heavy,
                                rmsd_top, prepare_receptor_with_charges,
                                receptor_max_abs_charge, crystal_dump_centroid,
                                make_holo_dimer, split_top, reprotonate_cofactor)

PROJECT = os.path.expanduser("~/conserved_site_project")
MUT_DIR = os.path.join(PROJECT, "07c_mut_docking_v3")
VIEWER_DIR = os.path.join(MUT_DIR, "viewer_files")
WT_DIR = os.path.join(PROJECT, "06c_docking_wt_v3")
LIG_DIR = os.path.join(PROJECT, "05b_ligand_v2")
STR_DIR = os.path.join(PROJECT, "03b_structure_v2")
AS_DIR = os.path.join(PROJECT, "02b_active_site_v2")
LOG_DIR = os.path.join(PROJECT, "logs")
LOG_FILE = os.path.join(PROJECT, "pipeline.log")
STAGE_LOG = os.path.join(LOG_DIR, "v3_07_mutants.log")
PYMOL = os.environ.get("PYMOL", "/opt/homebrew/bin/pymol")
VINA = os.environ.get("VINA", "/opt/homebrew/bin/vina")
OBABEL = os.environ.get("OBABEL", "/Users/ario/conserved_site_project/.venv/bin/obabel")

AA1to3 = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
          "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
          "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL"}


def log(msg):
    line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] [V3] STAGE7: {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")
    with open(STAGE_LOG, "a") as f:
        f.write(line + "\n")


def build_panel(selected, ref_seq):
    """Same panel logic as v2, but never includes G217W (audit decision)."""
    panel = []
    ala_targets = [p for p in selected if ref_seq[p-1] not in ("A", "G")]
    for p in ala_targets[:8]:
        wt = ref_seq[p-1]
        if wt == "A": continue
        panel.append({"id": f"{wt}{p}A", "mutations": [(p, "ALA")], "category": "ala_scan"})

    opp_map = {"R":"E", "C":"S", "H":"A", "N":"D", "Y":"F", "K":"E", "D":"K", "E":"K"}
    seen = {m["id"] for m in panel}
    for p in selected:
        if len([m for m in panel if m["category"] == "opposite"]) >= 7: break
        wt = ref_seq[p-1]
        if wt not in opp_map: continue
        new = opp_map[wt]
        mid = f"{wt}{p}{new}"
        if mid in seen: continue
        panel.append({"id": mid, "mutations": [(p, AA1to3[new])], "category": "opposite"})
        seen.add(mid)

    doubles = [
        ("C195A_H196A", [(195, "ALA"), (196, "ALA")], "double_dyad"),
        ("R175E_R176E", [(175, "GLU"), (176, "GLU")], "double_phosclamp"),
        ("C195S_H196N", [(195, "SER"), (196, "ASN")], "double_polar_neutral"),
        ("R215A_N226A", [(215, "ALA"), (226, "ALA")], "double_substrate_orient"),
        ("Y258F_F225Y", [(258, "PHE"), (225, "TYR")], "double_aromatic_swap"),
    ]
    for mid, muts, cat in doubles:
        panel.append({"id": mid, "mutations": muts, "category": cat})

    for p in [50, 175, 176, 215]:
        wt = ref_seq[p-1]
        mid = f"{wt}{p}A"
        if mid not in seen:
            panel.append({"id": mid, "mutations": [(p, "ALA")], "category": "arg_clamp"})
            seen.add(mid)

    wt170 = ref_seq[170-1]
    panel.append({"id": f"{wt170}170A", "mutations": [(170, "ALA")], "category": "control_surface"})
    return panel


def pymol_mutate_with_sculpt(in_pdb, out_pdb, mutations):
    """FIX 4: PyMOL Mutagenesis Wizard + sculpt-based light minimisation
    on side-chain. Falls back to first rotamer if sculpt unavailable."""
    cmds = [
        f"load {in_pdb}, prot",
        "wizard mutagenesis",
        "refresh_wizard",
    ]
    for resi, new_aa3 in mutations:
        cmds += [
            f'cmd.get_wizard().do_select("/prot//A/{resi}/")',
            f'cmd.get_wizard().set_mode("{new_aa3}")',
            "cmd.frame(1)",  # use first rotamer
            "cmd.get_wizard().apply()",
        ]
    cmds += ["set_wizard"]

    # After all mutations applied, run sculpt to relax the mutated side chains
    # Sculpt acts on whole structure, but we only care about side-chain neighbours.
    resi_csv = "+".join(str(r) for r, _ in mutations)
    cmds += [
        f'sel_mut = "resi {resi_csv}"',
        f'sel_env = f"byres ((resi {resi_csv}) around 6)"',
        "try:",
        "    cmd.sculpt_activate('prot', state=1)",
        "    for _ in range(3):",
        "        cmd.sculpt_iterate('prot', state=1, cycles=20)",
        "    cmd.sculpt_deactivate('prot')",
        "    print('SCULPT_OK')",
        "except Exception as e:",
        "    print('SCULPT_FAIL: ' + str(e))",
        f"save {out_pdb}, prot",
    ]

    pml_path = out_pdb.replace(".pdb", ".pml")
    with open(pml_path, "w") as f:
        f.write("\n".join(cmds))
    proc = subprocess.run([PYMOL, "-cq", pml_path], capture_output=True, text=True)
    return proc


def check_clashes(pdb_path, resi_list, threshold=1.8):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("p", pdb_path)
    target_atoms, other_atoms = [], []
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
        ta_c = ta.get_coord()
        for oa in other_atoms:
            d = np.linalg.norm(ta_c - oa.get_coord())
            if d < threshold:
                clashes.append((ri, ta.get_name(), oa.get_parent().id[1], oa.get_name(), float(d)))
                break
    return clashes


def build_complex(receptor_pdb, top_pose_pdb, out_complex):
    with open(out_complex, "w") as out:
        with open(receptor_pdb) as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM", "TER")):
                    out.write(line)
        out.write("TER\n")
        with open(top_pose_pdb) as f:
            for line in f:
                if line.startswith(("ATOM  ", "HETATM")):
                    out.write(line[:21] + "X" + line[22:])
        out.write("TER\nEND\n")


def dock(rec_pdbqt, lig_pdbqt, out_pdbqt, log_file, centroid):
    cx, cy, cz = centroid
    cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", "22", "--size_y", "22", "--size_z", "22",
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
    log("Stage 7 v3 starting")

    selected = json.load(open(os.path.join(AS_DIR, "selected_meta.json")))["selected"]

    import requests
    try:
        rseq = requests.get("https://rest.uniprot.org/uniprotkb/P04818.fasta", timeout=30).text
        ref_seq = "".join(rseq.split("\n")[1:])
    except Exception as e:
        log(f"uniprot fetch failed: {e}; falling back to v2 cached seq")
        # If network fails, derive from PDB chain A
        from Bio.PDB import PDBParser
        seq = ""
        s = PDBParser(QUIET=True).get_structure("p", os.path.join(STR_DIR, "protein_dimer_h.pdb"))
        from Bio.PDB.Polypeptide import three_to_one
        last_resi = 0
        for r in s[0]["A"]:
            if r.id[0] != " ": continue
            ri = r.id[1]
            try:
                aa = three_to_one(r.get_resname())
            except Exception:
                aa = "X"
            while last_resi + 1 < ri:
                seq += "X"; last_resi += 1
            seq += aa
            last_resi = ri
        ref_seq = "X" * (1) + seq  # offset so [p-1] indexes by 1
        # Actually we need ref_seq[p-1] to give residue at position p (1-indexed)
        # Build a position-indexed sequence
        pos_seq = ["X"] * 350
        s = PDBParser(QUIET=True).get_structure("p", os.path.join(STR_DIR, "protein_dimer_h.pdb"))
        for r in s[0]["A"]:
            if r.id[0] != " ": continue
            try:
                aa = three_to_one(r.get_resname())
                if 0 < r.id[1] <= len(pos_seq):
                    pos_seq[r.id[1]-1] = aa
            except Exception:
                pass
        ref_seq = "".join(pos_seq)
    log(f"ref seq len={len(ref_seq)}")

    panel = build_panel(selected, ref_seq)
    log(f"Built panel of {len(panel)} mutants (no G217W)")

    prot_h = os.path.join(STR_DIR, "protein_dimer_h.pdb")
    cof_a = os.path.join(WT_DIR, "cofactor_chainA_h_v3.pdb")
    cof_b = os.path.join(WT_DIR, "cofactor_chainB_h_v3.pdb")
    native_lig = os.path.join(STR_DIR, "ligand.pdb")
    centroid = crystal_dump_centroid(native_lig)
    log(f"centroid (crystal dUMP): {centroid.tolist()}")

    lig_pdbqt = os.path.join(LIG_DIR, "dump.pdbqt")
    wt_apo_aff = json.load(open(os.path.join(WT_DIR, "wt_apo.json")))["top_affinity"]
    wt_holo_aff = json.load(open(os.path.join(WT_DIR, "wt_holo.json")))["top_affinity"]
    log(f"WT apo aff = {wt_apo_aff:.2f}, WT holo aff = {wt_holo_aff:.2f}")

    rows = []
    skipped = []
    for mi, m in enumerate(panel):
        mid = m["id"]
        muts = m["mutations"]
        cat = m["category"]
        log(f"[{mi+1}/{len(panel)}] {mid} ({cat})")

        sub = os.path.join(MUT_DIR, mid)
        os.makedirs(sub, exist_ok=True)

        # 1) PyMOL mutate + sculpt
        mut_pdb = os.path.join(sub, f"{mid}_mut.pdb")
        proc = pymol_mutate_with_sculpt(prot_h, mut_pdb, muts)
        if proc.returncode != 0 or not os.path.exists(mut_pdb):
            log(f"  pymol failed rc={proc.returncode}: {proc.stderr[:200]}")
            skipped.append((mid, "pymol_failed"))
            continue
        sculpt_msg = "OK" if "SCULPT_OK" in proc.stdout else (
            "FAIL" if "SCULPT_FAIL" in proc.stdout else "UNKNOWN")
        log(f"  sculpt: {sculpt_msg}")

        # 2) Clash check (informational)
        resi_list = [r for r, _ in muts]
        clashes = check_clashes(mut_pdb, resi_list, threshold=1.8)
        if clashes:
            log(f"  {len(clashes)} clashes <1.8A: {clashes[:2]}")

        # 3) Add hydrogens
        mut_h = os.path.join(sub, f"{mid}_mut_h.pdb")
        subprocess.run([OBABEL, mut_pdb, "-O", mut_h, "-h"],
                       capture_output=True, text=True)

        # 4) Receptors with Gasteiger charges
        rec_apo = os.path.join(sub, f"{mid}_apo.pdbqt")
        rec_holo = os.path.join(sub, f"{mid}_holo.pdbqt")
        ok_apo, m_apo = prepare_receptor_with_charges(mut_h, rec_apo, f"{mid}_apo")
        if not ok_apo:
            log(f"  apo prep failed; skipping")
            skipped.append((mid, "apo_prep_failed"))
            continue
        holo_pdb = os.path.join(sub, f"{mid}_holo.pdb")
        make_holo_dimer(mut_h, cof_a, cof_b, holo_pdb)
        ok_holo, m_holo = prepare_receptor_with_charges(holo_pdb, rec_holo, f"{mid}_holo")
        if not ok_holo:
            log(f"  holo prep failed; using apo")
            rec_holo = rec_apo

        # 5) Dock both conditions
        for cond, rec in [("apo", rec_apo), ("holo", rec_holo)]:
            out_pdbqt = os.path.join(sub, f"{mid}_{cond}.pdbqt")
            log_file = os.path.join(sub, f"{mid}_{cond}.log")
            proc = dock(rec, lig_pdbqt, out_pdbqt, log_file, centroid)
            if proc.returncode != 0:
                log(f"  vina {mid} {cond} rc={proc.returncode}: {proc.stderr[:200]}")
                rows.append({"mutant": mid, "category": cat, "condition": cond,
                            "top_affinity": float("nan"), "mean_topk": float("nan"),
                            "n_modes": 0,
                            "rmsd_to_native": float("nan"),
                            "delta_vina_vs_wt": float("nan"),
                            "mis_docked": True, "n_clashes": len(clashes),
                            "error": proc.stderr[:200]})
                continue

            affs = parse_vina_pdbqt(out_pdbqt)
            n_modes = len(affs)
            top = affs[0] if affs else float("nan")
            mean_topk = float(np.mean(affs[:min(3, n_modes)])) if n_modes > 0 else float("nan")

            poses = parse_pdbqt_models(out_pdbqt)
            native = native_heavy(native_lig)
            rmsd = rmsd_top(poses[0], native) if poses else float("nan")
            mis_docked = (not math.isnan(rmsd)) and rmsd > 3.0

            wt_aff = wt_apo_aff if cond == "apo" else wt_holo_aff
            # FIX 3: positive = destabilising
            delta = (top - wt_aff) if not math.isnan(top) else float("nan")

            top_pdbqt = os.path.join(sub, f"{mid}_{cond}_top.pdbqt")
            split_top(out_pdbqt, top_pdbqt)
            top_pose_pdb = os.path.join(VIEWER_DIR, f"{mid}_{cond}_top_pose.pdb")
            subprocess.run([OBABEL, top_pdbqt, "-O", top_pose_pdb],
                           capture_output=True, text=True)
            complex_pdb = os.path.join(VIEWER_DIR, f"{mid}_{cond}_complex.pdb")
            recv_pdb = mut_h if cond == "apo" else holo_pdb
            build_complex(recv_pdb, top_pose_pdb, complex_pdb)

            rows.append({"mutant": mid, "category": cat, "condition": cond,
                        "top_affinity": top, "mean_topk": mean_topk,
                        "n_modes": n_modes,
                        "rmsd_to_native": rmsd,
                        "delta_vina_vs_wt": delta,
                        "mis_docked": mis_docked, "n_clashes": len(clashes),
                        "complex_pdb": complex_pdb, "top_pose_pdb": top_pose_pdb})
            log(f"  {mid} {cond}: top={top:.2f} delta={delta:+.2f} n={n_modes} "
                f"rmsd={rmsd:.2f} mis_docked={mis_docked}")

    # CSV
    fieldnames = ["mutant", "category", "condition", "top_affinity", "mean_topk",
                 "n_modes", "rmsd_to_native", "delta_vina_vs_wt",
                 "mis_docked", "n_clashes", "complex_pdb", "top_pose_pdb", "error"]
    csv_path = os.path.join(MUT_DIR, "mutant_results_v3.csv")
    with open(csv_path, "w", newline="") as f:
        f.write("# Sign convention: delta_vina_vs_wt = top_aff_mut - top_aff_wt; "
                "positive = destabilising. mean_topk = mean(affinities[:min(3, n_modes)]).\n")
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    log(f"wrote {csv_path} ({len(rows)} rows)")

    with open(os.path.join(MUT_DIR, "summary_v3.json"), "w") as f:
        json.dump({"n_panel": len(panel), "n_completed": len(rows),
                  "skipped": skipped, "wt_apo": wt_apo_aff, "wt_holo": wt_holo_aff,
                  "sign_convention": "delta_vina_vs_wt = top_aff_mut - top_aff_wt; positive = destabilising",
                  "mean_topk_rule": "mean(affinities[:min(3, n_modes)])",
                  "mis_docked_threshold_A": 3.0,
                  "rotamer_method": "PyMOL Mutagenesis Wizard frame=1 + sculpt 3x20 cycles"},
                  f, indent=2, default=str)

    log(f"Stage 7 v3 DONE: {len(rows)} runs, {len(skipped)} skipped")


if __name__ == "__main__":
    main()
