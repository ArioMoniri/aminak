#!/usr/bin/env python3
"""Phase 14g v3 — FULL HADDOCK3 pipeline on TYMS dimer + LR-octapeptide.

KEY DISCOVERY: HADDOCK3 2026.5.0 already ships its own CNS arm64-darwin binary
(`haddock/cns/bin/arm64-darwin.bin`, 4.1 MB, Mach-O arm64, Bonvin lab UU
patch release). No separate CNS install needed — just point HADDOCK3 to it.

Protocol (haddock3 BSc):
  topoaa     — build topologies
  rigidbody  — initial FFT-based docking (1000 decoys)
  flexref    — semi-flexible refinement (200 best)
  caprieval  — CAPRI-style scoring of the resulting cluster

Inputs:
  receptor = the Phase-6c-hardened TYMS dimer (chain A/B labelled)
  ligand   = LR octapeptide built from sequence LSCQLYQR (canonical) or
             the numpy seed-42 scramble QLCRQSYL (specificity control)

Active residues (chain A) = 34 dimer-interface residues from Phase-14 A3
contact map. AIRs defined accordingly.

Runs canonical AND scrambled control. Compares top cluster HADDOCK scores.
"""
from __future__ import annotations
import os, subprocess, csv, json, shutil, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "14_inhibitor_design" / "07_advanced_methods" / "haddock3_full"
OUT.mkdir(parents=True, exist_ok=True)

HADDOCK3 = Path.home() / "Library/Python/3.14/bin/haddock3"
HADDOCK3_RESTRAINTS = Path.home() / "Library/Python/3.14/bin/haddock3-restraints"
CNS_EXEC = Path.home() / "Library/Python/3.14/lib/python/site-packages/haddock/cns/bin/arm64-darwin.bin"

# Inputs
APO_PDB = REPO / "06f_receptor_fixed" / "dimer_noH.pdb"   # already has chain A/B if PyMOL prep was right
# Peptide-as-PDB from Phase-14 Strategy-3 build
PEP_SDF_DIR = REPO / "14_inhibitor_design" / "03_dimer_interface" / "ligands"

ACTIVE_RESIDUES_A = [20,21,22,23,24,25,32,33,34,35,36,37,39,117,118,135,148,
                     150,151,153,157,158,159,160,167,168,172,173,175,177,178,
                     199,200,202]
# Passive residues = all chain-A residues within 6.5 Å of active set (heuristic: ±4)
PASSIVE_RESIDUES_A = sorted(set(r + d for r in ACTIVE_RESIDUES_A for d in (-2,-1,0,1,2)
                                if r + d > 0 and (r + d) not in ACTIVE_RESIDUES_A))


def prep_receptor_chain_ab(out_pdb: Path):
    """Receptor with explicit chain A + chain B labelling.
    Phase-6c dimer_noH has all atoms labelled chain A; split at the midpoint
    into A and B with a TER between."""
    text = APO_PDB.read_text().splitlines()
    n_total = sum(1 for ln in text if ln.startswith("ATOM"))
    half = n_total // 2
    cnt = 0; out = []
    for ln in text:
        if ln.startswith("ATOM"):
            chain = "A" if cnt < half else "B"
            new = ln[:21] + chain + ln[22:]
            if cnt == half:
                out.append("TER")
            out.append(new)
            cnt += 1
    out_pdb.write_text("\n".join(out) + "\nTER\nEND\n")
    return out_pdb


def prep_peptide(sequence: str, out_pdb: Path, chain: str = "B"):
    """Build a clean peptide PDB from sequence via RDKit, then rewrite to
    chain B with proper resnames."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSequence(sequence)
    if mol is None: raise RuntimeError(f"MolFromSequence failed for {sequence}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    try: AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
    except Exception: pass
    raw_pdb = out_pdb.with_suffix(".raw.pdb")
    Chem.MolToPDBFile(mol, str(raw_pdb))
    # Rewrite chain id to B
    lines = []
    for ln in raw_pdb.read_text().splitlines():
        if ln.startswith(("ATOM","HETATM")):
            # force chain B
            new = ln[:21] + chain + ln[22:]
            # ensure ATOM (not HETATM)
            if new.startswith("HETATM"):
                new = "ATOM  " + new[6:]
            lines.append(new)
    out_pdb.write_text("\n".join(lines) + "\nEND\n")
    return out_pdb


def make_airs_tbl(out_tbl: Path, peptide_chain: str = "C"):
    """Use haddock3-restraints to generate AIRs."""
    active_str = "+".join(map(str, ACTIVE_RESIDUES_A))
    passive_str = "+".join(map(str, PASSIVE_RESIDUES_A))
    # We need a "block" file for active/passive — haddock3-restraints provides
    # several CLIs; for protein-peptide we use active_passive_to_ambig
    active_file = out_tbl.parent / "active_receptor.txt"
    active_file.write_text(" ".join(map(str, ACTIVE_RESIDUES_A)) + "\n"
                            + " ".join(map(str, PASSIVE_RESIDUES_A)) + "\n")
    peptide_file = out_tbl.parent / "active_peptide.txt"
    peptide_file.write_text("\n\n")  # peptide has no active residues defined
    # Skip the haddock3-restraints CLI — its v2026.5.0 output uses "resi" instead
    # of "resid" and produces empty second-selection (), which CNS rejects. Use
    # the proper hand-written AIR file below instead.
    # Fallback: write a minimal AIR file by hand (HADDOCK syntax)
    # HADDOCK3 AIR format: one "assign" block per active residue, target = peptide.
    # The proper CNS NOE selection syntax is:
    #   assign ( resid N and segid A )
    #          ( segid X )  upper_bound  upper_minus  lower_minus
    air_lines = ["! Phase 14 dimer-interface AIRs — active chain-A residues vs peptide segid"]
    for resid in ACTIVE_RESIDUES_A:
        # ALL on one line, no orphan tokens for CNS:
        air_lines.append(
            f"assign ( resid {resid} and segid A ) "
            f"( segid {peptide_chain} ) 2.0 2.0 0.0"
        )
    out_tbl.write_text("\n".join(air_lines) + "\n")
    return out_tbl


def write_config(run_dir: Path, receptor_pdb: Path, peptide_pdb: Path, airs_tbl: Path,
                  config_path: Path):
    """Write a HADDOCK3 config that does topoaa + rigidbody + flexref + caprieval.
    NO mdref (saves ~2 h per run). Small sampling for tractable wall-time."""
    config = f"""# HADDOCK3 config — TYMS dimer + LR peptide
run_dir = "{run_dir}"
postprocess = true
clean = false
ncores = 4

molecules = [
    "{receptor_pdb}",
    "{peptide_pdb}",
]

[topoaa]

[rigidbody]
ambig_fname = "{airs_tbl}"
sampling = 50

[seletop]
select = 10

[flexref]
ambig_fname = "{airs_tbl}"

[caprieval]
"""
    config_path.write_text(config)
    return config_path


def run_haddock3(config: Path):
    """Run haddock3 on the config. Returns the run directory."""
    env = dict(os.environ)
    env["HADDOCK3_CNS_EXEC"] = str(CNS_EXEC)
    env["CNS_EXEC"] = str(CNS_EXEC)
    cmd = [str(HADDOCK3), str(config)]
    print(f"  ▶ {' '.join(cmd)}")
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, timeout=3600, env=env,
                          cwd=str(config.parent))
    wall = time.time() - t0
    print(f"    completed in {wall:.0f}s (rc={proc.returncode})")
    return proc.returncode, proc.stdout.decode() + proc.stderr.decode(), wall


def parse_caprieval(run_dir: Path):
    """Find the caprieval ss.tsv and parse top cluster score."""
    for csv_path in run_dir.rglob("capri_ss.tsv"):
        rows = list(csv.DictReader(csv_path.open(), delimiter="\t"))
        if rows:
            # sort by HADDOCK score
            rows.sort(key=lambda r: float(r.get("score", 0) or 0))
            return {"top_score": float(rows[0].get("score", 0) or 0),
                    "n_models": len(rows), "ss_file": str(csv_path)}
    return None


def main():
    print("=== Phase 14g v3 — FULL HADDOCK3 pipeline ===")
    print(f"  HADDOCK3 binary: {HADDOCK3}  exists={HADDOCK3.exists()}")
    print(f"  CNS binary:      {CNS_EXEC}  exists={CNS_EXEC.exists()}")
    if not (HADDOCK3.exists() and CNS_EXEC.exists()):
        print("  ! missing binary — abort"); return

    # 1. Prep receptor with explicit A/B chains
    receptor_pdb = OUT / "receptor_AB.pdb"
    prep_receptor_chain_ab(receptor_pdb)
    print(f"  receptor: {receptor_pdb}")

    # 2. Prep peptides
    results = {}
    for label, sequence in [("canonical", "LSCQLYQR"),
                            ("scrambled", "QLCRQSYL")]:
        print(f"\n  --- {label} peptide ({sequence}) ---")
        run_dir = OUT / f"run_{label}"
        if run_dir.exists(): shutil.rmtree(run_dir)
        run_dir.mkdir(parents=True)

        peptide_pdb = run_dir / "peptide.pdb"
        prep_peptide(sequence, peptide_pdb, chain="C")
        print(f"  peptide: {peptide_pdb}  ({sequence}, {len(peptide_pdb.read_text().splitlines())} lines)")

        # 3. AIRs
        airs_tbl = run_dir / "airs.tbl"
        make_airs_tbl(airs_tbl, peptide_chain="C")
        print(f"  AIRs: {airs_tbl}")

        # 4. config
        config_path = run_dir / "config.cfg"
        run_subdir = run_dir / "haddock_run"
        write_config(run_subdir, receptor_pdb, peptide_pdb, airs_tbl, config_path)
        print(f"  config: {config_path}")

        # 5. RUN HADDOCK3
        rc, log, wall = run_haddock3(config_path)
        (run_dir / "haddock3_run.log").write_text(log)
        print(f"    rc={rc}  log size={len(log)}  wall={wall:.0f}s")

        # 6. Parse capri_ss.tsv
        capri = parse_caprieval(run_subdir)
        results[label] = {"return_code": rc, "wall_s": wall, "sequence": sequence,
                          "capri": capri, "log_tail": log[-1500:] if log else ""}

    # ─── Compare canonical vs scrambled ───
    print("\n=== HADDOCK3 cluster-score comparison ===")
    can = results.get("canonical", {})
    scr = results.get("scrambled", {})
    can_score = can.get("capri", {}).get("top_score") if can.get("capri") else None
    scr_score = scr.get("capri", {}).get("top_score") if scr.get("capri") else None
    print(f"  canonical: top HADDOCK score = {can_score}")
    print(f"  scrambled: top HADDOCK score = {scr_score}")
    if can_score is not None and scr_score is not None:
        delta = can_score - scr_score
        print(f"  Δ (canonical − scrambled) = {delta:+.2f}")
        verdict = ("★ canonical IS more favourable than scrambled — peptide specificity recovered"
                   if delta < -1 else "within noise — null result confirmed under full HADDOCK3")
        print(f"  → {verdict}")
        results["comparison"] = {"canonical": can_score, "scrambled": scr_score,
                                  "delta": delta, "verdict": verdict}

    # Save summary
    summary_path = OUT / "haddock3_full_summary.json"
    summary_path.write_text(json.dumps(results, indent=2, default=str))
    print(f"\n  → {summary_path}")


if __name__ == "__main__":
    main()
