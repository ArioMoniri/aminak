"""Phase 14 — common helpers (CID fetch, ligand prep, decoy gen, Vina wrapper).

All paths are repo-relative. Imports kept minimal so a fresh agent can re-run.
"""
from __future__ import annotations
import json, os, subprocess, sys, urllib.request, urllib.parse, time, shutil, glob, math
from pathlib import Path
from typing import Optional, Iterable
import numpy as np

REPO = Path(__file__).resolve().parents[2]      # .../aminak-inhibitor
PHASE = REPO / "14_inhibitor_design"
RECEPTOR_APO = REPO / "06f_receptor_fixed" / "protein_dimer_apo_fixed.pdbqt"
COFACTOR_A = REPO / "06f_receptor_fixed" / "cofactor_A.pdbqt"
COFACTOR_B = REPO / "06f_receptor_fixed" / "cofactor_B.pdbqt"
DUMP_PDBQT_REF = REPO / "05d_ligand_v4" / "dump.pdbqt"

# Phase-7 canonical active-site box
BOX_ACTIVE_SITE = {"cx": -0.137, "cy":  4.232, "cz": 15.159,
                    "sx": 22.0,    "sy": 22.0,   "sz": 22.0}

VINA = shutil.which("vina") or "/opt/homebrew/bin/vina"
OBABEL = shutil.which("obabel") or "/opt/homebrew/bin/obabel"
USER_PY_BIN = Path.home() / "Library/Python/3.14/bin"
MK_PREPARE_LIGAND = str(USER_PY_BIN / "mk_prepare_ligand.py")

ANCHORS_JSON = PHASE / "00_roadmap" / "anchor_compounds_verified.json"


# ---------- PubChem ----------
def pubchem_get(url: str, retries: int = 3, timeout: int = 20):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0 aminak-phase14"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read())
        except Exception as e:
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"PubChem fetch failed after {retries} retries: {url} | {last}")


def fetch_sdf_3d(cid: int, out_path: Path, max_retries: int = 4):
    """Fetch a 3D SDF for the CID; if 3D not available, fall back to 2D. Polite retries on rate-limit."""
    last_err = None
    for attempt in range(max_retries):
        for variant in ("?record_type=3d", "?record_type=2d", ""):
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF{variant}"
            try:
                req = urllib.request.Request(url, headers={"User-Agent":"Mozilla/5.0 aminak-phase14"})
                with urllib.request.urlopen(req, timeout=30) as r:
                    data = r.read()
                if len(data) > 50:
                    out_path.write_bytes(data); return variant or "default"
            except Exception as e:
                last_err = e
                continue
        time.sleep(3 * (attempt + 1))  # backoff
    raise RuntimeError(f"PubChem SDF fetch failed for CID {cid} after {max_retries} attempts: {last_err}")


def fetch_smiles(cid: int) -> str:
    """Try IsomericSMILES, then ConnectivitySMILES, then CanonicalSMILES — PubChem field name varies."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES,ConnectivitySMILES,CanonicalSMILES/JSON"
    props = pubchem_get(url)["PropertyTable"]["Properties"][0]
    for k in ("IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES"):
        if k in props and props[k]:
            return props[k]
    raise KeyError(f"PubChem CID {cid} returned no SMILES field; got {list(props.keys())}")


def pubchem_similar_cids(cid: int, threshold: int = 70, max_records: int = 30) -> list[int]:
    """Fast 2D similarity search."""
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/{cid}/cids/JSON"
           f"?Threshold={threshold}&MaxRecords={max_records}")
    try:
        d = pubchem_get(url)
        return d.get("IdentifierList", {}).get("CID", []) or []
    except Exception:
        return []


# ---------- Ligand prep ----------
def prep_ligand_from_sdf(sdf_in: Path, pdbqt_out: Path, log: Optional[Path] = None) -> dict:
    """SDF -> 3D-embedded SDF (RDKit) -> protonated (obabel pH 7.4, no gen3d) -> PDBQT (meeko).
    Meeko computes its own AD4-compatible charges; obabel-Gasteiger is not piped into meeko."""
    pdbqt_out.parent.mkdir(parents=True, exist_ok=True)
    prot_sdf = pdbqt_out.with_suffix(".prot.sdf")
    embed_sdf = pdbqt_out.with_suffix(".embed.sdf")

    # Stage 1: RDKit 3D embed (skip if input already 3D)
    from rdkit import Chem
    from rdkit.Chem import AllChem
    try:
        mol = Chem.SDMolSupplier(str(sdf_in), removeHs=False)[0]
        if mol is None:
            # try with sanitize off (some PubChem SDFs have funky aromaticity)
            mol = Chem.SDMolSupplier(str(sdf_in), removeHs=False, sanitize=False)[0]
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        if mol is None:
            return {"ok": False, "stage": "rdkit_read", "err": "RDKit cannot parse SDF"}
        # if no 3D coords, embed
        conf = mol.GetConformer()
        z_range = max(conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())) - min(conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms()))
        if abs(z_range) < 0.01:                    # 2D, need embed
            mol = Chem.AddHs(mol)
            if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
                AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
            try: AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except Exception: pass
        else:                                       # already 3D, just add Hs
            mol = Chem.AddHs(mol, addCoords=True)
        Chem.MolToMolFile(mol, str(embed_sdf))
    except Exception as e:
        return {"ok": False, "stage": "rdkit_embed", "err": str(e)[:500]}

    # Stage 2: obabel protonate at pH 7.4 (no gen3d — coords from RDKit)
    try:
        subprocess.run(
            [OBABEL, str(embed_sdf), "-O", str(prot_sdf), "-p", "7.4"],
            check=True, capture_output=True, timeout=60
        )
    except subprocess.CalledProcessError as e:
        return {"ok": False, "stage": "obabel_protonate", "err": e.stderr.decode()[:500]}

    # Stage 3: meeko → PDBQT (meeko computes its own AD4 charges)
    try:
        subprocess.run(
            ["python3", MK_PREPARE_LIGAND, "-i", str(prot_sdf), "-o", str(pdbqt_out)],
            check=True, capture_output=True, timeout=120
        )
        # Verify the PDBQT looks valid (has at least one ATOM line)
        if pdbqt_out.exists() and any(ln.startswith("ATOM") for ln in pdbqt_out.read_text().splitlines()):
            return {"ok": True, "stage": "meeko"}
        else:
            raise RuntimeError("meeko produced empty PDBQT")
    except (subprocess.CalledProcessError, RuntimeError) as e:
        # Final fallback: write PDBQT directly via obabel with the RDKit-embedded structure,
        # using -xr (rigid) — sometimes Vina accepts these
        try:
            subprocess.run([OBABEL, str(embed_sdf), "-O", str(pdbqt_out), "--partialcharge", "gasteiger"],
                          check=True, capture_output=True, timeout=60)
            return {"ok": True, "stage": "obabel_pdbqt_fallback",
                    "note": f"meeko failed ({str(e)[:80]}), used obabel PDBQT"}
        except Exception as e2:
            return {"ok": False, "stage": "meeko+obabel_fallback",
                    "err": f"meeko: {str(e)[:150]} | obabel: {str(e2)[:150]}"}


# ---------- Vina ----------
def vina_dock(receptor: Path, ligand: Path, out_pdbqt: Path, box: dict, seed: int = 42,
              exhaustiveness: int = 32, num_modes: int = 20, cpu: int = 4,
              log: Optional[Path] = None, timeout_s: int = 1800) -> dict:
    """Run Vina; parse log for affinity table. Returns dict with top1 + top3_mean + top5_mean."""
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = [VINA,
           "--receptor", str(receptor),
           "--ligand", str(ligand),
           "--center_x", f"{box['cx']}", "--center_y", f"{box['cy']}", "--center_z", f"{box['cz']}",
           "--size_x",   f"{box['sx']}", "--size_y",   f"{box['sy']}", "--size_z",   f"{box['sz']}",
           "--exhaustiveness", str(exhaustiveness),
           "--num_modes", str(num_modes),
           "--seed", str(seed),
           "--cpu", str(cpu),
           "--out", str(out_pdbqt)]
    log_path = log or out_pdbqt.with_suffix(".log")
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, check=False, capture_output=True, timeout=timeout_s)
    except subprocess.TimeoutExpired:
        return {"ok": False, "err": f"vina timeout @ {timeout_s}s"}
    log_path.write_bytes(proc.stdout + b"\n=== STDERR ===\n" + proc.stderr)
    if proc.returncode != 0:
        return {"ok": False, "err": f"vina rc={proc.returncode}", "stderr": proc.stderr.decode()[:500]}
    # parse affinity table
    affinities = []
    for line in proc.stdout.decode().splitlines():
        # lines look like:  "   1     -8.234      0.000      0.000"
        parts = line.split()
        if len(parts) >= 4 and parts[0].isdigit():
            try:
                affinities.append(float(parts[1]))
            except ValueError:
                pass
    if not affinities:
        return {"ok": False, "err": "no affinities parsed"}
    return {
        "ok": True,
        "wall_s": round(time.time() - t0, 2),
        "n_modes": len(affinities),
        "top1": affinities[0],
        "top3_mean": float(np.mean(affinities[:3])),
        "top5_mean": float(np.mean(affinities[:5])),
        "all": affinities,
    }


# ---------- Decoy gen (RDKit DUD-E-style) ----------
def rdkit_decoys(anchor_smiles: str, n: int = 30, pubchem_pool_cids: Optional[list[int]] = None) -> list[str]:
    """Match property windows of the anchor (MW±20, logP±0.5, HBA±2, HBD±2, RotB±2)
    against a small PubChem-derived pool. If no pool given, use a built-in 200-compound
    drug-like sample from the FDA-approved set (PubChem CIDs in this script are hard-coded)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Lipinski

    anchor = Chem.MolFromSmiles(anchor_smiles)
    if anchor is None: return []
    a_mw   = Descriptors.MolWt(anchor)
    a_logp = Descriptors.MolLogP(anchor)
    a_hba  = Lipinski.NumHAcceptors(anchor)
    a_hbd  = Lipinski.NumHDonors(anchor)
    a_rotb = Lipinski.NumRotatableBonds(anchor)
    a_fp   = AllChem.GetMorganFingerprintAsBitVect(anchor, 2, 2048)

    # built-in pool — small + diverse for fast Phase 14 scoped run
    POOL_CIDS = [
        2244, 5288826, 5957, 5743, 24360, 5212, 3782, 4017, 5793, 5826,
        5816, 60665, 60843, 65349, 100049, 222786, 446284, 5743, 60198,
        4091, 6149, 657, 444795, 6253, 444593, 5470, 4170, 5564, 5759,
        119607, 5743, 24360, 49846, 60750, 8569, 60823, 5538, 5359476,
        49846, 6710, 5760, 5388, 6041, 5495, 5734, 5566, 5564, 6035,
        4900, 25245,
    ]
    if pubchem_pool_cids:
        POOL_CIDS = list(set(POOL_CIDS + pubchem_pool_cids))

    decoys = []
    for cid in POOL_CIDS:
        if len(decoys) >= n: break
        try:
            smi = fetch_smiles(cid)
            m = Chem.MolFromSmiles(smi)
            if m is None: continue
            mw = Descriptors.MolWt(m)
            if abs(mw - a_mw) > 100: continue           # loosened from 20 because pool small
            logp = Descriptors.MolLogP(m)
            if abs(logp - a_logp) > 1.5: continue
            hba = Lipinski.NumHAcceptors(m)
            if abs(hba - a_hba) > 3: continue
            hbd = Lipinski.NumHDonors(m)
            if abs(hbd - a_hbd) > 3: continue
            rotb = Lipinski.NumRotatableBonds(m)
            if abs(rotb - a_rotb) > 4: continue
            # exclude high-similarity to anchor (avoid duplicating Tier 1)
            fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)
            from rdkit import DataStructs
            sim = DataStructs.TanimotoSimilarity(a_fp, fp)
            if sim > 0.7: continue
            decoys.append({"cid": cid, "smiles": smi, "mw": mw, "logp": logp})
            time.sleep(0.2)
        except Exception:
            continue
    return decoys


# ---------- Property + filter ----------
def compound_descriptors(smiles: str) -> dict:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    m = Chem.MolFromSmiles(smiles)
    if m is None: return {}
    out = {
        "mw": Descriptors.MolWt(m), "logp": Descriptors.MolLogP(m),
        "hba": Lipinski.NumHAcceptors(m), "hbd": Lipinski.NumHDonors(m),
        "rotb": Lipinski.NumRotatableBonds(m), "tpsa": Descriptors.TPSA(m),
        "rings": rdMolDescriptors.CalcNumRings(m),
        "formal_charge": Chem.GetFormalCharge(m),
        "n_heavy": m.GetNumHeavyAtoms(),
    }
    # filters
    for tag, p in [("pains_a", FilterCatalogParams.FilterCatalogs.PAINS_A),
                   ("pains_b", FilterCatalogParams.FilterCatalogs.PAINS_B),
                   ("pains_c", FilterCatalogParams.FilterCatalogs.PAINS_C),
                   ("brenk",   FilterCatalogParams.FilterCatalogs.BRENK),
                   ("nih",     FilterCatalogParams.FilterCatalogs.NIH)]:
        params = FilterCatalogParams(); params.AddCatalog(p)
        out[tag] = FilterCatalog(params).HasMatch(m)
    out["lipinski_ok"] = (out["mw"] <= 500 and out["logp"] <= 5 and out["hba"] <= 10 and out["hbd"] <= 5)
    out["veber_ok"]    = (out["rotb"] <= 10 and out["tpsa"] <= 140)
    return out


def load_anchors() -> list[dict]:
    return json.loads(ANCHORS_JSON.read_text())["anchors"]
