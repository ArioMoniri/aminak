"""Step 6 — Validation: local Ramachandran plots, per-residue DOPE profiles,
and a single quality_overview.png comparing all models.

Also writes SAVES_MANUAL.md describing how to upload to UCLA SAVES for the full PROCHECK suite.
"""
from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    MODELS_DIR, STEP3_DIR, STEP4_DIR, STEP6_DIR, png_ok, setup_logger,
)

LOG = setup_logger("step6_validate")

# Approximate Ramachandran allowed regions for general residues.
# Coordinates in degrees (phi, psi). Polygons used for "favoured" + "allowed".
GENERAL_FAVOURED = [
    # Right-handed alpha helix
    [(-100, -70), (-50, -70), (-40, -60), (-50, -25), (-100, -25)],
    # Beta sheet
    [(-180, 90), (-180, 180), (-50, 180), (-50, 90)],
    [(-180, -180), (-180, -120), (-50, -120), (-50, -180)],
    # Left-handed alpha
    [(40, 30), (90, 30), (90, 80), (40, 80)],
]
GENERAL_ALLOWED = [
    [(-180, -180), (-180, 180), (-30, 180), (-30, -180)],
    [(30, 0), (130, 0), (130, 110), (30, 110)],
]


def in_polygon(x: float, y: float, poly: list[tuple[float, float]]) -> bool:
    """Ray-casting point-in-polygon."""
    n = len(poly)
    inside = False
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi):
            inside = not inside
        j = i
    return inside


def classify(phi: float, psi: float) -> str:
    """Return 'favoured', 'allowed', or 'outlier'."""
    for poly in GENERAL_FAVOURED:
        if in_polygon(phi, psi, poly):
            return "favoured"
    for poly in GENERAL_ALLOWED:
        if in_polygon(phi, psi, poly):
            return "allowed"
    return "outlier"


def compute_phi_psi(pdb_path: Path) -> list[tuple[int, float, float, str]]:
    """Return list of (residue_number, phi_deg, psi_deg, residue_name)."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    ppb = PPBuilder()
    out: list[tuple[int, float, float, str]] = []
    for model in structure:
        for chain in model:
            for pp in ppb.build_peptides(chain):
                phi_psi_list = pp.get_phi_psi_list()
                for res, (phi, psi) in zip(pp, phi_psi_list):
                    if phi is None or psi is None:
                        continue
                    out.append((res.id[1],
                                math.degrees(phi),
                                math.degrees(psi),
                                res.get_resname()))
    return out


def render_ramachandran(phi_psi: list[tuple[int, float, float, str]],
                        png_out: Path, label: str) -> tuple[float, float, float]:
    """Render Ramachandran plot. Returns (favoured%, allowed%, outlier%)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection

    fig, ax = plt.subplots(figsize=(6.5, 6))
    # Favoured / allowed regions
    fav_patches = [Polygon(poly, closed=True) for poly in GENERAL_FAVOURED]
    all_patches = [Polygon(poly, closed=True) for poly in GENERAL_ALLOWED]
    ax.add_collection(PatchCollection(all_patches, facecolor="#dde8f4",
                                       edgecolor="#7090b8", linewidths=0.7,
                                       alpha=0.9, zorder=1))
    ax.add_collection(PatchCollection(fav_patches, facecolor="#9cc6ec",
                                       edgecolor="#3a6ea8", linewidths=0.8,
                                       alpha=0.95, zorder=2))

    n = len(phi_psi)
    n_fav = n_all = n_out = 0
    glycines = []
    prolines = []
    others = []
    for resn, phi, psi, name in phi_psi:
        cls = classify(phi, psi)
        if cls == "favoured":
            n_fav += 1
        elif cls == "allowed":
            n_all += 1
        else:
            n_out += 1
        if name == "GLY":
            glycines.append((phi, psi))
        elif name == "PRO":
            prolines.append((phi, psi))
        else:
            others.append((phi, psi))

    if others:
        xs, ys = zip(*others)
        ax.scatter(xs, ys, s=14, color="#222222", alpha=0.85, label="general", zorder=3)
    if glycines:
        xs, ys = zip(*glycines)
        ax.scatter(xs, ys, s=22, color="#cc3333", marker="^", alpha=0.9,
                   label="GLY", zorder=4)
    if prolines:
        xs, ys = zip(*prolines)
        ax.scatter(xs, ys, s=22, color="#229922", marker="s", alpha=0.9,
                   label="PRO", zorder=4)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-180, 181, 60))
    ax.axhline(0, color="grey", linewidth=0.4)
    ax.axvline(0, color="grey", linewidth=0.4)
    ax.set_xlabel("φ (degrees)")
    ax.set_ylabel("ψ (degrees)")
    pct_fav = 100.0 * n_fav / max(n, 1)
    pct_all = 100.0 * n_all / max(n, 1)
    pct_out = 100.0 * n_out / max(n, 1)
    ax.set_title(f"Ramachandran: {label}\n"
                 f"favoured {pct_fav:.1f}%   allowed {pct_all:.1f}%   "
                 f"outlier {pct_out:.1f}%   (n={n})")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(linestyle="--", linewidth=0.3, alpha=0.4)
    fig.tight_layout()
    fig.savefig(png_out, dpi=180)
    plt.close(fig)
    return pct_fav, pct_all, pct_out


def per_residue_dope(model_pdb: Path, alignment_path: Path,
                     sequence_code: str, knowns: list[str]) -> list[float] | None:
    """Use Modeller's complete_pdb + assess_dope to get per-residue DOPE.

    Returns list of per-residue DOPE energies, or None on failure.
    """
    cwd_before = os.getcwd()
    os.chdir(model_pdb.parent)
    try:
        from modeller import Environ, Selection, log
        from modeller.scripts import complete_pdb

        log.none()
        env = Environ()
        env.libs.topology.read(file="$(LIB)/top_heav.lib")
        env.libs.parameters.read(file="$(LIB)/par.lib")
        mdl = complete_pdb(env, str(model_pdb))
        sel = Selection(mdl)
        ene_path = model_pdb.with_suffix(".profile")
        sel.assess_dope(output="ENERGY_PROFILE NO_REPORT",
                        file=str(ene_path), normalize_profile=True,
                        smoothing_window=15)
        per_res: list[float] = []
        if ene_path.exists():
            for line in ene_path.read_text().splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                # Skip lines with too few columns
                if len(parts) < 3:
                    continue
                try:
                    val = float(parts[-1])
                    per_res.append(val)
                except (ValueError, IndexError):
                    continue
        return per_res
    except Exception as e:
        LOG.warning("per_residue_dope failed for %s: %s", model_pdb.name, e)
        return None
    finally:
        os.chdir(cwd_before)


def render_dope_profile(per_res: list[float], png_out: Path, label: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 4))
    x = list(range(1, len(per_res) + 1))
    ax.plot(x, per_res, color="#1f5e9c", linewidth=1.0)
    ax.fill_between(x, per_res, 0, where=[v < 0 for v in per_res],
                    color="#7fbfe2", alpha=0.5)
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.set_xlabel("Residue index (alignment position)")
    ax.set_ylabel("Normalized DOPE per residue")
    ax.set_title(f"Per-residue DOPE profile: {label}")
    ax.grid(linestyle="--", linewidth=0.4, alpha=0.6)
    fig.tight_layout()
    fig.savefig(png_out, dpi=180)
    plt.close(fig)


def render_quality_overview(scores_csv: Path, rmsd_csv: Path,
                            rama_stats: list[dict], png_out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Load scores
    scores_by_id: dict[int, dict] = {}
    with open(scores_csv) as fh:
        for row in csv.DictReader(fh):
            scores_by_id[int(row["model_id"])] = row
    rmsd_by_id: dict[int, dict] = {}
    with open(rmsd_csv) as fh:
        for row in csv.DictReader(fh):
            rmsd_by_id[int(row["model_id"])] = row
    rama_by_id: dict[int, dict] = {r["model_id"]: r for r in rama_stats}

    ids = sorted(scores_by_id.keys())
    dope_vals = [float(scores_by_id[i]["DOPE"]) for i in ids]
    molpdf = [float(scores_by_id[i]["molpdf"]) for i in ids]
    rmsds = [float(rmsd_by_id.get(i, {}).get("rmsd_to_crystal", "nan")) for i in ids]
    fav = [float(rama_by_id.get(i, {}).get("pct_favoured", 0.0)) for i in ids]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    ax = axes[0, 0]
    ax.bar(ids, dope_vals, color="#3a6ea8")
    ax.set_title("DOPE score per model (lower = better)")
    ax.set_xlabel("Model")
    ax.set_ylabel("DOPE")
    ax.grid(axis="y", linestyle="--", linewidth=0.3, alpha=0.5)

    ax = axes[0, 1]
    ax.bar(ids, molpdf, color="#a85a3a")
    ax.set_title("molpdf per model (lower = better)")
    ax.set_xlabel("Model")
    ax.set_ylabel("molpdf")
    ax.grid(axis="y", linestyle="--", linewidth=0.3, alpha=0.5)

    ax = axes[1, 0]
    ax.bar(ids, rmsds, color="#3a8a5e")
    ax.set_title("Cα RMSD to 1HVY chain A (Å, lower = closer)")
    ax.set_xlabel("Model")
    ax.set_ylabel("RMSD (Å)")
    ax.grid(axis="y", linestyle="--", linewidth=0.3, alpha=0.5)

    ax = axes[1, 1]
    ax.bar(ids, fav, color="#8a6e3a")
    ax.axhline(90, color="grey", linestyle="--", linewidth=0.7)
    ax.set_title("Ramachandran favoured % (target ≥ 90%)")
    ax.set_xlabel("Model")
    ax.set_ylabel("Favoured (%)")
    ax.set_ylim(0, 100)
    ax.grid(axis="y", linestyle="--", linewidth=0.3, alpha=0.5)

    fig.suptitle("Phase 6 Modeller homology models — quality overview",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(png_out, dpi=170)
    plt.close(fig)


def main() -> int:
    STEP6_DIR.mkdir(parents=True, exist_ok=True)
    knowns_meta = json.loads((STEP3_DIR / "knowns.json").read_text())
    knowns: list[str] = knowns_meta["knowns"]
    sequence_code: str = knowns_meta["sequence_code"]
    alignment_path = STEP3_DIR / "alignment.ali"

    models = sorted(MODELS_DIR.glob("target.B99990*.pdb"))
    LOG.info("Validating %d models", len(models))

    rama_stats_rows: list[dict] = []
    for i, m in enumerate(models, start=1):
        # Ramachandran
        try:
            phi_psi = compute_phi_psi(m)
        except Exception as e:
            LOG.warning("phi/psi failed for %s: %s", m.name, e)
            phi_psi = []
        png = STEP6_DIR / f"ramachandran_model{i:02d}.png"
        pct_fav, pct_all, pct_out = render_ramachandran(phi_psi, png, m.stem)
        rama_stats_rows.append({
            "model_id": i,
            "model_pdb": m.name,
            "n_residues_with_phi_psi": len(phi_psi),
            "pct_favoured": round(pct_fav, 2),
            "pct_allowed": round(pct_all, 2),
            "pct_outlier": round(pct_out, 2),
        })
        LOG.info("Rama model %02d: favoured=%.1f%% allowed=%.1f%% outlier=%.1f%%",
                 i, pct_fav, pct_all, pct_out)

        # Per-residue DOPE
        per_res = per_residue_dope(m, alignment_path, sequence_code, knowns)
        if per_res:
            png_dope = STEP6_DIR / f"dope_profile_model{i:02d}.png"
            try:
                render_dope_profile(per_res, png_dope, m.stem)
            except Exception as e:
                LOG.warning("DOPE profile render failed for %s: %s", m.name, e)
        else:
            LOG.warning("No per-residue DOPE for %s", m.name)

    # Save Ramachandran stats CSV
    csv_path = STEP6_DIR / "ramachandran_stats.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["model_id", "model_pdb",
                                            "n_residues_with_phi_psi",
                                            "pct_favoured", "pct_allowed",
                                            "pct_outlier"])
        w.writeheader()
        w.writerows(rama_stats_rows)
    LOG.info("Wrote %s", csv_path)

    # Quality overview
    overview = STEP6_DIR / "quality_overview.png"
    try:
        render_quality_overview(STEP4_DIR / "scores.csv",
                                Path(__file__).parent.parent.parent
                                / "10_modeller" / "05_comparison" / "rmsd_per_model.csv",
                                rama_stats_rows, overview)
    except Exception as e:
        LOG.warning("Quality overview render failed: %s", e)
    if not png_ok(overview):
        LOG.warning("quality_overview.png too small or missing")

    # SAVES_MANUAL.md
    saves_md = STEP6_DIR / "SAVES_MANUAL.md"
    saves_md.write_text(
        "# Manual UCLA SAVES validation\n\n"
        "For full PROCHECK / ERRAT / VERIFY3D / WHATCHECK validation, upload\n"
        "`10_modeller/04_modeller_run/models/best_model.pdb` to "
        "<https://saves.mbi.ucla.edu/>.\n\n"
        "The SAVES web tool has no programmatic API; the upload must be done\n"
        "manually OR via the project's computer-use channel.\n\n"
        "Local Ramachandran (this folder) is the open-source equivalent of\n"
        "PROCHECK's basic Ramachandran plot — it uses Biopython's φ/ψ\n"
        "computation and a hand-drawn favoured/allowed-region overlay.\n\n"
        "Steps (manual):\n"
        "1. Visit https://saves.mbi.ucla.edu/.\n"
        "2. Upload `best_model.pdb` (≈ a few hundred KB).\n"
        "3. Choose 'Run all programs' and submit.\n"
        "4. Download the PROCHECK Ramachandran PDF and ERRAT/VERIFY3D HTMLs;\n"
        "   place into this folder.\n"
    )
    LOG.info("Wrote %s", saves_md)
    LOG.info("STEP 6 OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
