"""Step 5 — PyMOL pairwise overlays and all-together overlay + RMSD CSV.

Calls headless pymol via subprocess (NEVER `import pymol` in this process).
Per model i: align to crystal, capture pairwise PNG, per-residue Cα RMSD bar plot, RMSD value.
All-together: load all 10 models + crystal, rainbow colour, align to crystal, render.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    MODELS_DIR, PYMOL_BIN, STEP1_DIR, STEP4_DIR, STEP5_DIR, png_ok, setup_logger,
)

LOG = setup_logger("step5_pymol_compare")


def list_models() -> list[Path]:
    return sorted(p for p in MODELS_DIR.glob("target.B*.pdb"))


def run_pymol_py(script_text: str, label: str) -> tuple[bool, str]:
    """Write a Python script that uses `from pymol import cmd` and run with pymol -cqr."""
    script_path = STEP5_DIR / f"_tmp_{label}.py"
    script_path.write_text(script_text)
    cmd = [PYMOL_BIN, "-cq", "-r", str(script_path)]
    LOG.info("PyMOL [%s]", label)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        LOG.warning("PyMOL %s rc=%d stderr=%s stdout-tail=%s",
                    label, proc.returncode, proc.stderr[:300], proc.stdout[-300:])
    return proc.returncode == 0, proc.stdout


def pairwise_script(model_pdb: Path, crystal_pdb: Path, png_out: Path,
                    rmsd_json_out: Path) -> str:
    return f"""from pymol import cmd
import json
cmd.load(r"{crystal_pdb}", "crystal")
cmd.load(r"{model_pdb}", "mdl")
cmd.hide("everything")
cmd.show("cartoon")
cmd.color("magenta", "crystal")
cmd.color("green", "mdl")
result = cmd.align("mdl and name CA", "crystal and name CA")
rmsd = float(result[0])
n_atoms = int(result[1])
with open(r"{rmsd_json_out}", "w") as fh:
    json.dump({{"rmsd": rmsd, "n_atoms_aligned": n_atoms}}, fh)
cmd.bg_color("white")
cmd.set("ray_opaque_background", "on")
cmd.orient()
cmd.zoom("crystal", buffer=2)
cmd.ray(1200, 900)
cmd.png(r"{png_out}", dpi=150)
cmd.delete("all")
"""


def heatmap_dump_script(model_pdb: Path, crystal_pdb: Path, csv_out: Path) -> str:
    return f"""from pymol import cmd
import csv, math
cmd.load(r"{crystal_pdb}", "crystal")
cmd.load(r"{model_pdb}", "mdl")
cmd.align("mdl and name CA", "crystal and name CA")
crystal_ca = {{}}
mdl_ca = {{}}
cmd.iterate_state(1, "crystal and name CA", "crystal_ca[int(resi)] = (x, y, z)",
                   space={{"crystal_ca": crystal_ca, "int": int}})
cmd.iterate_state(1, "mdl and name CA", "mdl_ca[int(resi)] = (x, y, z)",
                   space={{"mdl_ca": mdl_ca, "int": int}})
rows = []
for r in sorted(set(crystal_ca) & set(mdl_ca)):
    a = crystal_ca[r]
    b = mdl_ca[r]
    d = math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)
    rows.append((r, d))
with open(r"{csv_out}", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["resi", "ca_dist"])
    for r, d in rows:
        w.writerow([r, "{{:.4f}}".format(d)])
cmd.delete("all")
"""


def all_together_script(model_pdbs: list[Path], crystal_pdb: Path,
                        png_out: Path) -> str:
    loads = "\n".join(
        f'cmd.load(r"{p}", "mdl_{i+1:02d}")'
        for i, p in enumerate(model_pdbs)
    )
    aligns = "\n".join(
        f'cmd.align("mdl_{i+1:02d} and name CA", "crystal and name CA")'
        for i in range(len(model_pdbs))
    )
    n = len(model_pdbs)
    colors = ["red", "orange", "yellow", "limegreen", "green", "teal",
              "cyan", "blue", "purple", "magenta"]
    if n > len(colors):
        colors = colors * ((n // len(colors)) + 1)
    color_lines = "\n".join(
        f'cmd.color("{colors[i]}", "mdl_{i+1:02d}")' for i in range(n)
    )
    return f"""from pymol import cmd
cmd.load(r"{crystal_pdb}", "crystal")
{loads}
cmd.hide("everything")
cmd.show("cartoon")
cmd.color("white", "crystal")
{color_lines}
{aligns}
cmd.bg_color("white")
cmd.set("ray_opaque_background", "on")
cmd.set("cartoon_transparency", 0.05)
cmd.orient()
cmd.zoom("crystal", buffer=2)
cmd.ray(1600, 1200)
cmd.png(r"{png_out}", dpi=150)
"""


def render_heatmap_png(csv_path: Path, png_path: Path, model_label: str) -> None:
    import csv as _csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    resi: list[int] = []
    dist: list[float] = []
    with open(csv_path) as fh:
        r = _csv.DictReader(fh)
        for row in r:
            resi.append(int(row["resi"]))
            dist.append(float(row["ca_dist"]))
    if not resi:
        return
    fig, ax = plt.subplots(figsize=(12, 4.5))
    # color bars by magnitude
    colors = ["#3b82c4" if d < 1 else "#f6c945" if d < 2 else "#d6604d" for d in dist]
    ax.bar(resi, dist, color=colors, width=1.0)
    ax.axhline(2.0, color="grey", linestyle="--", linewidth=0.8, label="2 Å threshold")
    ax.axhline(1.0, color="grey", linestyle=":", linewidth=0.8, label="1 Å threshold")
    ax.set_xlabel("Residue number")
    ax.set_ylabel("Cα distance (Å)")
    ax.set_title(f"Per-residue Cα distance: {model_label} vs 1HVY chain A")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.6)
    fig.tight_layout()
    fig.savefig(png_path, dpi=220)
    plt.close(fig)


def main() -> int:
    STEP5_DIR.mkdir(parents=True, exist_ok=True)
    crystal = STEP1_DIR / "1hvy_chainA.pdb"
    if not crystal.exists():
        LOG.error("Missing crystal: %s", crystal)
        return 2

    models = list_models()
    if not models:
        LOG.error("No models found in %s", MODELS_DIR)
        return 3
    LOG.info("Comparing %d models against %s", len(models), crystal.name)

    rmsd_rows: list[dict] = []
    for i, m in enumerate(models, start=1):
        png_out = STEP5_DIR / f"pairwise_model{i:02d}.png"
        rmsd_json = STEP5_DIR / f"_pairwise_model{i:02d}.json"
        script = pairwise_script(m, crystal, png_out, rmsd_json)
        ok, _ = run_pymol_py(script, f"pair{i:02d}")
        if not ok or not png_ok(png_out):
            LOG.warning("Pairwise PNG failed/small for model %d", i)
        try:
            data = json.loads(rmsd_json.read_text())
            rmsd = float(data["rmsd"])
            n_atoms = int(data["n_atoms_aligned"])
        except Exception as e:
            LOG.warning("RMSD parse failed for %d: %s", i, e)
            rmsd, n_atoms = float("nan"), 0
        rmsd_rows.append({
            "model_id": i,
            "model_pdb": m.name,
            "rmsd_to_crystal": rmsd,
            "n_atoms_aligned": n_atoms,
        })
        LOG.info("Model %02d  RMSD=%.3f Å  n_atoms=%d", i, rmsd, n_atoms)

        heat_csv = STEP5_DIR / f"per_residue_dist_model{i:02d}.csv"
        h_script = heatmap_dump_script(m, crystal, heat_csv)
        ok2, _ = run_pymol_py(h_script, f"heat{i:02d}")
        if ok2 and heat_csv.exists():
            heat_png = STEP5_DIR / f"per_residue_heatmap_model{i:02d}.png"
            try:
                render_heatmap_png(heat_csv, heat_png, m.stem)
            except Exception as e:
                LOG.warning("Heatmap render failed for %d: %s", i, e)

    csv_path = STEP5_DIR / "rmsd_per_model.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["model_id", "model_pdb",
                                            "rmsd_to_crystal", "n_atoms_aligned"])
        w.writeheader()
        w.writerows(rmsd_rows)
    LOG.info("Wrote %s", csv_path)

    all_png = STEP5_DIR / "all_models_overlay.png"
    all_script = all_together_script(models, crystal, all_png)
    ok, _ = run_pymol_py(all_script, "all_overlay")
    if not png_ok(all_png):
        LOG.warning("All-models overlay PNG too small or missing")

    valid = [r for r in rmsd_rows if r["rmsd_to_crystal"] == r["rmsd_to_crystal"]]
    best_rmsd = min(valid, key=lambda r: r["rmsd_to_crystal"]) if valid else None

    best_dope_id = None
    best_dope_pdb = None
    best_dope_val = None
    scores_csv = STEP4_DIR / "scores.csv"
    if scores_csv.exists():
        with open(scores_csv) as fh:
            rows = list(csv.DictReader(fh))
            if rows:
                best = min(rows, key=lambda r: float(r["DOPE"]))
                best_dope_id = int(best["model_id"])
                best_dope_pdb = best["model_pdb"]
                best_dope_val = float(best["DOPE"])

    summary = {
        "best_by_rmsd": best_rmsd,
        "best_by_dope": {
            "model_id": best_dope_id,
            "model_pdb": best_dope_pdb,
            "DOPE": best_dope_val,
        },
        "n_models": len(rmsd_rows),
    }
    (STEP5_DIR / "best_summary.json").write_text(json.dumps(summary, indent=2))
    LOG.info("Best by RMSD: %s   Best by DOPE: %s",
             best_rmsd, summary["best_by_dope"])

    for p in STEP5_DIR.glob("_tmp_*.py"):
        p.unlink()
    LOG.info("STEP 5 OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
