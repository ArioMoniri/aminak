#!/usr/bin/env python3
"""
Build per-complex 3Dmol.js interactive viewer HTML files for every PDB in
07b_mut_docking_v2/viewer_files/  (and 06b_docking_wt_v2/wt_*.pdb if present),
plus an index page that lists them all.

Output:
    viewers/index.html
    viewers/<mut>_<cond>_complex.html
    viewers/<mut>_<cond>_top_pose.html
    viewers/wt_apo.html, viewers/wt_holo.html
    viewers/dump_ligand.html

Run from project root:
    python scripts/v2/build_viewers.py
"""
from __future__ import annotations
import os, sys, pathlib, json, html
import datetime

ROOT = pathlib.Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
VIEWERS = ROOT / "viewers"
VIEWERS.mkdir(exist_ok=True)
# Prefer v5 (in-place cofactor reprotonation, 0.000 Å heavy-atom RMSD, clash-free);
# fall back to v4 → v3 → v2.
CANDIDATES = [
    (ROOT / "07e_mut_docking_v5" / "viewer_files", ROOT / "06e_docking_wt_v5", "v5"),
    (ROOT / "07d_mut_docking_v4" / "viewer_files", ROOT / "06d_docking_wt_v4", "v4"),
    (ROOT / "07c_mut_docking_v3" / "viewer_files", ROOT / "06c_docking_wt_v3", "v3"),
    (ROOT / "07b_mut_docking_v2" / "viewer_files", ROOT / "06b_docking_wt_v2", "v2"),
]
MUTDIR = WTDIR = None
SOURCE_TAG = "?"
for m, w, tag in CANDIDATES:
    if m.exists() and (w / "wt_holo.pdbqt").exists():
        MUTDIR, WTDIR, SOURCE_TAG = m, w, tag
        break
if MUTDIR is None:
    MUTDIR, WTDIR, SOURCE_TAG = CANDIDATES[-1][0], CANDIDATES[-1][1], "v2 (fallback)"
LIGDIR = ROOT / "05b_ligand_v2"
STRUCTDIR = ROOT / "03b_structure_v2"
print(f"viewer source: {SOURCE_TAG} (MUTDIR={MUTDIR.parent.name}/{MUTDIR.name}, WTDIR={WTDIR.name})")

CDN_3DMOL = "https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"


def render_html(title: str, pdb_path: pathlib.Path, *,
                lig_resn: str = "UMP",
                cof_resn: str | None = "D16",
                description: str = "",
                width: int = 800, height: int = 560) -> str:
    """Generate a self-contained 3Dmol.js viewer HTML for one PDB."""
    if not pdb_path.exists():
        return f"<html><body>Missing: {pdb_path}</body></html>"
    pdb_text = pdb_path.read_text()
    # Escape backticks/backslashes inside JS template literal
    pdb_js = pdb_text.replace("\\", "\\\\").replace("`", "\\`")
    cof_block = ""
    if cof_resn:
        cof_block = f"viewer.setStyle({{resn: '{cof_resn}'}}, {{stick: {{colorscheme: 'cyanCarbon', radius: 0.20}}}});"
    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<title>{html.escape(title)}</title>
<script src="{CDN_3DMOL}"></script>
<style>
 body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif;
        margin: 0; padding: 0 16px 12px 16px; background:#0c1116; color:#e8eef5; }}
 h1 {{ font-size: 17px; font-weight: 600; margin: 12px 0 4px; color:#9ad5ff; }}
 .meta {{ font-size: 12px; color:#aab8c8; margin: 0 0 8px; line-height:1.4; }}
 .panel {{ background:#1a2230; border:1px solid #283344; border-radius:8px;
          padding:0; overflow:hidden; }}
 .legend {{ font-size: 11.5px; color:#cad5e6; padding: 6px 12px; border-top:1px solid #283344; }}
 .legend span {{ display:inline-block; padding: 2px 8px; border-radius: 4px;
                margin: 2px 4px 2px 0; }}
 .nav a {{ color:#9ad5ff; margin-right: 12px; font-size: 12px; text-decoration:none; }}
 .nav a:hover {{ text-decoration:underline; }}
</style>
</head><body>
<div class="nav" style="margin: 10px 0 4px;">
  <a href="index.html">← Index</a>
  <a href="https://github.com/ArioMoniri/aminak">GitHub</a>
</div>
<h1>{html.escape(title)}</h1>
<div class="meta">{html.escape(description)}</div>
<div class="panel">
 <div id="viewer" style="width:100%; max-width:{width}px; height:{height}px; margin:auto; position:relative;"></div>
 <div class="legend">
  <span style="background:#1d3a5a;color:#cfe6ff;">Receptor cartoon</span>
  <span style="background:#5a3a1d;color:#ffe6cf;">Active-site sticks</span>
  <span style="background:#3b1d5a;color:#dccfff;">dUMP (UMP)</span>
  {('<span style="background:#1d5a4a;color:#cfffe6;">Cofactor (D16/raltitrexed)</span>' if cof_resn else '')}
 </div>
</div>
<script>
const pdbText = `{pdb_js}`;
const viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "0x0c1116"}});
viewer.addModel(pdbText, "pdb");
// Receptor cartoon
viewer.setStyle({{}}, {{cartoon: {{colorscheme: 'spectrum', opacity: 0.80}}}});
// Active-site residues coloured by chain
const activeResi = [50,80,87,109,135,170,175,176,195,196,214,215,217,218,221,225,226,258];
viewer.setStyle({{resi: activeResi}}, {{stick: {{colorscheme: 'amino', radius: 0.18}},
                                       cartoon: {{colorscheme: 'spectrum'}}}});
// Ligand (dUMP)
viewer.setStyle({{resn: '{lig_resn}'}}, {{stick: {{colorscheme: 'magentaCarbon', radius: 0.22}}}});
{cof_block}
viewer.zoomTo({{resn: '{lig_resn}'}});
viewer.zoom(0.85);
viewer.render();
</script>
</body></html>
"""


def make_index(entries: list[tuple[str, str, str]]):
    """entries: list of (group, label, href)"""
    today = datetime.date.today().isoformat()
    cards = []
    by_group: dict[str, list] = {}
    for grp, label, href in entries:
        by_group.setdefault(grp, []).append((label, href))
    for grp in sorted(by_group):
        items = "".join(
            f'<li><a href="{html.escape(h)}">{html.escape(l)}</a></li>'
            for l, h in sorted(by_group[grp])
        )
        cards.append(f"<section><h2>{html.escape(grp)}</h2><ul>{items}</ul></section>")
    body = "\n".join(cards)
    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<title>TYMS / dUMP — Interactive 3D Viewers</title>
<style>
 body {{ font-family: system-ui, -apple-system, Segoe UI, sans-serif;
        margin: 0; padding: 24px; background:#0c1116; color:#e8eef5;
        max-width: 980px; margin: 0 auto; }}
 h1 {{ color:#9ad5ff; font-weight: 600; }}
 h2 {{ color:#a4d5b8; font-size: 16px; border-bottom:1px solid #283344;
       padding-bottom: 4px; margin-top: 28px; }}
 ul {{ columns: 2; column-gap: 24px; padding-left: 18px; }}
 li {{ margin: 4px 0; break-inside: avoid; }}
 a  {{ color:#cfe6ff; text-decoration: none; font-size: 13.5px; }}
 a:hover {{ text-decoration: underline; }}
 .note {{ background:#15202b; border-left: 3px solid #9ad5ff;
          padding: 10px 14px; font-size: 13px; line-height: 1.5;
          color: #cad5e6; border-radius: 0 6px 6px 0; }}
</style>
</head><body>
<h1>Interactive 3D viewers — TYMS / dUMP / 1HVY</h1>
<p class="note">
 All viewers are powered by <a href="https://3dmol.csb.pitt.edu">3Dmol.js</a>.
 Each PDB file is embedded inline (no external loading), so these pages work both locally
 and via GitHub Pages. Drag to rotate, scroll to zoom, right-drag to translate.
 Receptor is shown as cartoon (spectrum-coloured), active-site residues as sticks,
 dUMP substrate in magenta, and (where present) raltitrexed cofactor in cyan.
</p>
<p class="note" style="border-color:#a4d5b8;">
 <strong>Conditions:</strong> <code>apo</code> = receptor with cofactor pocket emptied
 (raltitrexed removed); <code>holo</code> = raltitrexed retained in chain A. Both runs use
 the dimer receptor (chains A + B) and identical Vina box & seed.
</p>
{body}
<p style="font-size: 11px; color: #6b7785; margin-top: 32px;">
 Generated {today} · Repo: <a href="https://github.com/ArioMoniri/aminak">ArioMoniri/aminak</a>
</p>
</body></html>
"""


def main():
    entries: list[tuple[str, str, str]] = []

    # --- WT and reference structures ---
    for fn, label, group, desc in [
        ("wt_apo_complex.pdb", "WT (apo) + top dUMP pose", "Wild-type",
         "Wild-type dimer receptor (chains A+B), cofactor pocket emptied. The top Vina pose for dUMP is concatenated."),
        ("wt_holo_complex.pdb", "WT (holo) + top dUMP pose", "Wild-type",
         "Wild-type dimer receptor with raltitrexed retained in chain A. Top Vina pose of dUMP."),
    ]:
        p = WTDIR / fn
        if p.exists():
            out = VIEWERS / f"{p.stem}.html"
            out.write_text(render_html(label, p, description=desc))
            entries.append((group, label, out.name))

    # --- structure prep snapshots ---
    for fn, label, group, desc in [
        ("protein_dimer_h.pdb", "Cleaned dimer (chains A+B, protonated)", "Reference structures",
         "Cleaned 1HVY dimer with hydrogens, CME43→CYS restoration, ligand and cofactor removed."),
        ("ligand.pdb", "Crystal dUMP (chain A, reference)", "Reference structures",
         "dUMP from 1HVY chain A used as the docking target and RMSD reference."),
    ]:
        p = STRUCTDIR / fn
        if p.exists():
            out = VIEWERS / f"struct_{p.stem}.html"
            out.write_text(render_html(label, p, description=desc, cof_resn=None))
            entries.append((group, label, out.name))

    # --- ligand multi-format previews ---
    p = LIGDIR / "dump.pdb"
    if p.exists():
        out = VIEWERS / "ligand_dump.html"
        out.write_text(render_html("dUMP — natural substrate", p,
                                   description="2'-deoxyuridine 5'-monophosphate. Available in PDB, MOL2, SDF, PDBQT under 05b_ligand_v2/.",
                                   cof_resn=None))
        entries.append(("Reference structures", "dUMP ligand (multi-format)", out.name))

    # --- mutant complexes ---
    if MUTDIR.exists():
        for pdb in sorted(MUTDIR.glob("*_complex.pdb")):
            stem = pdb.stem  # e.g. C195A_apo_complex
            parts = stem.replace("_complex", "").split("_")
            # parts may be [MUT, COND] or [MUT1, MUT2, COND]
            if len(parts) >= 2 and parts[-1] in ("apo", "holo"):
                cond = parts[-1]
                mut = "_".join(parts[:-1])
            else:
                cond, mut = "?", stem
            label = f"{mut} ({cond}) + top dUMP pose"
            desc = (f"Mutant {mut} ({'apo cofactor pocket' if cond == 'apo' else 'raltitrexed retained (holo)'}). "
                    f"Top Vina pose of dUMP shown alongside the mutated receptor.")
            out = VIEWERS / f"{stem}.html"
            out.write_text(render_html(label, pdb, description=desc))
            grp = "Single mutants" if "_" not in mut else "Double mutants"
            if mut.startswith("CTRL"):
                grp = "Control"
            entries.append((grp, label, out.name))

        for pdb in sorted(MUTDIR.glob("*_top_pose.pdb")):
            stem = pdb.stem
            parts = stem.replace("_top_pose", "").split("_")
            if len(parts) >= 2 and parts[-1] in ("apo", "holo"):
                cond = parts[-1]
                mut = "_".join(parts[:-1])
            else:
                cond, mut = "?", stem
            label = f"{mut} ({cond}) — dUMP top pose only"
            desc = "Just the top docked dUMP pose (no receptor) for fast comparison."
            out = VIEWERS / f"{stem}.html"
            out.write_text(render_html(label, pdb, description=desc, cof_resn=None))
            entries.append(("Top-pose only", label, out.name))

    # --- index page ---
    (VIEWERS / "index.html").write_text(make_index(entries))
    print(f"wrote {len(entries)} viewer pages + index.html → {VIEWERS}")


if __name__ == "__main__":
    main()
