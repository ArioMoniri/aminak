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
 .controls {{ padding: 6px 12px; border-top:1px solid #283344;
              font-size: 12px; color:#cad5e6; }}
 .controls button {{ background:#283b53; color:#cfe6ff; border:1px solid #3a4d65;
                     border-radius: 4px; padding: 4px 10px; margin: 2px 4px 2px 0;
                     cursor: pointer; font-size: 11.5px; }}
 .controls button:hover {{ background:#3a4d65; }}
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
  <span style="background:#1d3a5a;color:#cfe6ff;">Receptor cartoon + transparent surface</span>
  <span style="background:#5a3a1d;color:#ffe6cf;">Active-site sticks (yellow C)</span>
  <span style="background:#3b1d5a;color:#dccfff;">dUMP (UMP — magenta sticks)</span>
  {('<span style="background:#1d5a4a;color:#cfffe6;">Cofactor (D16/raltitrexed — cyan sticks)</span>' if cof_resn else '')}
 </div>
 <div class="controls">
  <strong>Toggle:</strong>
  <button onclick="hideSurface()">Hide surface</button>
  <button onclick="showSurface()">Show surface</button>
  <button onclick="cartoonOnly()">Cartoon only</button>
  <button onclick="zoomLigand()">Zoom to ligand</button>
  <button onclick="spinOn()">Spin</button>
  <button onclick="spinOff()">Stop spin</button>
 </div>
</div>
<script>
function _withV(fn) {{
  if (!window._v) {{ console.warn('viewer not ready'); return; }}
  try {{ fn(window._v); }} catch (e) {{ console.error(e); }}
}}
function hideSurface()  {{ _withV(function(v){{ v.removeAllSurfaces(); v.render(); }}); }}
function showSurface()  {{ _withV(function(v){{
  var p = v.addSurface($3Dmol.SurfaceType.MS, {{opacity:0.30, color:'lightgrey'}}, {{polymer:true}});
  if (p && p.then) p.then(function(){{ v.render(); }}); else v.render();
}}); }}
function cartoonOnly()  {{ _withV(function(v){{
  v.removeAllSurfaces();
  v.setStyle({{polymer:true}}, {{cartoon: {{colorscheme: 'spectrum', opacity: 1.0}}}});
  v.render();
}}); }}
function zoomLigand()   {{ _withV(function(v){{ v.zoomTo({{resn: '{lig_resn}'}}); v.zoom(0.85); v.render(); }}); }}
function spinOn()       {{ _withV(function(v){{ v.spin(true); }}); }}
function spinOff()      {{ _withV(function(v){{ v.spin(false); }}); }}
</script>
<script>
const pdbText = `{pdb_js}`;
const viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "0x0c1116"}});
window._v = viewer;   // expose for the toggle buttons
viewer.addModel(pdbText, "pdb");

// 1. Receptor: cartoon (always visible) + semi-transparent surface (only when toggle is on)
viewer.setStyle({{polymer: true}}, {{cartoon: {{colorscheme: 'spectrum', opacity: 0.85}}}});
viewer.addSurface($3Dmol.SurfaceType.MS,
                  {{opacity: 0.30, color: 'lightgrey'}},
                  {{polymer: true}});

// 2. Active-site residues — sticks, amino-acid colouring
const activeResi = [50,80,87,109,135,170,175,176,195,196,214,215,217,218,221,225,226,258];
viewer.setStyle({{resi: activeResi, polymer: true}}, {{stick: {{colorscheme: 'amino', radius: 0.18}},
                                                      cartoon: {{colorscheme: 'spectrum'}}}});
// Label active-site residues by Cα
viewer.addStyle({{resi: activeResi, atom: 'CA'}}, {{}});
// (3Dmol auto-labels via clicker on Cα; explicit labels on a handful of catalytic residues:)
const catalyticResi = [195, 196, 175, 176, 215, 226];
catalyticResi.forEach(function(r){{
  viewer.addLabel(String(r), {{
     fontSize: 11, fontColor: 'white', backgroundColor: '0x202830',
     backgroundOpacity: 0.75, alignment: 'center'
  }}, {{resi: r, atom: 'CA', polymer: true}});
}});

// 3. Ligand (dUMP) — fat magenta sticks (the user's feature request)
viewer.setStyle({{resn: '{lig_resn}'}}, {{stick: {{colorscheme: 'magentaCarbon', radius: 0.30}}}});
{cof_block}

// 4. Frame on the ligand
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

    # Featured viewers (highlighted at top)
    FEATURED_HREFS = [
        ("WT (apo) + dUMP — reference",                     "wt_apo_complex.html"),
        ("WT (holo) + dUMP — reference",                    "wt_holo_complex.html"),
        ("R215A_N226A holo — top destabiliser",             "R215A_N226A_holo_complex.html"),
        ("H196A holo — catalytic dyad probe",               "H196A_holo_complex.html"),
        ("R175E_R176E holo — phosphate clamp inversion",    "R175E_R176E_holo_complex.html"),
        ("T170A holo — distant-surface negative control",   "T170A_holo_complex.html"),
        ("Modeller best model (by DOPE)",                   "modeller_model03.html"),
        ("Modeller best model (by RMSD vs crystal)",        "modeller_model10.html"),
    ]
    feat_cards = []
    for label, href in FEATURED_HREFS:
        feat_cards.append(
            f'<a class="featured-card" href="{html.escape(href)}">'
            f'<div class="featured-title">{html.escape(label)}</div>'
            f'<div class="featured-href">{html.escape(href)}</div></a>'
        )
    cards.append(
        '<section><h2 style="color:#ffd479;">★ Featured viewers</h2>'
        '<div class="featured-grid">' + "".join(feat_cards) + '</div></section>'
    )

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
 .featured-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
                   gap: 12px; margin-top: 8px; }}
 .featured-card {{ background:#172230; border:1px solid #2c4866; border-radius: 8px;
                   padding: 12px 14px; text-decoration:none; transition: all 0.15s;
                   display:block; }}
 .featured-card:hover {{ border-color:#ffd479; transform: translateY(-1px);
                          background:#1c2a3c; }}
 .featured-title {{ color:#ffd479; font-weight:600; font-size: 13.5px; margin-bottom:3px; }}
 .featured-href  {{ color:#7a8da3; font-size: 11.5px; font-family: ui-monospace, monospace; }}
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
