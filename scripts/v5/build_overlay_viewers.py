#!/usr/bin/env python3
"""
Build Modeller↔crystal overlay viewers — each page loads BOTH the Modeller
model and the 1HVY crystal chain A into one 3Dmol scene, coloured
differently (model = green cartoon, crystal = magenta cartoon).

Output: viewers/modeller_overlay_model{01..10}.html
Plus an "all 10 overlay" page: viewers/modeller_overlay_all.html.
"""
from __future__ import annotations
import os, pathlib, html

ROOT    = pathlib.Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
CRYSTAL = ROOT / "10_modeller" / "01_clean_pdb" / "1hvy_chainA.pdb"
MODELS  = ROOT / "10_modeller" / "04_modeller_run" / "models"
OUT_DIR = ROOT / "viewers"

CDN = "https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"
GHPAGES = "https://ariomoniri.github.io/aminak"
CRYSTAL_TEXT = CRYSTAL.read_text() if CRYSTAL.exists() else ""

def pair_html(model_id: str, model_path: pathlib.Path) -> str:
    if not model_path.exists():
        return ""
    model_text = model_path.read_text()
    cry_js = CRYSTAL_TEXT.replace("\\", "\\\\").replace("`", "\\`")
    mdl_js = model_text.replace("\\", "\\\\").replace("`", "\\`")
    title = f"Modeller {model_id} ↔ 1HVY crystal — overlay"
    desc  = ("Side-by-side overlay of the Modeller homology model (green) and the 1HVY "
             "experimental crystal chain A (magenta). Drag-rotate, scroll-zoom, toggle "
             "surface, spin. The closer the two cartoons sit, the smaller the Cα RMSD.")
    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<title>{html.escape(title)}</title>
<script src="{CDN}"></script>
<style>
 body {{ font-family: system-ui, sans-serif; margin: 0; padding: 0 16px 12px; background:#0c1116; color:#e8eef5; }}
 h1 {{ font-size: 17px; color:#9ad5ff; margin: 12px 0 4px; }}
 .meta {{ font-size: 12px; color:#aab8c8; margin: 0 0 8px; }}
 .panel {{ background:#1a2230; border:1px solid #283344; border-radius:8px; overflow:hidden; }}
 .legend {{ font-size: 11.5px; color:#cad5e6; padding: 6px 12px; border-top:1px solid #283344; }}
 .legend span {{ display:inline-block; padding: 2px 8px; border-radius:4px; margin: 2px 4px 2px 0; }}
 .controls {{ padding: 6px 12px; border-top:1px solid #283344; font-size: 12px; color:#cad5e6; }}
 .controls button {{ background:#283b53; color:#cfe6ff; border:1px solid #3a4d65;
                     border-radius:4px; padding:4px 10px; margin:2px 4px 2px 0; cursor:pointer; font-size:11.5px; }}
 .controls button:hover {{ background:#3a4d65; }}
 .nav a {{ color:#9ad5ff; margin-right:12px; font-size:12px; text-decoration:none; }}
</style>
</head><body>
<div class="nav" style="margin:10px 0 4px;">
  <a href="index.html">← Index</a>
  <a href="{GHPAGES}">GitHub Pages</a>
  <a href="https://github.com/ArioMoniri/aminak/tree/main/10_modeller">Phase 6 sources</a>
</div>
<h1>{html.escape(title)}</h1>
<div class="meta">{html.escape(desc)}</div>
<div class="panel">
 <div id="viewer" style="width:100%; max-width:900px; height:600px; margin:auto; position:relative;"></div>
 <div class="legend">
  <span style="background:#5a1d3a; color:#ffcde6;">Crystal (1HVY chain A) — magenta cartoon</span>
  <span style="background:#1d5a3a; color:#cdffe6;">Modeller model — green cartoon</span>
  <span style="background:#3b1d5a; color:#dccfff;">Cα backbone</span>
 </div>
 <div class="controls">
  <strong>Toggle:</strong>
  <button onclick="window.setCartoon()">Cartoons only</button>
  <button onclick="window.setRibbon()">Ribbon</button>
  <button onclick="window.setBackbone()">Cα trace</button>
  <button onclick="window.zoomAll()">Zoom to fit</button>
  <button onclick="window.spinOn()">Spin</button>
  <button onclick="window.spinOff()">Stop spin</button>
 </div>
</div>
<script>
// Use window-bound viewer + functions so onclick handlers always resolve correctly.
window._v = $3Dmol.createViewer("viewer", {{backgroundColor: "0x0c1116"}});
const crystalText = `{cry_js}`;
const modelText   = `{mdl_js}`;
window._v.addModel(crystalText, "pdb");   // model index 0 — crystal
window._v.addModel(modelText, "pdb");     // model index 1 — Modeller model

function _withV(fn) {{
  if (!window._v) {{ console.warn("3Dmol viewer not ready"); return; }}
  try {{ fn(window._v); }} catch (e) {{ console.error("viewer error:", e); }}
}}
window.setCartoon = function() {{ _withV(function(v){{
  v.setStyle({{model:0}}, {{cartoon: {{color: '#d62895', opacity: 0.85}}}});
  v.setStyle({{model:1}}, {{cartoon: {{color: '#2ca02c', opacity: 0.85}}}});
  v.render();
}}); }};
window.setRibbon = function() {{ _withV(function(v){{
  v.setStyle({{model:0}}, {{cartoon: {{style:'ribbon', color:'#d62895'}}}});
  v.setStyle({{model:1}}, {{cartoon: {{style:'ribbon', color:'#2ca02c'}}}});
  v.render();
}}); }};
window.setBackbone = function() {{ _withV(function(v){{
  v.setStyle({{model:0}}, {{line: {{color:'#d62895', linewidth:2.5}}}});
  v.setStyle({{model:1}}, {{line: {{color:'#2ca02c', linewidth:2.5}}}});
  v.render();
}}); }};
window.zoomAll = function() {{ _withV(function(v){{ v.zoomTo(); v.zoom(0.9); v.render(); }}); }};
window.spinOn  = function() {{ _withV(function(v){{ v.spin(true); }}); }};
window.spinOff = function() {{ _withV(function(v){{ v.spin(false); }}); }};

// Initial render: cartoons + zoom to fit.
window.setCartoon();
window.zoomAll();
window._v.render();
</script>
</body></html>
"""

def make_all_overlay() -> str:
    if not CRYSTAL_TEXT: return ""
    # Build per-model JS embed
    blocks = []
    palette = ["#1f77b4","#2ca02c","#9467bd","#8c564b","#e377c2","#7f7f7f",
               "#bcbd22","#17becf","#aec7e8","#ffbb78"]
    cry_js = CRYSTAL_TEXT.replace("\\","\\\\").replace("`","\\`")
    for i in range(1, 11):
        mpath = MODELS / f"target.B999900{i:02d}.pdb"
        if not mpath.exists(): continue
        mtext = mpath.read_text().replace("\\","\\\\").replace("`","\\`")
        col = palette[(i-1)%len(palette)]
        blocks.append(f"""
const m{i}Text = `{mtext}`;
viewer.addModel(m{i}Text, "pdb");
viewer.setStyle({{model: {i}}}, {{cartoon: {{color: '{col}', opacity: 0.55}}}});""")
    inner = "".join(blocks)
    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"><title>All 10 Modeller models + 1HVY crystal — overlay</title>
<script src="{CDN}"></script>
<style>
 body {{ font-family: system-ui, sans-serif; margin: 0; padding: 0 16px 12px; background:#0c1116; color:#e8eef5; }}
 h1 {{ font-size: 17px; color:#9ad5ff; margin: 12px 0 4px; }}
 .panel {{ background:#1a2230; border:1px solid #283344; border-radius:8px; overflow:hidden; }}
 .legend {{ font-size: 11.5px; color:#cad5e6; padding: 8px 12px; border-top:1px solid #283344; }}
 .legend span {{ display:inline-block; padding: 2px 8px; border-radius:4px; margin:1px 2px; }}
 .nav a {{ color:#9ad5ff; margin-right:12px; font-size:12px; text-decoration:none; }}
</style></head><body>
<div class="nav"><a href="index.html">← Index</a></div>
<h1>All 10 Modeller models + 1HVY crystal — overlay</h1>
<div class="panel">
 <div id="viewer" style="width:100%; max-width:1100px; height:680px; margin:auto;"></div>
 <div class="legend">
  <span style="background:#5a1d3a; color:#ffcde6;">Crystal (1HVY chain A) — magenta cartoon</span>
  Models 1-10 in distinct colours; opacity 55% to make overlap visible.
 </div>
</div>
<script>
const viewer = $3Dmol.createViewer("viewer", {{backgroundColor:"0x0c1116"}});
const crystalText = `{cry_js}`;
viewer.addModel(crystalText, "pdb");
viewer.setStyle({{model:0}}, {{cartoon: {{color: '#d62895', opacity: 0.85}}}});
{inner}
viewer.zoomTo();
viewer.zoom(0.95);
viewer.render();
</script></body></html>"""

OUT_DIR.mkdir(parents=True, exist_ok=True)
n = 0
for i in range(1, 11):
    mpath = MODELS / f"target.B999900{i:02d}.pdb"
    out = OUT_DIR / f"modeller_overlay_model{i:02d}.html"
    out.write_text(pair_html(f"model {i}", mpath))
    n += 1
    print(f"  wrote {out.name} ({out.stat().st_size//1024} KB)")
all_out = OUT_DIR / "modeller_overlay_all.html"
all_out.write_text(make_all_overlay())
print(f"  wrote {all_out.name} ({all_out.stat().st_size//1024} KB)")
print(f"\n{n+1} overlay pages total.")
