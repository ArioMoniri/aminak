#!/usr/bin/env python3
"""
Generate docs/assets/repo-visualization.svg as a CLICKABLE treemap.
Every block is wrapped in <a xlink:href="GITHUB_URL"> so clicking opens
the corresponding folder on GitHub.
"""
import os, pathlib, html
import squarify

ROOT = pathlib.Path(".").resolve()
GH = "https://github.com/ArioMoniri/aminak/tree/main"
EXCLUDE = {".git",".venv","__pycache__",".DS_Store"}

def folder_size(p: pathlib.Path) -> int:
    total = 0
    for root, dirs, files in os.walk(p):
        dirs[:] = [d for d in dirs if d not in EXCLUDE]
        for f in files:
            try: total += os.path.getsize(os.path.join(root, f))
            except OSError: pass
    return total

GROUPS = {
    "Phase 1 (v1)":   ("#1f3b5e", ["01_msa","02_active_site","03_structure","04_pymol","05_ligand",
                                   "06_docking_wt","07_mut_docking","08_analysis","09_report"]),
    "Phase 2 (v2)":   ("#2b6f9c", ["01b_msa_v2","02b_active_site_v2","03b_structure_v2","04b_pymol_v2",
                                   "05b_ligand_v2","06b_docking_wt_v2","07b_mut_docking_v2",
                                   "08b_analysis_v2","09b_report_v2"]),
    "Phase 3 (v3)":   ("#2f8a6f", ["06c_docking_wt_v3","07c_mut_docking_v3","08c_analysis_v3","09c_report_v3"]),
    "Phase 4 (v4)":   ("#caa44a", ["03d_structure_v4","05d_ligand_v4","06d_docking_wt_v4","07d_mut_docking_v4",
                                   "08d_analysis_v4","09d_report_v4"]),
    "Phase 5 (v5)":   ("#b8593c", ["03e_structure_v5","06e_docking_wt_v5","07e_mut_docking_v5",
                                   "08e_analysis_v5","09e_report_v5"]),
    "Phase 6 (Modeller)": ("#b9408a", ["10_modeller"]),
    "Enhanced":     ("#a04dca", ["11_enhanced"]),
    "Scripts":      ("#446694", ["scripts"]),
    "Reviews":      ("#7e3b8a", ["reviews","reviews_v2","reviews_v3","reviews_v4","reviews_v5","reviews_phase6"]),
    "Viewers":      ("#d4a857", ["viewers"]),
    "Setup/Logs":   ("#888888", ["00_setup","logs"]),
}

items = []
for grp, (color, names) in GROUPS.items():
    for n in names:
        p = ROOT / n
        if p.exists():
            sz = folder_size(p)
            if sz > 0:
                items.append((grp, color, n, sz))

# Squarify wants normalized values
W, H = 1400, 800
sizes = [it[3] for it in items]
norm  = squarify.normalize_sizes(sizes, W, H)
rects = squarify.squarify(norm, 0, 0, W, H)

# Build SVG with <a> wrappers
parts = []
parts.append(f'<?xml version="1.0" encoding="utf-8"?>\n'
             f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
             f'viewBox="0 0 {W} {H+90}" width="100%" font-family="-apple-system, system-ui, Segoe UI, sans-serif">')
# Title
parts.append(f'<text x="{W/2}" y="30" text-anchor="middle" fill="#1f3b5e" '
             f'font-size="20" font-weight="700">aminak — clickable repository map</text>')
parts.append(f'<text x="{W/2}" y="54" text-anchor="middle" fill="#666" '
             f'font-size="13" font-style="italic">block area = folder size (MB) · click any block to open it on GitHub · '
             f'auto-refreshed on push by repo-visualizer</text>')

# Each rect wrapped in <a>
for (grp, color, name, sz), r in zip(items, rects):
    url = f"{GH}/{name}"
    x, y, w, h = r["x"], r["y"]+70, r["dx"], r["dy"]
    mb = sz/1024/1024
    label = f"{name}\\n{mb:.1f} MB"
    parts.append(f'<a xlink:href="{html.escape(url)}" target="_blank">'
                 f'<title>{html.escape(name)} — {mb:.1f} MB ({grp}) · click to open on GitHub</title>'
                 f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
                 f'fill="{color}" stroke="white" stroke-width="2" opacity="0.94">'
                 f'<set attributeName="opacity" to="1.0" begin="mouseover" end="mouseout"/></rect>')
    # Label only if block is big enough
    if w > 60 and h > 30:
        fs = max(9, min(14, int(min(w,h)/8)))
        cx, cy = x + w/2, y + h/2
        parts.append(f'<text x="{cx:.1f}" y="{cy-fs/2:.1f}" text-anchor="middle" '
                     f'fill="white" font-size="{fs}" font-weight="700" pointer-events="none">'
                     f'{html.escape(name)}</text>')
        parts.append(f'<text x="{cx:.1f}" y="{cy+fs:.1f}" text-anchor="middle" '
                     f'fill="white" font-size="{max(8, fs-2)}" pointer-events="none">'
                     f'{mb:.1f} MB</text>')
    parts.append('</a>')

# Legend
y0 = H + 75
parts.append(f'<text x="20" y="{y0}" fill="#1f3b5e" font-size="12" font-weight="700">Phase legend:</text>')
for i, (grp, (color, _)) in enumerate(GROUPS.items()):
    x0 = 130 + i*120
    parts.append(f'<rect x="{x0}" y="{y0-12}" width="14" height="14" fill="{color}" rx="2"/>')
    parts.append(f'<text x="{x0+20}" y="{y0}" fill="#222" font-size="11">{html.escape(grp)}</text>')
parts.append('</svg>')

out = ROOT / "docs" / "assets" / "repo-visualization.svg"
out.write_text("".join(parts))
print(f"wrote {out} ({out.stat().st_size} bytes, {len(items)} clickable blocks)")
