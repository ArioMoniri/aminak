#!/usr/bin/env python3
"""
Render rotating GIFs of WT (apo + holo) and one key mutant. PyMOL renders
~24 frames around the y-axis, then ImageMagick combines into a small,
README-embeddable looped GIF.

Output: 11_enhanced/gifs/*.gif
"""
import os, pathlib, subprocess, shutil

ROOT = pathlib.Path(os.environ.get("PROJECT_DIR", os.path.expanduser("~/conserved_site_project")))
PYMOL = "/opt/homebrew/bin/pymol"
OUT = ROOT / "11_enhanced" / "gifs"
OUT.mkdir(parents=True, exist_ok=True)
TMP = ROOT / "11_enhanced" / "_frames"
TMP.mkdir(parents=True, exist_ok=True)

N_FRAMES = 30
SIZE = (640, 480)

TARGETS = [
    ("wt_apo", ROOT/"06e_docking_wt_v5"/"wt_apo_complex.pdb",   "WT (apo) + dUMP"),
    ("wt_holo",ROOT/"06e_docking_wt_v5"/"wt_holo_complex.pdb",  "WT (holo) + dUMP"),
    ("R215A_N226A_holo", ROOT/"07e_mut_docking_v5"/"viewer_files"/"R215A_N226A_holo_complex.pdb",
     "R215A_N226A (top destabiliser)"),
]


def render_frames(name: str, pdb: pathlib.Path):
    frame_dir = TMP / name
    if frame_dir.exists(): shutil.rmtree(frame_dir)
    frame_dir.mkdir(parents=True)
    pml = TMP / f"{name}.pml"
    pml.write_text(f"""
load {pdb}, complex
hide everything
bg_color white
set ray_shadows, 0
set ambient, 0.42
set spec_reflect, 0.10
set surface_quality, 1
set transparency, 0.0
set cartoon_transparency, 0.0

select rec, polymer and chain A
show cartoon, rec
color white, rec
show surface, rec
color grey80, rec
set surface_color, grey80, rec

select lig, resn UMP
show sticks, lig
color magenta, lig and elem c
util.cnc lig
set stick_radius, 0.28, lig

select cof, resn D16
show sticks, cof
color cyan, cof and elem c
util.cnc cof
set stick_radius, 0.16, cof

select active, chain A and resi 50+80+109+135+175+176+195+196+214+215+217+218+225+226+258 and polymer
show sticks, active
color yellow, active and elem c
util.cnc active
set stick_radius, 0.16, active

orient lig
zoom (lig or active), 4
ray {SIZE[0]}, {SIZE[1]}

# Render N frames rotating about y-axis
python
import os
N = {N_FRAMES}
for i in range(N):
    cmd.rotate("y", 360.0/N, "complex")
    cmd.ray({SIZE[0]}, {SIZE[1]})
    cmd.png(os.path.join("{frame_dir}", "frame_%02d.png" % i), dpi=120, ray=0)
python end
""")
    subprocess.run([PYMOL, "-cq", str(pml)], capture_output=True, timeout=600)
    return frame_dir


def make_gif(frame_dir: pathlib.Path, out_gif: pathlib.Path):
    cmd = ["/opt/homebrew/bin/magick", "-delay", "8", "-loop", "0",
           str(frame_dir/"frame_*.png"), "-layers", "Optimize", str(out_gif)]
    # magick globs natively
    rc = subprocess.run(["/bin/sh", "-c",
                         f"/opt/homebrew/bin/magick -delay 8 -loop 0 "
                         f"'{frame_dir}/frame_*.png' -layers Optimize '{out_gif}'"],
                        capture_output=True, text=True, timeout=180)
    return rc


for name, pdb, _label in TARGETS:
    if not pdb.exists():
        print(f"[skip] {name}: source PDB missing ({pdb})")
        continue
    print(f"[render] {name} ...")
    fd = render_frames(name, pdb)
    if not any(fd.iterdir()):
        print(f"  no frames produced for {name}")
        continue
    out = OUT / f"{name}.gif"
    rc = make_gif(fd, out)
    if out.exists():
        print(f"  -> {out} ({out.stat().st_size//1024} KB)")
    else:
        print(f"  GIF build failed: {rc.stderr[:300]}")
print("\n== summary ==")
for g in sorted(OUT.glob("*.gif")):
    print(f"  {g.name:40s} {g.stat().st_size//1024:>5} KB")
