
load /Users/ario/conserved_site_project/07e_mut_docking_v5/viewer_files/Y258F_F225Y_holo_complex.pdb, complex
hide everything
bg_color white
set ray_shadows, 0
set ambient, 0.30
set spec_reflect, 0.15
set surface_quality, 1
set transparency, 0.45
set cartoon_transparency, 0.0

# Receptor surface (chain A core; greyscale)
select receptor, polymer and chain A
show surface, receptor
color grey75, receptor
# Underlying cartoon for reference
show cartoon, receptor
color grey90, receptor

# Ligand sticks (dUMP magenta) - both crystal and docked pose if present
select lig, resn UMP or resn UNK or resn UNL
show sticks, lig
color magenta, lig
set stick_radius, 0.20, lig

# Cofactor if present (D16) - cyan sticks behind
select cof, resn D16
show sticks, cof
color cyan, cof
set stick_radius, 0.16, cof

# Interacting residues (within 4.5 A of ligand)
select interacting, (chain A and resi 50+80+109+135+175+176+195+196+214+215+217+218+221+225+226+258) and polymer
show sticks, interacting
color yellow, interacting and elem c
util.cnc interacting
set stick_radius, 0.16, interacting

# Labels for interacting residues (CA atoms)
label interacting and name CA, "%s%s" % (resn, resi)
set label_size, 12
set label_color, black
set label_position, (0, 1.5, 0)
set label_font_id, 7

# Highlight mutation site(s) in vivid orange
select mutsite, chain A and resi 258+225 and polymer
color orange, mutsite and elem c
util.cnc mutsite
set stick_radius, 0.24, mutsite
label mutsite and name CA, "MUT %s%s" % (resn, resi)

# Centre & orient
orient lig
zoom (lig or interacting), 4
ray 1600, 1200
png /Users/ario/conserved_site_project/11_enhanced/pymol/Y258F_F225Y_holo_render.png, dpi=300

# Also: a "wide" view including dimer interface, semi-transparent
hide everything, receptor
show cartoon, receptor
color spectrum, receptor
show cartoon, polymer and chain B
color grey80, polymer and chain B
show sticks, lig
show sticks, cof
show sticks, interacting
show sticks, mutsite
orient lig
zoom complex, 5
ray 1600, 1200
png /Users/ario/conserved_site_project/11_enhanced/pymol/Y258F_F225Y_holo_render_wide.png, dpi=300
