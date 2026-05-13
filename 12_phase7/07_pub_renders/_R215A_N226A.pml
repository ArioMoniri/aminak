
load /Users/ario/conserved_site_project/07e_mut_docking_v5/viewer_files/R215A_N226A_holo_complex.pdb, cpx
remove resn HOH
bg_color white

# detect ligand residue id
hide everything
# Receptor
select receptor, polymer.protein
show cartoon, receptor
color grey80, receptor
set cartoon_transparency, 0.55, receptor

# Active site (within 5 A of UMP/dUMP)
select ligand, resn UMP+DUP+UPN+DUM
show sticks, ligand
color cyan, ligand
util.cnc ligand

# Cofactor D16
select cof, resn D16+CB3+THF+MTX
show sticks, cof
color teal, cof
util.cnc cof

# Active site residues within 5 A of ligand
select active5, byres (receptor within 5 of ligand)
show sticks, active5
color slate, active5
util.cnc active5
label active5 and name CA, "%s%s*" % (oneletter, resi)

# Mutated residue (sticks, pink)
select mut_res, (resi 215 or resi 226) and receptor
show sticks, mut_res
color magenta, mut_res
util.cnc mut_res

# H-bonds and close contacts: distances <3.5 A between ligand and active5
distance hbonds, (ligand and not elem H), (active5 and not elem H), 3.5
color yellow, hbonds
hide labels, hbonds

# Title via PyMOL set_title
set label_size, 14
set label_color, black
set label_outline_color, white

# View
orient ligand, animate=0
zoom ligand, 7
set ray_shadows, 0
set ambient, 0.35
set ray_trace_mode, 1
set antialias, 2
viewport 1600, 1200

# Header strip via cmd.set_title (only shows in some pymol forks; use a label too)
pseudoatom title_anchor, pos=[0,0,0]
label title_anchor, "R215A_N226A | dVina = +0.77 kcal/mol"
hide labels, title_anchor

ray 1600, 1200
png /Users/ario/conserved_site_project/12_phase7/07_pub_renders/R215A_N226A_pub.png
