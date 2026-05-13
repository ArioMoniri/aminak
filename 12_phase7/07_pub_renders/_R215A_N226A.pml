
load /Users/ario/conserved_site_project/07e_mut_docking_v5/viewer_files/R215A_N226A_holo_complex.pdb, cpx
remove resn HOH

# White publication background (ray_trace_mode 1 with cartoon outline,
# but force opaque white background and disable shadows)
bg_color white
set opaque_background, on
set ray_opaque_background, on
set ray_shadows, 0
set ray_trace_mode, 1
set ray_trace_color, black
set antialias, 2
set ambient, 0.45
set spec_reflect, 0.10
set cartoon_transparency, 0.50
hide everything

# Suppress hydrogens everywhere — protein + ligand + cofactor — this was
# the dominant clutter source in the previous renders
remove elem H

# Receptor cartoon
select receptor, polymer.protein
show cartoon, receptor
color grey70, receptor

# Substrate ligand (dUMP) — bright green, distinct from cofactor
select ligand, resn UMP+DUP+UPN+DUM
show sticks, ligand
color limegreen, ligand
util.cnc ligand
set stick_radius, 0.30, ligand

# Cofactor (raltitrexed D16) — orange-magenta, clearly different from ligand
select cof, resn D16+CB3+THF+MTX
show sticks, cof
color hotpink, cof
util.cnc cof
set stick_radius, 0.22, cof

# Active-site residues within 5 A of dUMP — wheat sticks (no cofactor neighbours)
select active5, byres (receptor within 5 of ligand)
show sticks, active5
color wheat, active5
util.cnc active5
set stick_radius, 0.16, active5
label active5 and name CA, "%s%s*" % (oneletter, resi)

# Mutated residue (overlays active5; bright orange to differentiate from cofactor pink)
select mut_res, (resi 215 or resi 226) and receptor
show sticks, mut_res
color orange, mut_res
util.cnc mut_res
set stick_radius, 0.30, mut_res

# H-bonds: ONLY polar contacts ligand <-> protein side-chain (N + O atoms),
# distance <3.5 A.  Excludes intra-protein and ligand <-> cofactor lines.
distance hb_lig, (ligand and (elem N or elem O)), (active5 and sidechain and (elem N or elem O)), 3.5
color yellow, hb_lig
hide labels, hb_lig

# Label aesthetics
set label_size, 12
set label_color, black
set label_outline_color, white
set label_font_id, 7

# View centered on the substrate
orient ligand, animate=0
zoom ligand, 7
viewport 1600, 1200

# Title overlay (top-left)
pseudoatom title_anchor, pos=[0,0,0]
label title_anchor, "R215A_N226A | dVina = +0.77 kcal/mol"
hide labels, title_anchor

ray 1600, 1200
png /Users/ario/conserved_site_project/12_phase7/07_pub_renders/R215A_N226A_pub.png
