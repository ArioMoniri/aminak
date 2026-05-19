
reinitialize
load /Users/ario/Downloads/aminak-inhibitor/06f_receptor_fixed/dimer_noH.pdb, prot
load /Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/04_allosteric/apo_for_fpocket_out/pockets/pocket18_atm.pdb, pocket18
bg_color white
hide everything
show cartoon, prot
color grey85, prot
# cavity-18 residues
select cav, byres prot within 4.0 of pocket18
# OPAQUE wheat surface
show surface, cav
color wheat, cav
set surface_quality, 2
set surface_color, wheat, cav
set transparency, 0.0, cav
set ambient, 0.55
set light_count, 4
set spec_count, 0
set spec_reflect, 0
set ray_shadow, 0
set two_sided_lighting, 0
# loop 181-197 ∩ cavity in salmon
select cavloop, cav and resi 189-197
color salmon, cavloop
# 2 key sticks for orientation, well-placed labels with offset
select hot, prot and chain B and resi 196+200
show sticks, hot and not name H*
util.cnc hot
color forest, hot and elem C
set stick_radius, 0.24
label hot and resi 200 and name CA, "F200"
label hot and resi 196 and name CA, "H196"
set label_position, (5, 3, 0)
set label_size, 22
set label_color, black
set label_outline_color, white
set label_font_id, 7
set label_shadow_mode, 2
set ray_opaque_background, 1
# orient so cavity is FRONT-facing (zoom on the surface centroid + small rotation)
center cav
zoom cav, 4
ray 1600, 1200
png /Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/presentation/figures/cavity18_carve_v3.png, dpi=160
