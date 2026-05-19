
reinitialize
load /Users/ario/Downloads/aminak-inhibitor/06f_receptor_fixed/dimer_noH.pdb, prot
load /Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/04_allosteric/apo_for_fpocket_out/pockets/pocket18_atm.pdb, pocket18
bg_color white
hide everything, pocket18
show cartoon, prot
color grey80, prot
set cartoon_transparency, 0.4
# residues that belong to cavity 18
select cav, byres prot within 4.0 of pocket18
show surface, cav
set transparency, 0.10, cav
color wheat, cav
# the allosteric-loop subset (181-197)
select loop, cav and resi 189-197
color salmon, loop
# show sticks of the most-conserved contact residues only
select hot, prot and resi 54+87+190+191+196+200+201 and chain B
show sticks, hot and not name H*
color forest, hot and elem C
# only label two key residues to avoid stacking
label hot and resi 196 and name CA, "H196"
label hot and resi 200 and name CA, "F200"
set label_position, (4, 4, 4)
set label_size, 18
set label_color, black
set label_outline_color, white
set label_font_id, 7
set ray_opaque_background, 1
orient cav
zoom cav, 4
ray 1600, 1200
png /Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/presentation/figures/cavity18_carve_clean.png, dpi=160
