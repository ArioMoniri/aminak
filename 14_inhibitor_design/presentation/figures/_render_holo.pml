
reinitialize
load /Users/ario/Downloads/aminak-inhibitor/03_structure/1hvy.pdb, full
bg_color white
hide everything
show cartoon, full and polymer
color grey80, full and chain A
color grey60, full and chain B
# dUMP (UMP residue)
select dump, full and resn UMP
show sticks, dump and not name H*
util.cnc dump
color cyan, dump and elem C
# raltitrexed cofactor (resname D16)
select cofac, full and resn D16
show sticks, cofac and not name H*
util.cnc cofac
color hotpink, cofac and elem C
# catalytic residues
select cat, full and chain A and resi 195+196+215 and polymer
show sticks, cat and not name H*
color salmon, cat and elem C
label cat and resi 195 and name CA, "Cys195"
label dump and name N3, "dUMP"
label cofac and name N1 and chain A, "D16"
set label_size, 16
set label_color, black
set label_outline_color, white
set label_font_id, 7
set ray_opaque_background, 1
orient dump or cofac
zoom dump or cofac, 5
ray 1600, 1200
png /Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/presentation/figures/holo_dump_cofactor_clean.png, dpi=160
