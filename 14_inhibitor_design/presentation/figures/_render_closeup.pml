
reinitialize
load /Users/ario/Downloads/aminak-inhibitor/06f_receptor_fixed/dimer_noH.pdb, prot
bg_color white
hide everything
remove chain B
show cartoon, prot
color grey70, prot
select cat, prot and resi 195+196+175+176+215+226+258
show sticks, cat and not name H*
util.cnc cat
color salmon, cat and elem C
# only 3 labels — the ones the legend will refer to
label cat and resi 195 and name CA, "C195"
label cat and resi 215 and name CA, "R215"
label cat and resi 258 and name CA, "Y258"
set label_position, (3, 3, 3)
set label_size, 18
set label_color, black
set label_outline_color, white
set label_font_id, 7
set ray_opaque_background, 1
orient cat
zoom cat, 6
ray 1600, 1200
png /Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/presentation/figures/dimer_activesite_clean.png, dpi=160
