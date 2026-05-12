
load /Users/ario/conserved_site_project/03_structure/protein_h.pdb, prot
load /Users/ario/conserved_site_project/03_structure/ligand.pdb, lig
hide everything
show cartoon, prot
color gray70, prot
show sticks, lig
color cyan, lig
select active, prot and resi 80+195+196+217+218+225+226+258
show sticks, active
color orange, active
label active and name CA, "%s%s" % (resn, resi)
set label_size, 16
set label_color, black
bg_color white
orient active
zoom active, 8
ray 1600, 1200
png 02_closeup.png
