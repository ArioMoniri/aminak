
load /Users/ario/conserved_site_project/03_structure/protein_h.pdb, prot
load /Users/ario/conserved_site_project/03_structure/ligand.pdb, lig
hide everything
show surface, prot
color gray80, prot
show sticks, lig
color cyan, lig
select active, prot and resi 80+195+196+217+218+225+226+258
show sticks, active
color orange, active
set transparency, 0.35
bg_color white
orient
zoom prot, 4
ray 1600, 1200
png 01_overview.png
