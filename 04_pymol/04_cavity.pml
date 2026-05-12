
load /Users/ario/conserved_site_project/03_structure/protein_h.pdb, prot
load /Users/ario/conserved_site_project/03_structure/ligand.pdb, lig
hide everything
show surface, prot
color gray80, prot
set surface_cavity_mode, 1
set surface_cavity_radius, 7
set surface_cavity_cutoff, -7
show sticks, lig
color cyan, lig
bg_color white
orient
zoom prot, 4
ray 1600, 1200
png 04_cavity.png
