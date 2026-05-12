
load /Users/ario/conserved_site_project/03_structure/protein_cons.pdb, prot
load /Users/ario/conserved_site_project/03_structure/ligand.pdb, lig
hide everything
show cartoon, prot
spectrum b, blue_white_red, prot
show sticks, lig
color green, lig
bg_color white
orient
zoom prot, 4
ray 1600, 1200
png 03_conservation.png
