
load /Users/ario/conserved_site_project/04b_pymol_v2/protein_dimer_jsd_bfactor.pdb, prot
hide everything
bg_color white
show surface, prot
spectrum b, blue_white_red, prot, 0, 30
load /Users/ario/conserved_site_project/03b_structure_v2/ligand.pdb, ump
show sticks, ump
color cyan, ump
orient prot
zoom prot, 5
ray 1600, 1200
png conservation_surface.png
