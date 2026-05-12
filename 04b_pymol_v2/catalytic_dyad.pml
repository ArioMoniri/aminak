
load /Users/ario/conserved_site_project/04b_pymol_v2/protein_dimer_jsd_bfactor.pdb, prot
load /Users/ario/conserved_site_project/03b_structure_v2/ligand.pdb, ump
load /Users/ario/conserved_site_project/03b_structure_v2/cofactor_chainA_h.pdb, d16A
load /Users/ario/conserved_site_project/03b_structure_v2/cofactor_chainB_h.pdb, d16B
hide everything
bg_color white
spectrum b, blue_white_red, prot, 0, 30
show cartoon, prot
show sticks, ump
color cyan, ump
show sticks, d16A
color magenta, d16A
show sticks, d16B
color salmon, d16B
select active_A, prot and chain A and resi 50+109+175+176+195+196+214+215+225+226
select active_B, prot and chain B and resi 50+109+175+176+195+196+214+215+225+226
show sticks, active_A
show sticks, active_B

hide sticks, active_B
select dyad, prot and chain A and resi 195+196+175+176+215+226
show sticks, dyad
label dyad and name CA, "%s%s" % (resn, resi)
orient (dyad or ump or d16A)
zoom (dyad or ump or d16A), 3
ray 1600, 1200
png catalytic_dyad.png
