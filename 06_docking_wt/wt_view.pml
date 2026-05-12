
load /Users/ario/conserved_site_project/03_structure/protein_h.pdb, prot
load /Users/ario/conserved_site_project/06_docking_wt/wt_top_pose.pdb, pose
load /Users/ario/conserved_site_project/03_structure/ligand.pdb, native
hide everything
show cartoon, prot
color gray70, prot
select active, prot and resi 80+195+196+217+218+225+226+258
show sticks, active
color orange, active
show sticks, pose
color magenta, pose
show sticks, native
color cyan, native
bg_color white
orient pose
zoom pose, 8
ray 1600, 1200
png wt_topdock.png
