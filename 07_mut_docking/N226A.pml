
load /Users/ario/conserved_site_project/07_mut_docking/N226A_h.pdb, prot
load /Users/ario/conserved_site_project/07_mut_docking/N226A_top.pdb, pose
hide everything
show cartoon, prot
color gray70, prot
select mut, prot and resi 226
show sticks, mut
color red, mut
show sticks, pose
color magenta, pose
bg_color white
orient pose
zoom pose, 8
ray 1200, 900
png /Users/ario/conserved_site_project/07_mut_docking/N226A.png
