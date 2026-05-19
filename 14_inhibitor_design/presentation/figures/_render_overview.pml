
reinitialize
load /Users/ario/Downloads/aminak-inhibitor/06f_receptor_fixed/dimer_noH.pdb, prot
bg_color white
hide everything
show cartoon, prot
color grey80, prot and chain A
color grey60, prot and chain B
# catalytic + clamp residues highlighted as small spheres (NO labels — they collided in v1)
select cat_a, prot and chain A and resi 195+196+175+176+215
show spheres, cat_a and name CA
color salmon, cat_a and name CA
set sphere_scale, 1.5, cat_a and name CA
# single off-image legend label
set ray_opaque_background, 1
orient prot
zoom prot, 8
turn x, 10
ray 1600, 1200
png /Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/presentation/figures/dimer_overview_clean.png, dpi=160
