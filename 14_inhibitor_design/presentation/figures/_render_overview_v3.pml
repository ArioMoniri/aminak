
reinitialize
load /Users/ario/Downloads/aminak-inhibitor/06f_receptor_fixed/dimer_noH.pdb, prot
bg_color white
hide everything
show cartoon, prot
color grey80, prot and chain A
color grey60, prot and chain B
# catalytic + clamp sticks (chain A only, thin so they don't crowd)
select cat, prot and chain A and resi 195+196+215
show sticks, cat and not name H*
color salmon, cat and elem C
util.cnc cat
set stick_radius, 0.20
set ray_opaque_background, 1
orient prot
zoom prot, 10
turn x, 12
ray 1600, 1200
png /Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/presentation/figures/dimer_overview_v3.png, dpi=160
