
reinitialize
load /Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/04_allosteric/apo_for_fpocket.pdb, recep
load /Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/04_allosteric/poses/cav2_CID5564_pose1.pdb, lig
# styling
bg_color white
hide everything
show cartoon, recep
color grey80, recep and chain A
color grey50, recep and chain B
set cartoon_transparency, 0.20
# cavity residues — within 6 Å of the ligand
select cav_res, byres (recep within 6.0 of lig)
show sticks, cav_res and not name H*
util.cnc cav_res
color yelloworange, cav_res and elem C
# ligand
show sticks, lig
color cyan, lig and elem C
color red, lig and elem O
color blue, lig and elem N
color yellow, lig and elem S
color magenta, lig and elem F
color forest, lig and elem Cl
# polar contacts (H-bond < 3.5 Å between N/O atoms)
distance hbond, (recep and (elem N+O) and not elem H), (lig and (elem N+O)), 3.5
hide labels, hbond
color yellow, hbond
set dash_color, yellow
set dash_width, 2.5
# label key residues
label cav_res and name CA, "%s/%s%d" % (chain, resn, resi)
set label_size, 11
set label_color, black
# camera
orient lig
zoom lig, 12
ray 1200, 1200
png /Users/ario/Downloads/aminak-inhibitor/14_inhibitor_design/04_allosteric/poses/cav2_CID5564.png, dpi=140
quit
