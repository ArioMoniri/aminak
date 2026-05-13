
bg_color white
load /Users/ario/conserved_site_project/12_phase7/03_alphafold/AF-P04818-F1-model_v6.pdb, af
load /Users/ario/conserved_site_project/10b_modeller_refined/02_refined_models/refined_B99990003.pdb, modeller
load /Users/ario/conserved_site_project/10_modeller/01_clean_pdb/1hvy_chainA.pdb, hvy
remove resn HOH
hide everything
show cartoon
color cyan, af
color forest, modeller
color magenta, hvy
super af and name CA, hvy and name CA
super modeller and name CA, hvy and name CA
set ray_shadows, 0
set cartoon_transparency, 0.0
set ambient, 0.3
viewport 1600, 1200
orient hvy
zoom hvy, 5
ray 1600, 1200
png /Users/ario/conserved_site_project/12_phase7/03_alphafold/triple_overlay.png
