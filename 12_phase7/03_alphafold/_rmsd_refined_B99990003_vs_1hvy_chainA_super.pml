
load /Users/ario/conserved_site_project/10_modeller/01_clean_pdb/1hvy_chainA.pdb, target
load /Users/ario/conserved_site_project/10b_modeller_refined/02_refined_models/refined_B99990003.pdb, mobile
remove resn HOH
super mobile and name CA, target and name CA, object=aln_conv
super mobile and name CA, target and name CA, object=aln_raw, cycles=0
