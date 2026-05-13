
load /Users/ario/conserved_site_project/10_modeller/01_clean_pdb/1hvy_chainA.pdb, target
load /Users/ario/conserved_site_project/10b_modeller_refined/02_refined_models/refined_B99990010.pdb, mobile
remove resn HOH
align mobile and name CA, target and name CA, object=aln_conv
align mobile and name CA, target and name CA, object=aln_raw, cycles=0
