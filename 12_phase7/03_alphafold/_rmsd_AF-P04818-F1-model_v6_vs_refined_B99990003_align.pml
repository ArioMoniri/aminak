
load /Users/ario/conserved_site_project/10b_modeller_refined/02_refined_models/refined_B99990003.pdb, target
load /Users/ario/conserved_site_project/12_phase7/03_alphafold/AF-P04818-F1-model_v6.pdb, mobile
remove resn HOH
align mobile and name CA, target and name CA, object=aln_conv
align mobile and name CA, target and name CA, object=aln_raw, cycles=0
