load /Users/ario/conserved_site_project/03b_structure_v2/protein_dimer_h.pdb, prot
wizard mutagenesis
refresh_wizard
cmd.get_wizard().do_select("/prot//A/226/")
cmd.get_wizard().set_mode("ASP")
cmd.frame(1)
cmd.get_wizard().apply()
set_wizard
sel_mut = "resi 226"
sel_env = f"byres ((resi 226) around 6)"
try:
    cmd.sculpt_activate('prot', state=1)
    for _ in range(3):
        cmd.sculpt_iterate('prot', state=1, cycles=20)
    cmd.sculpt_deactivate('prot')
    print('SCULPT_OK')
except Exception as e:
    print('SCULPT_FAIL: ' + str(e))
save /Users/ario/conserved_site_project/07c_mut_docking_v3/N226D/N226D_mut.pdb, prot