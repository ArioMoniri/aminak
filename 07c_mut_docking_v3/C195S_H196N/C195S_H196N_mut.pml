load /Users/ario/conserved_site_project/03b_structure_v2/protein_dimer_h.pdb, prot
wizard mutagenesis
refresh_wizard
cmd.get_wizard().do_select("/prot//A/195/")
cmd.get_wizard().set_mode("SER")
cmd.frame(1)
cmd.get_wizard().apply()
cmd.get_wizard().do_select("/prot//A/196/")
cmd.get_wizard().set_mode("ASN")
cmd.frame(1)
cmd.get_wizard().apply()
set_wizard
sel_mut = "resi 195+196"
sel_env = f"byres ((resi 195+196) around 6)"
try:
    cmd.sculpt_activate('prot', state=1)
    for _ in range(3):
        cmd.sculpt_iterate('prot', state=1, cycles=20)
    cmd.sculpt_deactivate('prot')
    print('SCULPT_OK')
except Exception as e:
    print('SCULPT_FAIL: ' + str(e))
save /Users/ario/conserved_site_project/07c_mut_docking_v3/C195S_H196N/C195S_H196N_mut.pdb, prot