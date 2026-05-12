load /Users/ario/conserved_site_project/03b_structure_v2/protein_dimer_h.pdb, prot
wizard mutagenesis
refresh_wizard
cmd.get_wizard().set_mode("ALA")
cmd.get_wizard().do_select("/prot//A/196/")
best_strain = 1e30
best_idx = 0
n_rot = 0
try:
    while True:
        cmd.get_wizard().set_rotamer(n_rot)
        try:
            s = cmd.get_wizard().get_strain()
        except Exception:
            s = 1e30
        if s < best_strain:
            best_strain = s
            best_idx = n_rot
        n_rot += 1
        if n_rot > 30:
            break
except Exception:
    pass
try:
    cmd.get_wizard().set_rotamer(best_idx)
except Exception:
    pass
print("ROTAMER_PICK resi=196 aa=ALA best_idx=" + str(best_idx) + " strain=" + str(best_strain))
cmd.get_wizard().apply()
set_wizard
save /Users/ario/conserved_site_project/07b_mut_docking_v2/H196A/H196A_mut.pdb, prot