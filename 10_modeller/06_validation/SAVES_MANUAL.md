# Manual UCLA SAVES validation

For full PROCHECK / ERRAT / VERIFY3D / WHATCHECK validation, upload
`10_modeller/04_modeller_run/models/best_model.pdb` to <https://saves.mbi.ucla.edu/>.

The SAVES web tool has no programmatic API; the upload must be done
manually OR via the project's computer-use channel.

Local Ramachandran (this folder) is the open-source equivalent of
PROCHECK's basic Ramachandran plot — it uses Biopython's φ/ψ
computation and a hand-drawn favoured/allowed-region overlay.

Steps (manual):
1. Visit https://saves.mbi.ucla.edu/.
2. Upload `best_model.pdb` (≈ a few hundred KB).
3. Choose 'Run all programs' and submit.
4. Download the PROCHECK Ramachandran PDF and ERRAT/VERIFY3D HTMLs;
   place into this folder.
