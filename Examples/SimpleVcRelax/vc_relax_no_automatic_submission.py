import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Utilities
import sscha.Relax


# Load the dynamical matrix
dyn = CC.Phonons.Phonons("../ensemble_data_test/dyn", nqirr=1)

# Load the ensemble (0 K)
ensemble = sscha.Ensemble.Ensemble(dyn, T0=0)
ensemble.load("../ensemble_data_test/", population=2, N=1000)

# Prepare the minimization
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

# Prepare the auxiliary functions to save the minimization output
ioinfo = sscha.Utilities.IOInfo()
ioinfo.SetupSaving("minim_1")

# Prepare the vc-relax
relax = sscha.Relax.SSCHA(minim, max_pop=1)
relax.setup_custom_functions(custom_function_post=ioinfo.CFP_SaveAll)

# Run the vc-relax
relax.vc_relax(target_press=0, static_bulk_modulus=100, 
               restart_from_ens=True)

# Save the final dynamical matrix
relax.minim.dyn.save_qe("final_dyn")
