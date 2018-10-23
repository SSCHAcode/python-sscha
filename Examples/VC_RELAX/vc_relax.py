"""
Here we use a water toy model to test the automatic
relaxation toghether with the unic cell relaxation.

"""


import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sscha.Utilities

from CP2K_toy_model_calculator import CP2K_water_calculator
import os

dynmat = CC.Phonons.Phonons("../ensemble_data_test/dyn")

calc = CP2K_water_calculator()
print calc.implemented_properties
atms = dynmat.structure.get_ase_atoms()
atms.set_calculator(calc)

# Prepare the saving of the frequencies
freq_saving = sscha.Utilities.IOInfo()
freq_saving.SetupSaving("data/Total_freqs.dat")

# Prepare the Relax
ensemble = sscha.Ensemble.Ensemble(dynmat, 0)
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

# Setup a small step (we are using few configurations)
minim.min_step_dyn = 0.05
minim.min_step_struc = 0.05
minim.meaningful_factor = 1

# With few configurations it is possible to have imaginary frequencies
# We deactivate the preconditioning and set up the nonlinear change of variable
# This will smooth the minimization
minim.precond_dyn = True

# We setup the SSCHA relaxation.
N_CONFIGS = 400
MAX_ITERATIONS = 20
relax = sscha.Relax.SSCHA(minim, calc, N_CONFIGS, MAX_ITERATIONS)

# Avoid to save the ensemble at each new population
relax.save_ensemble = False

# Setup the method to save the frequencies at each step
# to check what the minimization is doing
relax.setup_custom_functions(custom_function_post=freq_saving.CFP_SaveFrequencies)


# Start the variable cell relaxation
# Use a trial static bulk modulus of 50 GPa to run the minimization
relax.vc_relax(static_bulk_modulus=50)

# Save the final result
relax.minim.dyn.save_qe("final_dyn")

