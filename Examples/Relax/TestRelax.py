"""
Here we use a water toy model to test the automatic
relaxation of the expert mode sscha.

This is performed with the EMT calculator. 
This is just for fun
"""


import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble, sscha.SchaMinimizer, sscha.Relax
import sscha.Utilities

from CP2K_toy_model_calculator import CP2K_water_calculator
import os

dynmat = CC.Phonons.Phonons("final_dyn")

calc = CP2K_water_calculator()
atms = dynmat.structure.get_ase_atoms()
atms.set_calculator(calc)

# Prepare the saving of the frequencies
freq_saving = sscha.Utilities.IOInfo()
freq_saving.SetupSaving("data/Total_freqs.dat")

# Prepare the Relax
ensemble = sscha.Ensemble.Ensemble(dynmat, 0)
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

# Setup a small step (we are using few configurations)
minim.min_step_dyn = 0.5
minim.min_step_struc = 0.1
minim.gradi_op = "all"
minim.meaningful_factor = 0.1

# With few configurations it is possible to have imaginary frequencies
# We deactivate the preconditioning and set up the nonlinear change of variable
# This will smooth the minimization
#minim.precond_dyn = True
#minim.root_representation = "root4"

# We setup the SSCHA relaxation.
N_CONFIGS = 1000
MAX_ITERATIONS = 40
relax = sscha.Relax.SSCHA(minim, calc, N_CONFIGS, MAX_ITERATIONS)

# Avoid to save the ensemble at each new population
relax.save_ensemble = False

# Setup the method to save the frequencies at each step
# to check what the minimization is doing
relax.setup_custom_functions(custom_function_post=freq_saving.CFP_SaveFrequencies)

# Start the relaxation
relax.relax()

# Save the final result
relax.minim.dyn.save_qe("final_dyn_new")

