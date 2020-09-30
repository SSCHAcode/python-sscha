from __future__ import print_function
from __future__ import division

# Import the numeric library
import numpy as np

# Import the cellconstructor and the SSCHA
import cellconstructor as CC
import cellconstructor.Phonons
import sscha, sscha.Ensemble, sscha.SchaMinimizer


# ========== HERE THE INPUT PARAMETERS ============
DATA_DIR = "../ensemble_data_test" # The position of the ensemble
DYN_NAME = "../ensemble_data_test/dyn" # The name of the dynamical matrix
NQIRR = 1 # How many irreducible q points? (it must match the number of dynamical matrices files)
N_RANDOM = 1000 # The number of random configurations
POPULATION = 2 # The ID of the population (it identifies the ensemble inside DATA_DIR)
OUTPUT_DYN = "final_dyn" # Where to save the output result

# The temperature of the simulation
Tg = 0 # Temperature at which the ensemble was generated
T = 0 # Temperature at which you want to perform the simulation

FIX_ATOM_TYPE = "O" # During the minimization, fix only oxygen atoms

# The minimization steps
MINIM_STEP_DYN = 0.01
MINIM_STEP_STRUC = 0.01

KONG_LIU_RATIO = 0.5 # The Kong-Liu effective sample size ratio after which the minimization is stopped


# ===================================================
INFO = """

We optimize the nuclear wave-function fixing an atomic species.
In the example, we fix {} atoms. 
Both the dynamical matrix and centroids degrees of freedom are constrained.

This is performed by using a custom function
that manipulates the gradient at each minimization step.

""".format(FIX_ATOM_TYPE)
print(INFO)



# Here the script starts.

# We load the dynamical matrix
dyn = CC.Phonons.Phonons(DYN_NAME, nqirr = NQIRR)

# We impose the acoustic sum rule and symmetrization
dyn.Symmetrize()



# ================== THE CONSTRAINING FUNCTION ====================

# Get an array that is True for each atom to be fixed.
fix_struct = np.array([x_at == FIX_ATOM_TYPE for x_at in dyn.structure.atoms])

# Obtain an array that can be applied also on the dyn gradient
# By repeating each value 3 times (the xyz cartesian coordinates)
fix_dyn = np.tile(fix_struct, (3, 1)).T.ravel()

all_grads = []
def fix_atoms(gradient_dyn, gradient_struct):
    # gradient_struct is a (n_at, 3) shaped array,
    # that contains the forces on atoms at each step.
    # We put to zero this force on the atoms we want
    #gradient_struct[fix_struct, :] = 0
    
    # Now the gradient violates translations
    # Compute and subtract the average force
    av_force = np.mean(gradient_struct, axis = 0)
    gradient_struct[:,:] -= np.tile(av_force, (dyn.structure.N_atoms, 1))

    all_grads.append(gradient_struct)

    
    # gradient_dyn is a (nq, 3*n_at, 3*n_at) shaped array.
    # nq is the q points, then there is the dynamical matrix.
    # Here cartesian and atomic indices are contracted.
    # For this reason we created a fix_dyn mask.
    gradient_dyn[:, :, fix_dyn] = 0
    gradient_dyn[:, fix_dyn, :] = 0

    # We now must impose the acoustic sum rule
    # We violated it because we fixed some atoms.
    # It can be imposed in the gamma point
    CC.symmetries.CustomASR(gradient_dyn[0, :, :])
# =================================================================
    

# We load the ensemble
print("Loading the ensemble...")
ensemble = sscha.Ensemble.Ensemble(dyn, Tg, dyn.GetSupercell())
ensemble.load(DATA_DIR, POPULATION, N_RANDOM)

# Update the temperature (if T is different from Tg)
ensemble.update_weights(dyn, T)

# Load the SSCHA minimizer
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)

# Setup the minimization parameters
minim.min_step_dyn = MINIM_STEP_DYN
minim.min_step_struc = MINIM_STEP_STRUC
minim.kong_liu_ratio = KONG_LIU_RATIO

# Setup the minimization
minim.init()
minim.print_info()


# Run the minimization, by constraining the atoms
minim.run( custom_function_gradient = fix_atoms )

# Finalize and save the results
minim.finalize()

minim.dyn.save_qe(OUTPUT_DYN)

minim.plot_results()

