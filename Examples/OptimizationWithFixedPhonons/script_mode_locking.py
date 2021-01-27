import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble
import sscha.SchaMinimizer
import sscha.Utilities

import numpy as np

# Standard input
DYN_NAME = "../ensemble_data_test/dyn"
FINAL_DYN = "final_dyn"
DATA_DIR = "../ensemble_data_test"
NQIRR = 1
POPULATION = 2
N_RANDOM = 1000
Tg = 0
T = 0

# Minimization data
MIN_STEP_DYN = 0.1
MIN_STEP_STRUC = 0.1
KONG_LIU_RATIO = 0.5

SAVE_FREQUENCIES_FILE = "frequencies.dat"

# Move only modes between 30 and 36 (ordered by energy)
FREE_MODE_START = 30
FREE_MODE_END = 36
PROJECT_DYN = True # Lock the dynamical matrix
PROJECT_STRUCTURE = True # Lock also the displacements

def main():
    # Load the dynamical matrix
    dyn = CC.Phonons.Phonons(DYN_NAME, NQIRR)
    dyn.Symmetrize()
    
    # Load the ensemble
    print("Loading the ensemble...")
    ensemble = sscha.Ensemble.Ensemble(dyn, Tg, dyn.GetSupercell())
    ensemble.load(DATA_DIR, POPULATION, N_RANDOM)

    # Prepare the minimizer
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
    minim.min_step_dyn = MIN_STEP_DYN
    minim.min_step_struc = MIN_STEP_STRUC
    minim.kong_liu_ratio = KONG_LIU_RATIO

    minim.init()
    minim.print_info()

    # Prepare the polarization vectors to define the subspace
    mode_lock = sscha.Utilities.ModeProjection(dyn)
    mode_lock.SetupFreeModes(FREE_MODE_START, FREE_MODE_END)
    
    # We configure also a custom function to save the frequencies
    # on a file
    save_freqs = sscha.Utilities.IOInfo()
    save_freqs.SetupSaving(SAVE_FREQUENCIES_FILE)
    

    # Run the minimization by locking all the modes
    # except those inside FREE_MODE_START and FREE_MODE_END
    minim.run(custom_function_gradient = mode_lock.CFG_ProjectOnModes,
              custom_function_post = save_freqs.CFP_SaveFrequencies)

    minim.finalize()
    minim.dyn.save_qe(FINAL_DYN)

    



if __name__ == "__main__":
    main()
