from __future__ import print_function
from __future__ import division

import cellconstructor as CC
import cellconstructor.Phonons

import spglib

# Import the modules of the force field
import fforces as ff
import fforces.Calculator
#from unit_cell_snte_calc import UnitCellCalculator

# Import the modules to run the sscha
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities

import numpy as np
import sys,os

# Load the dynamical matrix for the force field
ff_dyn = CC.Phonons.Phonons("ffield_dynq", 3)

# Setup the forcefield with the correct parameters
ff_calculator = ff.Calculator.ToyModelCalculator(ff_dyn)
ff_calculator.type_cal = "pbtex"
ff_calculator.p3 = 0.036475
ff_calculator.p4 = -0.022
ff_calculator.p4x = -0.014

TMIN = 150
TMAX = 250
DT = 10

N_CONFIGS = 10000

def simulate(T, n_confs):
    dir_name = "T_{}_N_{}".format(T, n_confs)
    
    dyn_final_name = os.path.join(dir_name, "SnTe_final")
    ensemble_final_name = os.path.join(dir_name, "ensemble")
    
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    if not os.path.exists(ensemble_final_name):
        os.makedirs(ensemble_final_name)        

    # Load the original dynamical matrix
    sscha_dyn = CC.Phonons.Phonons("SnTe_sscha", 3)
    sscha_dyn.ForcePositiveDefinite()
    sscha_dyn.Symmetrize()

    # Setup the minimization
    print("Generate the ensemble...")
    ens = sscha.Ensemble.Ensemble(sscha_dyn, T, sscha_dyn.GetSupercell())
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
    minim.min_step_dyn = 0.2
    minim.root_representation = "root2"
    #minim.precond_dyn = False
    minim.minim_struct = False

    # Perform the automatic relaxation
    print("Prepare relaxation...")
    relax = sscha.Relax.SSCHA(minim, ff_calculator, N_configs = n_confs, max_pop = 5)
    print("Start the relaxation...")
    relax.relax()

    print("New relaxation")
    relax.max_pop = 12
    relax.minim.min_step_dyn = 1
    relax.relax()

    print("Saving results in {} (population 1) and {}".format(ensemble_final_name, dyn_final_name))
    relax.minim.dyn.save_qe(dyn_final_name)
    relax.minim.ensemble.save_bin(ensemble_final_name, 1)



if __name__ == "__main__":
    T = TMIN
    while T <= TMAX:
        print("Simulating T = {} | N = {}".format(T, N_CONFIGS))
        simulate(T, N_CONFIGS)
        T += DT
    
