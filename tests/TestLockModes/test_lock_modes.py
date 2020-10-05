# -*- coding: utf-8 -*-

"""
This example shows how to use the expert mode
to allow the sscha relaxation only in a small subset of modes

Here we choose to optimize only the stretching vibrons of the ice XI structure
reported in the example.
"""
from __future__ import print_function
from __future__ import division

import sys, os

import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble, sscha.Utilities
import sscha.SchaMinimizer

def test_lock_the_modes():
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    # Load the starting dynamical matrix
    dyn = CC.Phonons.Phonons("../../Examples/ensemble_data_test/dyn")

    # Load the ensemble
    N_RAND = 1000
    T0 = 0
    POP = 2
    EQ_ENERGY = -144.40680397
    MEANINGFUL = 0.001

    ens = sscha.Ensemble.Ensemble(dyn, T0)
    ens.load("../../Examples/ensemble_data_test", POP, N_RAND)

    # Setup the constraint on the modes
    N_VIBRONS = 4

    # Here we setup the constraint
    mode_cons = sscha.Utilities.ModeProjection(dyn)
    mode_cons.SetupFreeModes(dyn.structure.N_atoms - N_VIBRONS,
                             dyn.structure.N_atoms)

    # We setup also an I/O utility to save the frequencies as a function of the minimization.
    # In this way we are able to see if it is minimizing only the selected vibrons
    IO_freq = sscha.Utilities.IOInfo()
    IO_freq.SetupSaving("frequencies.dat") # The frequencies will be saved at each step of the minimization


    # Now we can setup the standard minimizer
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
    minim.eq_energy = EQ_ENERGY
    minim.min_step_dyn = 0.5
    minim.meaningful_factor = MEANINGFUL

    minim.init(True)

    w_start, p_start = dyn.DiagonalizeSupercell()

    # The constraint must be specified when the run method is called
    #minim.run()
    minim.run(custom_function_gradient = mode_cons.CFG_ProjectOnModes,
              custom_function_post = IO_freq.CFP_SaveFrequencies)


    minim.finalize()
    minim.plot_results("minim.dat", plot = False)

    # Check the frequencies
    w_end, p_end = minim.dyn.DiagonalizeSupercell()


    tol = 1e-7
    print("\n".join(["{:3d} {:16.8f}   {:16.8f} => {:16.8f}".format(i+1, w_start[i], w_end[i], np.abs(w_start[i] - w_end[i])) for i in range(len(w_start))]))
    print("n_vibrons : ", N_VIBRONS)

    n_changes = 0
    for i in range(3, len(w_end)):
        dist = np.min(np.abs(w_start[i] - w_end))
        if dist > tol:
            n_changes += 1

    assert n_changes == N_VIBRONS, "Error, n_changes = {} | expected {}".format(n_changes, N_VIBRONS)


if __name__ == "__main__":
    test_lock_the_modes()
    
