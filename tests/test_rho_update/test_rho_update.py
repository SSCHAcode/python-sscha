# -*- coding: utf-8 -*-

"""
"""
from __future__ import print_function
from __future__ import division

import sys, os

import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble, sscha.Utilities
import sscha.SchaMinimizer

def test_rho_update():
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    # Load the starting dynamical matrix
    dyn = CC.Phonons.Phonons("../../Examples/ensemble_data_test/dyn")

    # Load the ensemble
    N_RAND = 100
    T0 = 0
    POP = 2
    EQ_ENERGY = -144.40680397
    MEANINGFUL = 0.001

    ens = sscha.Ensemble.Ensemble(dyn, T0)
    ens.load("../../Examples/ensemble_data_test", POP, N_RAND)

    final_dyn = CC.Phonons.Phonons("../../Examples/ensemble_data_test/dyn1_population2", full_name = True)
    
    ens.update_weights(final_dyn, T0)

    rho_true = np.loadtxt("rho.dat")
    assert np.max(np.abs(rho_true - ens.rho)) < 1e-6, "Error, the weights are wrong."
    

if __name__ == "__main__":
    test_rho_update()
    
