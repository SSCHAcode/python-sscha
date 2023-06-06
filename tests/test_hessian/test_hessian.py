# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division

import sys, os
import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons
import cellconstructor.Timer

import sscha, sscha.Ensemble
import sscha.SchaMinimizer

"""
This test makes a simple relaxation of the sample ensemble
provided within this distribution
"""

def test_hessian(verbose = False):
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    rho_final = np.loadtxt("../test_simple_relax/final_rho.dat")
    DATA_PATH = "../../Examples/ensemble_data_test/"

    dyn_start = CC.Phonons.Phonons(os.path.join(DATA_PATH, "dyn"))
    dyn_target = CC.Phonons.Phonons(os.path.join(DATA_PATH, "dyn1_population2"), full_name = True)

    # Perform the minimization
    ens = sscha.Ensemble.Ensemble(dyn_start, 0, dyn_start.GetSupercell())
    ens.load(DATA_PATH, 2, 1000)

    #Update the ensemble
    ens.update_weights(dyn_target, 0)

    # Compute the hessian
    timer = CC.Timer.Timer(active=True)
    hessian = timer.execute_timed_function(ens.get_free_energy_hessian, include_v4 = True)

    # Load the reference hessian
    #hessian.save_qe("good_reference")
    hessian_ref = CC.Phonons.Phonons("good_reference")
    assert np.allclose(hessian_ref.dynmats[0], hessian.dynmats[0], atol = 1e-8)

    if verbose:
        timer.print_report()

if __name__ == "__main__":
    test_hessian(True)

