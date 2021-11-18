# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import division

import sys, os
import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble
import sscha.SchaMinimizer

"""
This test makes a simple relaxation of the sample ensemble
provided within this distribution
"""

def test_gradient_comparison(verbose = False):
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    DATA_PATH = "../../Examples/ensemble_data_test/"

    dyn_start = CC.Phonons.Phonons(os.path.join(DATA_PATH, "dyn"))


    # Perform the minimization
    ens = sscha.Ensemble.Ensemble(dyn_start, 0, dyn_start.GetSupercell())
    ens.load(DATA_PATH, 2, 1000)

    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
    minim.minim_struct = False
    minim.min_step_dyn = 0.5
    minim.min_step_struc = 0.5
    minim.meaningful_factor = 1e-10
    minim.neglect_symmetries = True
    minim.max_ka = 5

    class CG:
        def __init__(self):
            self.ka = 0
        def compare_gradients(self, dyn_grad, struct_grad, minim):
            ka = self.ka
            
            if not os.path.exists("grad_{}.dat".format(ka)):
                np.savetxt("grad_{}.dat".format(ka), dyn_grad[0,:,:])
                np.savetxt("dyn_{}.dat".format(ka), minim.dyn.dynmats[0])
            else:
                correct_grad = np.loadtxt("grad_{}.dat".format(ka))
                correct_dyn = np.loadtxt("dyn_{}.dat".format(ka))
                diff = np.max(np.abs(dyn_grad - correct_grad))
                diffd = np.max(np.abs(minim.dyn.dynmats[0] - correct_dyn))
                print("KA = {} | Delta g  = {} | Delta d = {}".format(ka, diff, diffd))
            self.ka += 1

    cg = CG()
    
    minim.init()
    minim.run(custom_function_gradient = cg.compare_gradients)
    minim.finalize()

    

if __name__ == "__main__":
    test_gradient_comparison(True)
