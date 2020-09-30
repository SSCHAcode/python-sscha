from __future__ import print_function
import sscha, sscha.Ensemble
import cellconstructor as CC
import cellconstructor.Phonons

import numpy as np

import pytest
import sys, os

def test_gradient(verbose = False):
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)


    data_dir = "../../Examples/ensemble_data_test"
    dyn = CC.Phonons.Phonons("{}/dyn".format(data_dir))

    ens = sscha.Ensemble.Ensemble(dyn, 0, (1,1,1))
    # Load the ensemble
    ens.load(data_dir, 2, 1000)

    # Get the gradient
    fc1, err_fc1 = ens.get_preconditioned_gradient(return_error = True)
    fc2, err_fc2 = ens.get_preconditioned_gradient(return_error = True,
                                                   fast_grad = True)


    # Compute manually the gradient
    UpsMat = dyn.GetUpsilonMatrix( T = 0)
    v = ens.u_disps.dot(UpsMat)
    eforc = (ens.forces - ens.sscha_forces).reshape(ens.N, 3*dyn.structure.N_atoms)

    fc_new = np.einsum("ia, ib, i->ab", v, eforc, ens.rho) / np.sum(ens.rho)
    fc_new += fc_new.T
    fc_new /= 2


    error12 = np.max(np.abs(fc1- fc2))
    error13 = np.max(np.abs(fc1 - fc_new))
    error23 = np.max(np.abs(fc2 - fc_new))

    if verbose:
        print("Error 12:", error12)
        print("Error 13:", error13)
        print("Error 23:", error23)
    assert error12 < 1e-12
    assert error13 < 1e-12
    assert error23 < 1e-12



if __name__ == "__main__":
    test_gradient(verbose = True)    
