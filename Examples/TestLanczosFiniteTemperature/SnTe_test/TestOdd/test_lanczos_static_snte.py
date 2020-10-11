from __future__ import print_function
from __future__ import division

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.DynamicalLanczos, sscha.Ensemble
import sscha.Parallel, sscha.SchaMinimizer
from sscha.Parallel import pprint as print
import numpy as np
import time

import sys, os

import ase, ase.calculators.emt


def test_dynamic_lanczos_snte(verbose = False):
    # Change to the local directory
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    # Load the ensemble
    T = 250
    dyn = CC.Phonons.Phonons("../SnTe_sscha", 3)
    ens = sscha.Ensemble.Ensemble(dyn, T, dyn.GetSupercell())

    ens.load_bin("../ensemble", 1)

    # Perform the SSCHA minimization
    minim = sscha.SchaMinimizer.SSCHA_Minimizer(ens)
    minim.meaningful_factor = 0.01

    minim.init()
    minim.run(verbose = 0)
    minim.finalize()

    # Get the free energy hessian (Bianco way)
    sscha_hessian = minim.ensemble.get_free_energy_hessian()
    if verbose:
        sscha_hessian.save_qe("SnTe_hessian")

    # For each eigenvector, setup the Lanczos calculation
    w, pols = sscha_hessian.DiagonalizeSupercell()

    # Remove the translations
    ss = dyn.structure.generate_supercell(dyn.GetSupercell())
    trans = CC.Methods.get_translations(pols, ss.get_masses_array())

    w = w[~trans]
    pols = pols[:, ~trans]

    # Setup the Lanczos
    lanczos = sscha.DynamicalLanczos.Lanczos(minim.ensemble, mode = 2)
    lanczos.ignore_v4 = True
    lanczos.ignore_v3 = False
    lanczos.init()

    # Run the lanczos
    n_modes = len(w)
    N_ITERS = 50

    for i, w_i in enumerate(w):
        # Setup the perturbation along the mode
        lanczos.reset()
        lanczos.prepare_perturbation(pols[:,i])

        # Setup the saving directory if the test is in verbose mode
        save_dir = None
        if verbose:
            save_dir = "lanczos_mode_{}".format(i+3)
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)

        # Run the Lanczos calculation
        lanczos.run_FT(N_ITERS, save_dir, verbose)

        # Get the green function
        g_w = lanczos.get_green_function_continued_fraction([0], use_terminator = False)[0]
        w_lanc = np.sign(g_w) / np.sqrt(np.abs(g_w))
        
        assert np.abs(w_lanc - w_i) * CC.Units.RY_TO_CM < 10, "Error, frequencies {} do not match (B: {} cm-1| L: {} cm-1)".format(i, w_lanc* CC.Units.RY_TO_CM, w_i * CC.Units.RY_TO_CM)
    
    


if __name__ == "__main__":
    test_dynamic_lanczos_snte(verbose = True)
