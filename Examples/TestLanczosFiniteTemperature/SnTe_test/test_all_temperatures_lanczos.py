from __future__ import print_function
from __future__ import division

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.DynamicalLanczos, sscha.Ensemble
import sscha.Parallel
from sscha.Parallel import pprint as print
import numpy as np
import time

import sys, os

import ase, ase.calculators.emt


SMEARING=0.05 / CC.Units.RY_TO_CM
TERMINATOR = True
W_MIN = 0
W_MAX = 300
N_W = 10000

def test_all_temperature_lanczos(verbose = False):
    # Change to the local directory
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    files = [f for f in os.listdir(".") if f.startswith("T_") and os.path.isdir(f)]
    files.sort()
    
    T = [int(x.split("_")[1]) for x in files]
    N_confs = [int(x.split("_")[-1]) for x in files]

    for i in range(len(files)):
        get_lanczos(T[i], N_confs[i], verbose)


def get_lanczos(T, N_CONF, verbose):
    directory = "T_{}_N_{}".format(T, N_CONF)
    lanc_dir = os.path.join(directory, "lanczos")
    hessian_file = os.path.join(directory, "SnTe_hessian")
    save_file = os.path.join(directory, "green_function.dat")

    if not os.path.exists(lanc_dir):
        if sscha.Parallel.am_i_the_master():
            os.makedirs(lanc_dir)

    # Load the ensemble and set it to point toward the SSCHA result
    dyn = CC.Phonons.Phonons(os.path.join(directory, "SnTe_final"), 3)
    ens = sscha.Ensemble.Ensemble(dyn, T, dyn.GetSupercell())

    ens.load_bin(os.path.join(directory, "ensemble"), 1)
    #ens.update_weights(dyn, T)

    # Get the free energy hessian (Bianco way)
    sscha_hessian = ens.get_free_energy_hessian()
    if verbose and sscha.Parallel.am_i_the_master():
        sscha_hessian.save_qe(hessian_file)

    # For each eigenvector, setup the Lanczos calculation
    w, pols = sscha_hessian.DiagonalizeSupercell()

    # Remove the translations
    ss = dyn.structure.generate_supercell(dyn.GetSupercell())
    trans = CC.Methods.get_translations(pols, ss.get_masses_array())

    w = w[~trans]
    pols = pols[:, ~trans]

    # Setup the Lanczos
    lanczos = sscha.DynamicalLanczos.Lanczos(ens, mode = 2)
    lanczos.ignore_v4 = True
    lanczos.ignore_v3 = False
    lanczos.init()

    # Run the lanczos
    N_ITERS = 30

    # Setup the perturbation along the lowest energy mode
    lanczos.prepare_perturbation(pols[:,0])

    # Run the Lanczos calculation
    lanczos.run_FT(N_ITERS, lanc_dir, verbose)

    # get the full green function
    w_array = np.linspace(W_MIN, W_MAX, N_W)
    w_ry = w_array / CC.Units.RY_TO_CM

    gf = lanczos.get_green_function_continued_fraction(w_array,
                                                       use_terminator = TERMINATOR,
                                                       smearing = SMEARING)

    #spect = lanczos.get_spectral_function_from_Lenmann(w_array, smearing = SMEARING,
    #                                                   use_arnoldi = False)

    # Save the data
    data = np.zeros((N_W, 4), dtype = np.double)
    data[:,0] = w_array
    data[:,1] = -np.imag(gf)
    data[:,2] = np.real(gf)
    #data[:,3] = spect

    if sscha.Parallel.am_i_the_master():
        np.savetxt(save_file, data, header = "freq [cm-1]; -Im G(w); Re G(w); Spectral (w)\nGreen function is the response on the lowest energy mode")
    
    


if __name__ == "__main__":
    test_all_temperature_lanczos(verbose = True)
