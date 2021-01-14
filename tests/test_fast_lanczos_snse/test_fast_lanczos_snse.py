from __future__ import print_function

import numpy as np

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.DynamicalLanczos
import sscha.Ensemble

import scipy, scipy.sparse

from sscha.Parallel import pprint as print

import sys, os

def test_lanczos_symmetries(temperature = 250, N = 10000):
    """
    Test if the ensemble uwrapping works as the symmetries in the polarization
    space.
    """
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    PATH_TO_DYN = "../../Examples/TestLanczosFiniteTemperature/SnTe_test/T_{}_N_{}".format(temperature, N)

    if not os.path.exists(PATH_TO_DYN):
        raise IOError("Error, the path {} does not exist, please change temperature or N".format(PATH_TO_DYN))

    # Load the ensemble
    dyn = CC.Phonons.Phonons(os.path.join(PATH_TO_DYN, "SnTe_final"), 3)
    ens = sscha.Ensemble.Ensemble(dyn, temperature, dyn.GetSupercell())
    ens.load_bin(os.path.join(PATH_TO_DYN, "ensemble"), 1)

    dirname = "SC_T_{:d}_N_{:d}".format(temperature, N)
    if not os.path.exists(dirname):
        os.makedirs(dirname)


    # Get only the first configuration (this is the best way to test everything)
    first_config = np.zeros(ens.N, dtype = bool)
    first_config[:1] = True
    new_ens = ens.split(first_config)

    TEST_V3 = True
    TEST_V4 = False

    lanc1 = sscha.DynamicalLanczos.Lanczos(new_ens, unwrap_symmetries=True)
    lanc1.ignore_v3 = not TEST_V3
    lanc1.ignore_v4 = not TEST_V4
    lanc1.init()

    lanc2 = sscha.DynamicalLanczos.Lanczos(new_ens, unwrap_symmetries=False)
    lanc2.ignore_v3 = not TEST_V3
    lanc2.ignore_v4 = not TEST_V4
    lanc2.init()

    # Get a random vector
    v1 = np.random.uniform(size = lanc1.psi.shape)
    lanc1.psi[:] = v1
    lanc2.psi[:] = v1

    # Apply L
    print("APPLY UNWRAPPED.")
    sys.stdout.flush()
    w1 = lanc1.apply_anharmonic_FT()

    print("ALL THE UNWRAPPED FORCES AGAIN:")
    print(lanc2.Y)
    
    print("APPLY DIRECT SYMMETRIZATION")
    sys.stdout.flush()
    w2 = lanc2.apply_anharmonic_FT()

    disp = np.max(np.abs(w2 - w1))
    assert disp < 1e-12, "Error, the symmetrization in LANCZOS does not work properly."
    
    

def test_lanczos_snse(temperature = 250, N = 10000):
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    PATH_TO_DYN = "../../Examples/TestLanczosFiniteTemperature/SnTe_test/T_{}_N_{}".format(temperature, N)

    if not os.path.exists(PATH_TO_DYN):
        raise IOError("Error, the path {} does not exist, please change temperature or N".format(PATH_TO_DYN))

    # Load the ensemble
    dyn = CC.Phonons.Phonons(os.path.join(PATH_TO_DYN, "SnTe_final"), 3)
    ens = sscha.Ensemble.Ensemble(dyn, temperature, dyn.GetSupercell())
    ens.load_bin(os.path.join(PATH_TO_DYN, "ensemble"), 1)

    dirname = "SC_T_{:d}_N_{:d}".format(temperature, N)
    if not os.path.exists(dirname):
        os.makedirs(dirname)


    # Get only the first 1000 configurations
    first_1000_configs = np.zeros(ens.N, dtype = bool)
    first_1000_configs[:200] = True
    new_ens = ens.split(first_1000_configs)

    # Get the hessian in the standard way
    hessian = new_ens.get_free_energy_hessian(include_v4 = True, use_symmetries = True)
    hessian.save_qe(os.path.join(dirname, "hessian_v4_"))
    dyn.save_qe(os.path.join(dirname, "sscha"))
    

    # Prepare the Lanczos
    lanczos = sscha.DynamicalLanczos.Lanczos(new_ens, unwrap_symmetries=False)
    lanczos.ignore_v3 = False
    lanczos.ignore_v4 = False
    lanczos.init()

    # print("Computing L matrix...")
    # L = sscha.DynamicalLanczos.get_full_L_matrix(lanczos, transpose = False)
    # np.save(os.path.join(dirname, "L_good_final.npy"), L)
    # print("Saved L matrix.")
    # exit()
    
    # # Try to diagonalize using scipy to check eigenvalues and eigenvectors:
    # print("Going into scipy")
    # eigvals, eigvects = scipy.sparse.linalg.eigs(lanczos.L_linop, sigma = 0)

    # w = np.sign(np.real(eigvals)) * np.sqrt(np.abs(eigvals)) * CC.Units.RY_TO_CM
    # print("Frequencies: {} cm-1".format(w))
    # exit()

    # # Get the green function with lanczos
    # for i in range(3):
    #     dirnew = os.path.join(dirname, "mode_{:02d}".format(i))
    #     if not os.path.exists(dirnew):
    #         os.makedirs(dirnew)

    #     lanczos.reset()
    #     lanczos.prepare_mode(i)
    #     lanczos.run_FT(10, save_dir= dirnew, verbose=  True, n_rep_orth = 0)
    
    # exit()

    lanczos.prepare_mode(10)
    lanczos.run_FT(30, save_dir= dirname, verbose=  True, n_rep_orth = 1)
    gf = lanczos.get_green_function_continued_fraction(np.array([0]), smearing = 0, use_terminator= False)
    w2 = np.real(1 / gf[0])

    # Get frequency:
    print("Frequency: {}".format(np.sign(w2) * np.sqrt(np.abs(w2)) * CC.Units.RY_TO_CM))

    lanczos.save_status("LanczosSnTe_v4.npz")

    # Test the frequency
    w, pols = hessian.DiagonalizeSupercell()
    assert np.abs(w[3] - np.sqrt(np.abs(w2))) * CC.Units.RY_TO_CM < 1e-3, "Error, the lanczos w -> 0 does not match the hessian matrix with Bianco algorithm" 

    # # Get the free energy hessian
    # hessian = lanczos.run_biconjugate_gradient(algorithm = "bicgstab", use_preconditioning = True, tol = 1e-12)

    # np.savetxt("hessian_lanczos.dat", hessian)
    # w, pols = np.linalg.eigh(hessian)

    # print("Frequencies:")
    # print("\n".join(["{} cm-1".format(np.sign(x) * np.sqrt(np.abs(x)) * CC.Units.RY_TO_CM) for x in w]))

if __name__ == "__main__":
    test_lanczos_symmetries()
    test_lanczos_snse()
    
