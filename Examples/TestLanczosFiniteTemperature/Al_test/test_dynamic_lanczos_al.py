from __future__ import print_function
from __future__ import division

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.DynamicalLanczos, sscha.Ensemble
import sscha.Parallel
from sscha.Parallel import pprint as print
import numpy as np
import time

import ase, ase.calculators.emt


def compute_difference_sections(psi1, psi2, n_modes):
    """
    Divide the psi1 in sections in order to spot where the difference arise
    """

    start_Y = n_modes
    start_A = start_Y + (n_modes * (n_modes + 1)) // 2
    end_A = start_A + (n_modes * (n_modes + 1)) // 2

    assert end_A == len(psi1)
    assert end_A == len(psi2)


    diff_R = np.max(np.abs(psi1[:start_Y] - psi2[:start_Y]))
    diff_Y = np.max(np.abs(psi1[start_Y:start_A] - psi2[start_Y:start_A]))
    diff_A = np.max(np.abs(psi1[start_A:] - psi2[start_A:]))

    return diff_R, diff_Y, diff_A


def test_dynamic_lanczos_al():
    # Change to the local directory

    dyn = CC.Phonons.Phonons("../Al_222_", 4)
    ens = sscha.Ensemble.Ensemble(dyn, 300, dyn.GetSupercell())

    ens.generate(1000)
    ens.compute_ensemble(ase.calculators.emt.EMT(), compute_stress = False)


    lanczos = sscha.DynamicalLanczos.Lanczos(ens, mode = 2)
    lanczos.ignore_v4 = True
    lanczos.ignore_harmonic = True
    lanczos.ignore_v3 = False
    lanczos.prepare_symmetrization(no_sym = True)

    # Print the symmetrization info
    for i in range(lanczos.n_modes):
        print("Mode {} | freq = {:10.3} cm-1 | Degeneracy = ".format(i, lanczos.w[i] * CC.Units.RY_TO_CM), lanczos.degenerate_space[i])


    # Perform the apply_L_operator by using the MPI
    new_random = np.random.uniform(size = lanczos.psi.shape)

    # Put only in R something:
    #new_random[lanczos.n_modes:] = 0 
    #new_random[:lanczos.n_modes] = 0
    #new_random[lanczos.n_modes + 10:] = 0

    new_random /= np.sqrt(new_random.dot(new_random))
    lanczos.psi[:] = new_random

    new_psi =lanczos.L_linop.matvec(lanczos.psi).copy()

    lanczos.psi[:] = new_random
    new_psi_t =lanczos.L_linop.rmatvec(lanczos.psi).copy()
    #lanczos.psi.dot(lanczos.L_linop).copy()

    #np.save("MPI_psi.npy", new_psi)


    print("Prova...")

    print("Computing the full L operator")
    t1 = time.time()
    L_op = lanczos.get_full_L_operator_FT(verbose = True)
    np.save("L.npy", L_op)
    lanczos.psi[:] = new_random
    new_psi_2 =lanczos.L_linop.matvec(lanczos.psi).copy()

    lanczos.psi[:] = new_random
    new_psi_t2 =lanczos.L_linop.rmatvec(lanczos.psi).copy()
    t2 = time.time()

    print("The time to compute and apply the full L:", t2 - t1, "s")

    diff_v = np.abs(new_psi - new_psi_2)
    np.save("diff_direct.npy", diff_v)
    np.save("psiC_direct.npy", new_psi)
    np.save("psiPy_direct.npy", new_psi_2)
    diff = np.max(diff_v)
    print("Direct difference is {}".format(diff))
    diff_R, diff_Y, diff_A = compute_difference_sections(new_psi, new_psi_2, lanczos.n_modes)
    print("    R diff:", diff_R)
    print("    Y diff:", diff_Y)
    print("    A diff:", diff_A)


    diff_v = np.abs(new_psi_t - new_psi_t2)
    np.save("diff_adjoint.npy", diff_v)
    np.save("psiC_adjoint.npy", new_psi_t)
    np.save("psiPy_adjoint.npy", new_psi_t2)
    diff = np.max(diff_v)

    print("Adjoint difference is {}".format(diff))
    diff_R, diff_Y, diff_A = compute_difference_sections(new_psi_t, new_psi_t2, lanczos.n_modes)
    print("    R diff:", diff_R)
    print("    Y diff:", diff_Y)
    print("    A diff:", diff_A)



if __name__ == "__main__":
    test_dynamic_lanczos_al()
