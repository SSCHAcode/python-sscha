import time
import pytest
import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import ase, ase.calculators, ase.calculators.emt
import sys, os
import numpy as np

@pytest.mark.julia 
def test_upsilon_fourier(verbose = False):
    # Fix the seed to assure reproducibility 
    np.random.seed(0)

    # Set the current working directory
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    supercell = (4, 4, 4)
    # Load gold but build a crazy dynamical matrix just to test a low symmetry group
    # R3m (without inversion)
    struct = CC.Structure.Structure(2)
    a_param = 4
    struct.unit_cell = np.eye(3) * a_param
    struct.atoms[0] = "Au"
    struct.atoms[1] = "Ag"
    struct.coords[1, :] = np.ones(3) * a_param / 2 + 0.2
    struct.build_masses()
    
    calculator = ase.calculators.emt.EMT()

    # Get a dynamical matrix
    dynmat = CC.Phonons.compute_phonons_finite_displacements(
        struct,
        calculator, 
        supercell = supercell)
 
    dynmat.AdjustQStar()
    dynmat.Symmetrize()
    dynmat.ForcePositiveDefinite()
    
    nq = len(dynmat.q_tot)
    
    temperature = 300.0
    
    # Get the Upsilon Matrix
    Y_supercell = dynmat.GetUpsilonMatrix(T = temperature)

    # Get the polarization vectors
    w_q = np.zeros((6, nq), dtype = np.float64, order = "F")
    pols_q = np.zeros((6, 6, nq), dtype = np.complex128, order = "F")

    for i in range(nq):
        w_q[:, i], pols_q[:, :, i] = dynmat.DyagDinQ(i)

    masses = dynmat.structure.get_masses_array()

    Y_q_julia = sscha.Ensemble._wrapper_julia_get_upsilon_q(
        w_q, 
        pols_q,
        masses, 
        temperature
    )

    # Perform the fourier transform
    Y_fourier = CC.Phonons.GetDynQFromFCSupercell(
        Y_supercell,
        np.array(dynmat.q_tot),
        dynmat.structure,
        dynmat.structure.generate_supercell(supercell)
    )

    # Test the result
    for i in range(nq):
        if verbose:
            print()
            print("Testing q = %d" % i)
            print("Y fourier")
            print(Y_fourier[i, :, :])
            print("Y julia")
            print(Y_q_julia[:, :, i])
            print()
        assert np.max(np.abs(Y_fourier[i, :, :] - Y_q_julia[:, :, i])) < 1e-8, "Error in the fourier transform of the Upsilon matrix at q = %d" % i


if __name__ == "__main__":
    test_upsilon_fourier(verbose = True)