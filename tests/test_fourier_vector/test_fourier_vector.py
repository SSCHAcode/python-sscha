import time
import pytest
import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import ase, ase.calculators, ase.calculators.emt
import sys, os
import numpy as np

@pytest.mark.julia 
def test_fourier_vector(verbose = False):
    # Fix the seed to assure reproducibility 
    np.random.seed(0)

    # Set the current working directory
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    temperature = 300.0
    supercell = (2,2,2)

    # Load gold but build a crazy dynamical matrix just to test a low symmetry group
    # R3m (without inversion)
    struct = CC.Structure.Structure(2)
    a_param = 4
    struct.unit_cell = np.eye(3) * a_param
    struct.atoms[0] = "Au"
    struct.atoms[1] = "Ag"
    struct.coords[1, :] = np.ones(3) * a_param / 2 + 0.2
    struct.build_masses()

    super_struct, itau = struct.generate_supercell(supercell, get_itau = True)
    itau += 1

    nat_sc = super_struct.N_atoms
    n_random = 1
    vector_sc = np.random.normal(size = (n_random, 3*nat_sc))

    calculator = ase.calculators.emt.EMT()

    # Get a dynamical matrix
    dynmat = CC.Phonons.compute_phonons_finite_displacements(
        struct,
        calculator, 
        supercell=supercell)
 
    dynmat.AdjustQStar()
    dynmat.Symmetrize()
    dynmat.ForcePositiveDefinite()

    q_grid = np.array(dynmat.q_tot)

    # Get the lattice
    R_lat = np.zeros( (nat_sc, 3), dtype = np.float64)
    for i in range(nat_sc):
            R_lat[i,:] = super_struct.coords[i, :] - \
                struct.coords[itau[i] - 1, :]

    vector_q = sscha.Ensemble._wrapper_julia_vector_r2q(
            vector_sc,
            q_grid, 
            itau, 
            R_lat)

    # Return back
    vector_sc_new = sscha.Ensemble._wrapper_julia_vector_q2r(
            vector_q,
            q_grid, 
            itau, 
            R_lat)

    if verbose:
        print("Original Vector:")
        print(vector_sc)
        print("New vector:")
        print(vector_sc_new)

    # Get the Y matrix in q space
    w_r, p_r, w_q, pols_q = dynmat.DiagonalizeSupercell(return_qmodes=True)
    Y_q = sscha.Ensemble._wrapper_julia_get_upsilon_q(
        w_q,
        pols_q,
        struct.get_masses_array(),
        temperature
        )

    # Get Y in real space
    Y_r = dynmat.GetUpsilonMatrix(T = temperature) 

    # Multiply uYu in real space and compare with 
    # The same in fourier space
    Yv_q = sscha.Ensemble._wrapper_julia_matrix_vector_fourier(
            Y_q,
            vector_q)

    # Fourier transform back
    Yv_sc_ft = sscha.Ensemble._wrapper_julia_vector_q2r(
            Yv_q,
            q_grid,
            itau,
            R_lat)

    # Perform the matrix vector multiplication in real space
    Yv_sc_original = vector_sc.dot(Y_r)

    if verbose:
        print("Y q space:")
        print(Y_q)
        print()
        print("Y real space:")
        print(Y_r)
        print()
        print("Yv q:")
        print(Yv_q)
        print()
        print("Yv r ft:")
        print(Yv_sc_ft)
        print("Yv r:")
        print(Yv_sc_original)
            
    assert np.allclose(Yv_sc_original, Yv_sc_ft), "Error in the matrix vector multiplication in q space"





    # Check if they are the same
    assert np.max(np.abs(vector_sc_new - vector_sc)) < 1e-8,  "Error in the fourier transform of the vectors."

    
if __name__ == "__main__":
    test_fourier_vector(verbose=True)

