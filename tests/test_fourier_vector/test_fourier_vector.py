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

    # Load gold but build a crazy dynamical matrix just to test a low symmetry group
    # R3m (without inversion)
    struct = CC.Structure.Structure(2)
    a_param = 4
    struct.unit_cell = np.eye(3) * a_param
    struct.atoms[0] = "Au"
    struct.atoms[1] = "Ag"
    struct.coords[1, :] = np.ones(3) * a_param / 2 + 0.2
    struct.build_masses()

    supercell = (4,4,4)
    super_struct, itau = struct.generate_supercell(supercell, get_itau = True)
    itau += 1

    # Define the supercell
    q_grid = np.array(
            CC.symmetries.GetQGrid(struct.unit_cell,
                supercell)
            )

    nat_sc = super_struct.N_atoms
    n_random = 20
    vector_sc = np.random.normal(size = (n_random, 3*nat_sc))

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

    # Check if they are the same
    assert np.max(np.abs(vector_sc_new - vector_sc)) < 1e-8,  "Error in the fourier transform of the vectors."

    
if __name__ == "__main__":
    test_fourier_vector(verbose=True)

