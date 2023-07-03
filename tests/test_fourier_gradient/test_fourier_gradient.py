import time
import pytest
import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import ase, ase.calculators, ase.calculators.emt
import sys, os
import numpy as np

@pytest.mark.julia 
def test_fourier(verbose = False):
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
    
    calculator = ase.calculators.emt.EMT()

    # Get a dynamical matrix
    dynmat = CC.Phonons.compute_phonons_finite_displacements(
        struct,
        calculator, 
        supercell = (3,3,3))
 
    dynmat.AdjustQStar()
    dynmat.Symmetrize()
    dynmat.ForcePositiveDefinite()
    

    # Compute the gradient of the ensemble
    ensemble = sscha.Ensemble.Ensemble(dynmat, 0)
    ensemble.generate(50)
    ensemble.compute_ensemble(calculator)

    # Get the gradient in fourier
    if verbose:
        # Warm up
        ensemble.get_fourier_gradient()

    t1 = time.time()
    grad, err = ensemble.get_fourier_gradient()
    t2 = time.time()

    assert np.max(np.abs(np.imag(grad[0, :, :]))) < 1e-8, "Error, imaginary part non zero"

    # Get the standard gradient
    grad_standard = ensemble.get_preconditioned_gradient()
    t3 = time.time()

    if verbose:
        print("Time fourier:", t2 - t1)
        print("Time standard:", t3 - t2)

    
    for iq, q in enumerate(ensemble.current_dyn.q_tot):
        assert np.max(np.abs(grad[iq, :, :] - grad_standard[iq, :, :])) < 1e-8, "Error, the fourier gradient does not match with benchmark"




if __name__ == "__main__":
    test_fourier(verbose=True)

