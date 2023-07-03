import pytest
import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import ase, ase.calculators, ase.calculators.emt
import sys, os
import numpy as np

@pytest.mark.julia 
def test_fourier():
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
        supercell = (1,1,1))
    
    dynmat.Symmetrize()
    dynmat.ForcePositiveDefinite()
    

    # Compute the gradient of the ensemble
    ensemble = sscha.Ensemble.Ensemble(dynmat, 0)
    ensemble.generate(10)
    ensemble.compute_ensemble(calculator)

    # Get the gradient in fourier
    grad, err = ensemble.get_fourier_gradient()

    assert np.max(np.abs(np.imag(grad))) < 1e-8, "Error, imaginary part non zero"

    # Get the standard gradient
    grad_standard = ensemble.get_preconditioned_gradient()

    print("new gradient:")
    print(np.real(grad[0, :, :]))
    print("old gradient:")
    print(grad_standard)

    assert np.max(np.abs(grad - grad_standard)) < 1e-8, "Error, the fourier gradient does not match with benchmark"


if __name__ == "__main__":
    test_fourier()

