import sys, os
import numpy as np
import cellconstructor as CC
import cellconstructor.Phonons
import ase, ase.calculators, ase.calculators.emt

import sscha, sscha.Ensemble
import pytest
import numpy as np

@pytest.mark.julia
def test_forces(verbose = False):
    np.random.seed(0)
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)

    temperature = 300

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
        supercell = (4,4,4))
 
    dynmat.AdjustQStar()
    dynmat.Symmetrize()
    dynmat.ForcePositiveDefinite()
    

    # Compute the gradient of the ensemble
    ensemble = sscha.Ensemble.Ensemble(dynmat, temperature)
    ensemble.generate(50)
    ensemble.compute_ensemble(calculator)

    # Get the average forces
    av_forces_fourier = ensemble.get_fourier_forces(False)
    av_forces = ensemble.get_average_forces(False).ravel()
 
    if verbose:
        print("Average forces:")
        print(av_forces)
        print("Average forces fourier:")
        print(av_forces_fourier)

    assert np.allclose(av_forces, av_forces_fourier, atol= 1e-8)

    # Update the temperature
    ensemble.update_weights_fourier(dynmat, temperature + 10)

    # New forces
    av_forces_fourier = ensemble.get_fourier_forces(False)

    ensemble.update_weights(dynmat, temperature + 10)
    av_forces = ensemble.get_average_forces(False).ravel()

    if verbose:
        print()
        print("AFTER UPDATE TEMPERATURE")
        print("Average forces:")
        print(av_forces)
        print("Average forces fourier:")
        print(av_forces_fourier)

    assert np.allclose(av_forces, av_forces_fourier, atol= 1e-8)

    
if __name__ == "__main__":
    test_forces(True)
