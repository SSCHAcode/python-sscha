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

    temperature = 300
    n_configs = 100

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
        supercell = (4,4,4))
 
    dynmat.AdjustQStar()
    dynmat.Symmetrize()
    dynmat.ForcePositiveDefinite()
    

    # Compute the gradient of the ensemble
    ensemble = sscha.Ensemble.Ensemble(dynmat, temperature)
    ensemble.generate(n_configs)
    
    #Edit the dynamical matrix
    new_dyn = dynmat.Copy()
    #new_dyn.dynmats[1][0,0] += 0.0003
    new_dyn.structure.coords[1, :] += 0.01
    new_dyn.Symmetrize()
    new_dyn.ForcePositiveDefinite()

    # Warm up
    if verbose:
        ensemble.update_weights_fourier(new_dyn, temperature)

    timer = CC.Timer.Timer(active=True)

    timer.execute_timed_function(ensemble.update_weights_fourier, new_dyn, temperature)

    # get the standard weights
    weights_fourier = ensemble.rho.copy()
    
    timer.execute_timed_function(ensemble.update_weights, new_dyn, temperature)
    weights_standard = ensemble.rho.copy()

    if verbose:
        print("New weights:")
        print(weights_fourier)
        print("Old weights:")
        print(weights_standard)
        timer.print_report(verbosity_limit=0)
        
    

    assert np.allclose(weights_fourier, weights_standard)
   


if __name__ == "__main__":
    test_fourier(verbose=True)

