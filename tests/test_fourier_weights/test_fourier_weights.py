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
    new_temperature = 310
    n_configs = 2
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
        supercell=supercell)
 
    dynmat.AdjustQStar()
    dynmat.Symmetrize()
    dynmat.ForcePositiveDefinite()
    

    # Compute the gradient of the ensemble
    ensemble = sscha.Ensemble.Ensemble(dynmat, temperature)
    ensemble.generate(n_configs)
    
    #Edit the dynamical matrix
    new_dyn = dynmat.Copy()
    new_dyn.dynmats[-1][:,:] += np.random.normal(
            size=(3 * struct.N_atoms, 3 * struct.N_atoms)) * 1e-4
    new_dyn.structure.coords[1, :] += 0.01
    new_dyn.Symmetrize()
    new_dyn.ForcePositiveDefinite()

    # Warm up
    if verbose:
        ensemble.update_weights_fourier(new_dyn, new_temperature)

    timer = CC.Timer.Timer(active=True)

    timer.execute_timed_function(ensemble.update_weights_fourier, new_dyn, new_temperature)
    # print the displacements
    disp_new_qspace = ensemble.u_disps_qspace.copy()
    disp_new = sscha.Ensemble._wrapper_julia_vector_q2r(
            disp_new_qspace,
            ensemble.q_grid,
            ensemble.itau,
            ensemble.r_lat)

    # Get also sscha forces and energies
    energies = ensemble.sscha_energies.copy()
    forces = ensemble.sscha_forces.copy()

    # get the standard weights
    weights_fourier = ensemble.rho.copy()
    
    timer.execute_timed_function(ensemble.update_weights, new_dyn, new_temperature)
    weights_standard = ensemble.rho.copy()
    disp_old = ensemble.u_disps.copy()
    energies = ensemble.sscha_energies.copy()
    forces = ensemble.sscha_forces.copy()

    if verbose:
        print("New displacements:")
        print(disp_new)
        print("Old displacements:")
        print(disp_old)
        print("New weights:")
        print(weights_fourier)
        print("Old weights:")
        print(weights_standard)
        
    # Check if the displacements are the same
    assert np.allclose(disp_new, disp_old)
    assert np.allclose(energies, energies)
    assert np.allclose(forces, forces)
    assert np.allclose(weights_fourier, weights_standard)
   
    if verbose:
        timer.print_report(verbosity_limit=0)

    # Here, compute the gradient with the two methods
    # and check if they are the same

    gradient_old = ensemble.get_preconditioned_gradient() 
    gradient_new, err = ensemble.get_fourier_gradient()

    if verbose:
        print("Old gradient:")
        print(gradient_old)
        print("New gradient:")
        print(gradient_new)
        print("Delta:")
        print(gradient_old - gradient_new)
        print("Max delta:")
        print(np.max(np.abs(gradient_old - gradient_new)))

    assert np.allclose(gradient_old, gradient_new)



if __name__ == "__main__":
    test_fourier(verbose=True)

