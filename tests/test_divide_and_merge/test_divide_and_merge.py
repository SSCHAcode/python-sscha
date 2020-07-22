from __future__ import print_function
import cellconstructor as CC
import cellconstructor.Phonons

import numpy as np

import sscha, sscha.Ensemble
import sys, os

ENS_DIR = "../../Examples/ensemble_data_test"
N_RANDOM = 100
T = 0
__EPS__ = 1e-8

def test_divide_and_merge():
    # Setup the random seed to be able to reproduce the same calculations
    # In different runs
    np.random.seed(0)
    
    total_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(total_path)
    
    # Load the original dynamical matrix
    dyn_start = CC.Phonons.Phonons(os.path.join(ENS_DIR, "dyn"))
    dyn_end = CC.Phonons.Phonons(os.path.join(ENS_DIR, "dyn_population2_"))

    # Generate a random ensemble
    ens_original = sscha.Ensemble.Ensemble(dyn_start, T)
    ens_original.generate(N_RANDOM)
    ens_original.has_stress = False

    # Fake the force computation on half of the ensemble
    ens_original.force_computed[: N_RANDOM // 2] = True
    ens_original.update_weights(dyn_end, T)
    
    # Split the ensemble
    new_ensemble = ens_original.get_noncomputed()

    # Check if the weight is correct
    indexes = np.arange(N_RANDOM)[~ens_original.force_computed]
    g_inds = np.arange(N_RANDOM)[ens_original.force_computed]

    for i, ind in enumerate(indexes):
        TXT="INDEX: {} ({}) old = {} new = {}".format(i, ind,
                                                      new_ensemble.rho[i],
                                                      ens_original.rho[ind])

        assert np.abs(new_ensemble.rho[i] - ens_original.rho[ind]) < __EPS__, TXT

    print("Total weight: {}".format(np.sum(ens_original.rho) / N_RANDOM))

        
    # Remove the non computed ensemble
    all_rhos = ens_original.rho.copy()

    ens_original.remove_noncomputed()
    print("First part: {}".format(np.sum(ens_original.rho) / N_RANDOM))
    print("Second part: {}".format(np.sum(new_ensemble.rho) / N_RANDOM))

    assert np.max(np.abs(ens_original.rho[:] - all_rhos[g_inds])) < __EPS__, "Error while testing the remove_noncomputed function"

    # Test the saving and loading of the two half
    ens_original.save("data", population = 1)
    new_ensemble.save("data", population = 2)

    # Merge them back
    ens_original.merge(new_ensemble)

    assert np.max(np.abs(ens_original.rho - all_rhos)) < __EPS__, "Error during merging the ensembles"

    # Test the merging when loading
    other_ensemble = sscha.Ensemble.Ensemble(dyn_start, T)
    try:
        other_ensemble.load("data", population = 2, N = N_RANDOM //2,  load_noncomputed_ensemble = False)
        assert False, "Error, the previous loading should result in an handled exceptio"
    except IOError:
        other_ensemble.load("data", population = 2, N = N_RANDOM // 2, load_noncomputed_ensemble = True)

    other2 = sscha.Ensemble.Ensemble(dyn_start, T)
    other2.load("data", population = 1, N = N_RANDOM // 2, load_noncomputed_ensemble = False)

    other2.merge(other_ensemble)
    other2.update_weights(dyn_end, T)

    assert np.max(np.abs(other2.rho - all_rhos)) < __EPS__, "Error during merging the loaded ensemble"    
    

if __name__ == "__main__":
    test_divide_and_merge()
    
                                 
                                 
