import numpy as np

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble

SUPERCELL = (1,1,1)
T = 0
N_RANDOM = 1000
POPULATION = 2

# Load the dynamical matrix
dyn = CC.Phonons.Phonons("../ensemble_data_test/dyn", 1)

# Load the ensemble
ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL)
ens.load("../ensemble_data_test", POPULATION, N_RANDOM)

# Add the symmetrization of the ensemble
ens.unwrap_symmetries()

# Compute the static bubble
q_gamma = np.zeros(3, dtype = np.float64)
static_bubble = ens.get_dynamical_bubble(q_gamma, 0)

# Add the bubble to the dynamical matrix
odd_dyn = dyn.Copy()
odd_dyn.dynmats[0] += static_bubble

# Save the new dynamical matrix
odd_dyn.save_qe("static_dyn")

w_dyn, p = odd_dyn.DyagDinQ(0)

# Get the standard static correction
odd_corr = ens.get_odd_correction()
odd_dyn.dynmats[0] = odd_corr

w_odd, p = odd_dyn.DyagDinQ(0)

# Compare the frequencies
print "\n".join(["%16.8f  | %16.8f cm-1" % (w_dyn[x] * CC.Phonons.RY_TO_CM, w_odd[x] * CC.Phonons.RY_TO_CM) for x in range(len(w_dyn))])


