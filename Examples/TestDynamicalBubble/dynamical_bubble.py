import sys
import numpy as np

import cellconstructor as CC
import cellconstructor.Phonons

from matplotlib.pyplot import *


import sscha, sscha.Ensemble
import sys

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

print ""

# Now get the spectral function
w_s = np.linspace(0, 3500, 10000)
nat = dyn.structure.N_atoms
spectral_function = np.zeros(len(w_s), dtype = np.float64)
I = np.eye(3*nat, dtype  = np.complex128)


self_energies = ens.get_dynamical_bubble(q_gamma, w_s/CC.Phonons.RY_TO_CM, smearing = 10/CC.Phonons.RY_TO_CM)

# Get the mass array
m = dyn.structure.get_masses_array()
m = np.tile(m, (3, 1)).T.ravel() # Stretch in a 3*nat array
for i,_w_ in enumerate(w_s):
    w = _w_ / CC.Phonons.RY_TO_CM

    # Get the dynamical force constant
    new_dyn = dyn.dynmats[0] + self_energies[i]

    # Convert in a dynamical matrix
    new_dyn /= np.sqrt(np.einsum("a,b", m, m))
    
    g_minus_one = w**2 * I - new_dyn

    green_matrix = np.linalg.inv(g_minus_one)
    # - Im Tr(G)
    spectral_function[i] = np.einsum("aa", -np.imag(green_matrix))

    if i % 5 == 0:
        sys.stdout.write("\rStatus %d %%" %( i * 100 / len(w_s)))
print ""
    
# Plot the spectral function
figure(dpi = 170)
title("Spectral function")
xlabel("Frequency [cm-1]")
ylabel("Spectral function")
plot(w_s, spectral_function)
tight_layout()
show()
