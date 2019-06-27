import sys
import numpy as np
from numpy import *

import cellconstructor as CC
import cellconstructor.Phonons

from matplotlib.pyplot import *


import sscha, sscha.Ensemble
import sscha.Dynamical

SUPERCELL = (1,1,1)
T = 0
N_RANDOM = 1000
POPULATION = 2
N_W = 100
W_MIN = 2800
W_MAX = 3400
GAMMA = 25 / CC.Phonons.RY_TO_CM

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
# w_s = np.linspace(0, 3500, 10000)
# nat = dyn.structure.N_atoms
# spectral_function = np.zeros(len(w_s), dtype = np.float64)
# I = np.eye(3*nat, dtype  = np.complex128)


# self_energies = ens.get_dynamical_bubble(q_gamma, w_s/CC.Phonons.RY_TO_CM, smearing = 10/CC.Phonons.RY_TO_CM)

# # Get the mass array
# m = dyn.structure.get_masses_array()
# m = np.tile(m, (3, 1)).T.ravel() # Stretch in a 3*nat array
# for i,_w_ in enumerate(w_s):
#     w = _w_ / CC.Phonons.RY_TO_CM

#     # Get the dynamical force constant
#     new_dyn = dyn.dynmats[0] + self_energies[i]

#     # Convert in a dynamical matrix
#     new_dyn /= np.sqrt(np.einsum("a,b", m, m))
    
#     g_minus_one = w**2 * I - new_dyn

#     green_matrix = np.linalg.inv(g_minus_one)
#     # - Im Tr(G)
#     spectral_function[i] = np.einsum("aa", -np.imag(green_matrix))
 
#     if i % 5 == 0:
#         sys.stdout.write("\rStatus %d %%" %( i * 100 / len(w_s)))
# print ""

# Get the dynamical response
self_energies = []
w_s = linspace(W_MIN, W_MAX, N_W) / CC.Phonons.RY_TO_CM
print "Computing the v3 with dynamical correction..."
self_energies = ens.get_odd_correction(include_v4 = False,
                                       frequencies = w_s,
                                       smearing = GAMMA,
                                       return_only_correction = True)
print ""
print "Computing the v4 with the dynamical corrections..."
self_energiesv4 = ens.get_odd_correction(include_v4 = True,
                                       frequencies = w_s,
                                       smearing = GAMMA,
                                       return_only_correction = True)
print "Computing the spectral function..."
spectral_function = sscha.Dynamical.get_spectral_function(dyn, (1,1,1),
                                                       self_energies, w_s)
spectral_functionv4 = sscha.Dynamical.get_spectral_function(dyn, (1,1,1),
                                                            self_energiesv4, w_s)

# Get the frequencies
w_array = []
new_dyn = dyn.Copy()
for i in range(len(w_s)):
    new_dyn.dynmats[0] = dyn.dynmats[0] + self_energies[i]
    w, pols= new_dyn.DyagDinQ(0)
    w_array.append(w)
w_array = np.array(w_array)

# Plot the spectral function
figure(dpi = 170)
title("Spectral function")
xlabel("Frequency [cm-1]")
ylabel("Spectral function")
plot(w_s * CC.Phonons.RY_TO_CM, spectral_function, label="only bubble")
plot(w_s * CC.Phonons.RY_TO_CM, spectral_functionv4, label = "whole spectrum")
legend()
tight_layout()
savefig("bubble_vs_v4_dynamical.eps")

figure(dpi = 170)
title("Evolution of the poles")
xlabel("Perturbation frequency [cm-1]")
ylabel("Poles [cm-1]")
for i in range(dyn.structure.N_atoms*3):
    plot(w_s * CC.Phonons.RY_TO_CM, w_array[:, i] * CC.Phonons.RY_TO_CM, color = "b")
plot(w_s * CC.Phonons.RY_TO_CM, w_s * CC.Phonons.RY_TO_CM, color = "r")
tight_layout()
show()
