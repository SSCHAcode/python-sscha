from __future__ import print_function

"""
In this simple example we test the Lanczos procedure on ice. To see if we get
some meaningful result or not.
"""
import sys, os
import time
import numpy as np
import matplotlib.pyplot as plt

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble, sscha.DynamicalLanczos

T = 0
SUPERCELL = (1,1,1)
DATADIR = "../ensemble_data_test"
POPULATION = 2
NRANDOM = 1000

# Where to store the progress?
SAVE_DIR = "data"

# The frequencies/smearing for the dynamical responce
W_START = 0
W_END = 4000 / CC.Phonons.RY_TO_CM
NW = 5000
SMEARING = 5 / CC.Phonons.RY_TO_CM


# The number of eigenvalues to return
N_VALS = 16
N_ITERS = 33

# If the data dir does not exist, create it
if not os.path.exists(SAVE_DIR):
    os.makedirs(SAVE_DIR)

# Load the original dynamical matrix
dyn = CC.Phonons.Phonons("%s/dyn" % DATADIR)

# Load the original ensemble
ens = sscha.Ensemble.Ensemble(dyn, T, SUPERCELL)
ens.load(DATADIR, POPULATION, NRANDOM)

# Compute the Lanczos matrix
LanczosResponce = sscha.DynamicalLanczos.Lanczos(ens)

# Ignore for now v3 and v4
LanczosResponce.ignore_v3 = True
LanczosResponce.ignore_v4 = True

# Prepare the lanczos algorithm with a random vector
nat = dyn.structure.N_atoms
random_vector = np.random.uniform(size = 3*nat)
random_vector /= np.sqrt(random_vector.dot(random_vector))

LanczosResponce.prepare_perturbation(random_vector)

print ()
print (LanczosResponce.psi)
print ()

print("Preparation compleated.")
print("Running the Lanczos algorithm...")

t1 = time.time()
#LanczosResponce.run_full_diag(N_VALS, n_iter = N_ITERS)
LanczosResponce.run(N_ITERS, SAVE_DIR)
t2 = time.time()

print("Lanczos ended. Time elapsed = %.4f s" % (t2-t1))

print ()
print (LanczosResponce.psi)
print ()

LanczosResponce.psi *= 0
LanczosResponce.psi[0] = 1

LanczosResponce.apply_full_L()
print(LanczosResponce.psi)

# print()
# print("Eigenvalues found:", LanczosResponce.eigvals)

# # Now get the self-energy
# self_energy = LanczosResponce.GetFullSelfEnergy()

# print("Saving the self energy to 'self_energy.dat'")
# np.savetxt("self_energy.dat", self_energy)

# Now get the spectral function
w_array = np.linspace(W_START, W_END, NW)
spectral_function =LanczosResponce.get_spectral_function_from_Lenmann(w_array, SMEARING, True)
# LanczosResponce.GetSupercellSpectralFunctionFromEig(w_array, SMEARING)

print("Saving the spectral function to 'spectral_function.dat'")
np.savetxt("spectral_function.dat", np.transpose([w_array, spectral_function]), header = "W [Ry]; Spectral function")

# Plot the spectral function
plt.plot(w_array * CC.Phonons.RY_TO_CM, spectral_function)
plt.title("Spectral function")
plt.xlabel("Frequency [cm-1]")
plt.ylabel("Spectral function")
plt.tight_layout()
plt.show()
