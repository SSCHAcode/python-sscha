from __future__ import print_function

import sys
import numpy as np
from numpy import *

import cellconstructor as CC
import cellconstructor.Phonons

from matplotlib.pyplot import *

import sscha, sscha.Ensemble
import sscha.Dynamical

INFO = """
This very simple script computes the Raman
responce in using the standard dynamical implementation.
The computation is done neglecting the v4.

We test the standard raffaello definition whith the w dependent
lambda matrix against the Lanczos implementation.
"""

T = 0
N_RANDOM = 1000
POPULATION = 2
N_W = 100
W_MIN = 2800
W_MAX = 3400
GAMMA = 25 / CC.Units.RY_TO_CM
NQIRR = 1
DYNMAT_FILE = "../ensemble_data_test/dyn"
DATA_DIR = "../ensemble_data_test"
N_LANCZOS_STEPS = 40

# Load the dynamical file
print(INFO)
dyn = CC.Phonons.Phonons(DYNMAT_FILE, NQIRR)

# Load the ensemble
ens = sscha.Ensemble.Ensemble(dyn, T, dyn.GetSupercell())
ens.load(DATA_DIR, POPULATION, N_RANDOM)

# Prepare the w array
w_array = np.linspace(W_MIN, W_MAX, N_W)
w_ry = w_array / CC.Units.RY_TO_CM

# Ok, now get the Raman with the lanczos
lanczos = sscha.DynamicalLanczos.Lanczos(ens)
lanczos.ignore_v3 = False
lanczos.ignore_v4 = True
lanczos.prepare_symmetrization()

# Get the raman vector for the green function
pol_x = np.array([1,0,0], dtype = np.double)
raman_v = dyn.GetRamanVector(pol_x, pol_x)
# Setup the raman vector as first vector for the Lanczos
lanczos.prepare_perturbation(raman_v, masses_on_multiply = False)
# Run the lanczos algorithm
lanczos.run(N_LANCZOS_STEP)

# Get the green function using the continued fraction
gf = lanczos.get_green_function_continued_fraction(w_ry,
                                                   use_terminator = True,
                                                   last_average = 5,
                                                   smearing = 0)
# Get the absorb
raman_lanczos = -np.imag(gf)

# Get the lanczos spectral function
spectral_lanczos = lanczos.get_spectral_function_from_Lenmann(w_ry, GAMMA)


# Now we can do the same calculation using the real dynamical
# Unwrap the ensemble to use symmetries
ens.unwrap_symmetries()
self_energies = ens.get_odd_correction(include_v4 = False,
                                       frequencies = w_ry,
                                       smearing = GAMMA,
                                       return_only_correction = True)
raman_standard = sscha.Dynamical.GetRamanResponce(dyn.GetSupercell(),
                                                  self_energies,
                                                  w_ry,
                                                  pol_x,
                                                  pol_x)

# Get also the spectral function from the standard
spectral_standard = sscha.Dynamical.get_spectral_function(dyn, dyn.GetSupercell(), self_energies, w_ry)

# Prepare the file to be saved
np.savetxt("raman_comparison.dat", np.array([w_array, raman_lanczos, raman_standard]).T, header= "Frequency [cm-1], Raman (Lanczos), Raman (Standard)")
np.savetxt("spectral_comparison.dat", np.array([w_array, spectral_lanczos, spectral_standard]).T, "Frequency [cm-1], Spectral function (Lanczos), Spectral function (Standard)")
print("Done.")
