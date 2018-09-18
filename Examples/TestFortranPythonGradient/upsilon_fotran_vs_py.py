"""
Test the get_upsilon_matrix
Fortran Vs Python
"""

import numpy as np

import cellconstructor as CC
import cellconstructor.Phonons

import sscha, sscha.Ensemble
import SCHAModules

DYN="../ensemble_data_test/dyn"
N_RAND=500
POPULATION=2
T = 0
DATA_DIR = "../ensemble_data_test"


# Load the dynamical matrix
dyn = CC.Phonons.Phonons(DYN)

# Load the ensemble
print "Loading the ensemble..."
ens = sscha.Ensemble.Ensemble(dyn, T)
ens.load(DATA_DIR, POPULATION, N_RAND)
print "Ensemble loaded."


# Get the upsilon matrix using the python subroutine
ups = dyn.GetUpsilonMatrix(T)

# Save it to a file
np.savetxt("py_ups.dat", ups, header = "Upsilon matrix generated with python")


# Use the fortran subroutine to compute the upsilon matrix
w, pols = dyn.DyagDinQ(0)
w /= 2
trans = CC.Methods.get_translations(pols, dyn.structure.get_masses_array()).astype(np.intc) # Convert bool to int
mass = np.array(dyn.structure.masses.values())
mass *= 2

ityp = dyn.structure.get_ityp() + 1 #+1 needed for python -> fortran array indexing

ups_f = SCHAModules.get_upsilon_matrix(w, pols, trans, mass, ityp, T)

# Save it to a file
np.savetxt("f_usp.dat", ups_f, header ="Upsilon matrix generated with fortran")

# Compute the distance between the two matrix
diff = ups - ups_f
print "Distance between them is:", np.einsum("ab,ba", diff, diff)


